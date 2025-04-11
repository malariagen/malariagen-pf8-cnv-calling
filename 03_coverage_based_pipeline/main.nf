#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FilterIntervalsList }    from "./modules/FilterIntervalsList.nf"
include { FitPredictContigPloidy } from "./modules/FitPredictContigPloidy.nf"
include { FitPredictGermlineCNV }  from "./modules/FitPredictGermlineCNV.nf"
include { 
    PredictContigPloidy as PredictContigPloidy;
    PredictContigPloidy as PredictContigPloidyFromPretrainedModel
} from "./modules/PredictContigPloidy.nf"
include {
    PredictGermlineCNV as PredictGermlineCNV;
    PredictGermlineCNV as PredictGermlineCNVFromPretrainedModel
} from "./modules/PredictGermlineCNV.nf"

include {
    CollateCalls as CollateTestCalls;
    CollateCalls as CollateTrainCalls;
    CollateCalls as CollatePretrainedModelCalls;
} from "./modules/CollateCalls.nf"
include {
    PostprocessCalls as PostprocessTestCalls;
    PostprocessCalls as PostprocessTrainCalls;
    PostprocessCalls as PostprocessPretrainedModelCalls;
} from "./modules/PostprocessCalls.nf"
include {
    CheckPloidy as CheckTestPloidy;
    CheckPloidy as CheckTrainPloidy;
    CheckPloidy as CheckPretrainedModelPloidy;
} from "./modules/CheckPloidy.nf"


workflow {
    
    // ===== STEP 1.1: CREATE CHANNELS USING AUXILIARY FILES ===== //

    Channel.fromPath( params.manifest ) // From samples manifest file.
        | splitCsv( sep: "\t", header: true )
        | set { sample_meta_ch }
        // One row/dict/meta per sample:
        //  - SAMPLE_ID             : Sample ID (str)
        //  - PATH_TO_READCOUNTS    : Path to readcounts file (path)
        //  - IN_TRAINING_BLACKLIST : Blacklisted from training set? True/False
        //  - IN_TRAINING_SET       : Used for model training? True/False
        //  - CLUSTER               : Training cluster identity

    Channel.fromPath([ // Collect auxiliary files into single channel.
            params.intervals,
            params.annotated_intervals,
            params.blacklist_intervals,
            params.contig_ploidy_priors
        ])
        | collect
        | map { it ->
            def meta = [
                intervals            : it[0],
                annotated_intervals  : it[1],
                blacklist_intervals  : it[2],
                contig_ploidy_priors : it[3]
            ]

            [ meta ]
        }
        | set { aux_files_meta_ch } // [ meta ]




    // ===== STEP 1.2: CREATE CHANNEL FOR TRAINING SAMPLES ===== //

    // FILTER AWAY SAMPLES THAT ALREADY HAVE A PRETRAINED MODEL AND ARE NOT IN THE TRAINING SET. 
    sample_meta_ch                        // [ sample_meta ] for each sample
        | filter { it.PATH_TO_PRETRAINED_CONTIG_PLOIDY_MODEL == "." }
        | filter { it.PATH_TO_PRETRAINED_GCNV_MODEL == "." }
        | filter { it.IN_TRAINING_SET == "True"}
        | set { training_sample_meta_ch } // [ sample_meta ] for each sample

    // THEN, SPLIT THE SAMPLES UP BASED ON THE `CLUSTER` KEY. 
    training_sample_meta_ch                           // [ sample_meta ] for each sample
        | map { it -> tuple(it.CLUSTER, it) }         // [ cluster, sample_meta ] for each sample
        | groupTuple
        | set { training_sample_meta_per_cluster_ch } // [ cluster, partition_of_sample_metas ] for each cluster
    // NUMBER OF SAMPLES IN EACH PARTITION_OF_SAMPLE_METAS IS DICTATED IN `02_sample_clustering`, BY DEFAULT 500. 




    // ===== STEP 1.3: CREATE CHANNEL FOR TEST SAMPLES ===== //

    // FILTER AWAY SAMPLES THAT ALREADY HAVE A PRETRAINED MODEL AND ARE IN THE TRAINING SET. 
    sample_meta_ch                    // [ sample_meta ] for each sample
        | filter { it.PATH_TO_PRETRAINED_CONTIG_PLOIDY_MODEL == "." }
        | filter { it.PATH_TO_PRETRAINED_GCNV_MODEL == "." }
        | filter { it.IN_TRAINING_SET == "False"}
        | set { test_sample_meta_ch } // [ sample_meta ] for each sample

    // THEN, SPLIT THE SAMPLES UP BASED ON THE `CLUSTER` KEY. 
    test_sample_meta_ch                           // [ sample_meta ] for each sample
        | map { it -> tuple(it.CLUSTER, it) }     // [ cluster, sample_meta ] for each sample
        | groupTuple ( size : params.test_sample_max_partition_size , remainder : true ) // Partitioning test samples for parallellisation
        | set { test_sample_meta_per_cluster_ch } // [ cluster, partition_of_sample_metas ] for each partition for each cluster
    // NUMBER OF SAMPLES IN EACH PARTITION_OF_SAMPLE_METAS IS DICTATED BY `params.test_sample_max_partition_size`, BY DEFAULT 200. 




    // ===== STEP 1.4: CREATE CHANNEL FOR SAMPLES WITH A PRETRAINED MODEL ===== //

    // FILTER AWAY SAMPLES WITHOUT A PRETRAINED MODEL. 
    sample_meta_ch                                // [ sample_meta ] for each sample
        | filter { it.PATH_TO_PRETRAINED_CONTIG_PLOIDY_MODEL != "." }
        | filter { it.PATH_TO_PRETRAINED_GCNV_MODEL != "." }
        | set { pretrained_model_sample_meta_ch } // [ sample_meta ] for each sample

    // THEN, SPLIT THE SAMPLES UP BASED ON THE `CLUSTER` KEY. 
    pretrained_model_sample_meta_ch                           // [ sample_meta ] for each sample
        | map { it -> tuple(it.CLUSTER, it) }                 // [ cluster, sample_meta ] for each sample
        | groupTuple ( size : params.pretrained_model_sample_max_partition_size, remainder : true ) // Partitioning samples with pretrained models for parallellisation
        | set { pretrained_model_sample_meta_per_cluster_ch } // [ cluster, partition_of_sample_metas ] for each partition for each cluster
    // NUMBER OF SAMPLES IN EACH PARTITION_OF_SAMPLE_METAS IS DICTATED BY `params.pretrained_model_sample_max_partition_size`, BY DEFAULT 10. 




    // ===== STEP 1.5: CREATE FILTERED INTERVALS LIST AND ADD TO AUXILIARY FILES CHANNEL ===== //

    // GATHER ALL TRAINING SAMPLES AND USE THEM FOR FILTERING OUT BAD INTERVALS.
    training_sample_meta_ch                                   // [ sample_meta ] for each sample
        | collect                                             // [ partition_of_sample_metas ]
        | merge ( aux_files_meta_ch ) { a, b -> tuple(a, b) } // [ partition_of_sample_metas, aux_files_meta ]
        | FilterIntervalsList                                 // [ path("03_filtered.interval_list"), aux_files_meta ]
        | map { it ->
            def meta = [
                filtered_intervals   : it[0],
                intervals            : it[1].intervals,
                annotated_intervals  : it[1].annotated_intervals,
                blacklist_intervals  : it[1].blacklist_intervals,
                contig_ploidy_priors : it[1].contig_ploidy_priors
                
            ] // Added filtered_intervals file to the `aux_files_meta_ch`. 

            [ meta ]
        }
        | set { aux_files_meta_ch } // [ aux_files_meta ]




    // ===== STEP 2.1 FIT CONTIG PLOIDY MODELS USING TRAINING SAMPLES AND MAKE PLOIDY CALLS ===== //

    // Context: Predict the ploidy of each contig in each sample, i.e., the number of copies/chromosomes in the contig. 

    // FEED EACH TRAINING CLUSTER INTO MODEL TRAINING ALONG WITH AUXILIARY FILES.
    training_sample_meta_per_cluster_ch // [ cluster, partition_of_sample_metas ] for each cluster
        | combine ( aux_files_meta_ch ) // [ cluster, partition_of_sample_metas, aux_files_meta ] for each cluster
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1],
                filtered_intervals       : it[2].filtered_intervals,
                annotated_intervals      : it[2].annotated_intervals,
                blacklist_intervals      : it[2].blacklist_intervals,
                contig_ploidy_priors     : it[2].contig_ploidy_priors
            ]

            [ cluster, meta ]
        }
        // FIT THE CONTIG PLOIDY MODEL FOR EACH CLUSTER AND ALSO MAKE CONTIG PLOIDY CALLS.
        | FitPredictContigPloidy // [ cluster, meta, model, calls ] for each cluster
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[2],
                contig_ploidy_calls      : it[3]
            ]

            [ cluster, meta ] // Added contig ploidy model and calls to meta.
        }
        | set { training_samples_contig_ploidy_ch } // [ cluster, meta ] for each cluster
    
    // PUBLISH THE TRAINED MODELS TO THE RESULTS DIRECTORY. 
    training_samples_contig_ploidy_ch
        | map { it ->
            def tempDir = file(it[1].contig_ploidy_model)
            tempDir.copyTo( params.results_dir + "/models/" + tempDir.getBaseName() )
        }




    // ===== STEP 2.2 FIT GCNV MODELS USING TRAINING SAMPLES AND MAKE GCNV CALLS ===== //

    // Context: Use the ploidy prediction to normalise the coverage and get a copy ratio. 

    // FEED EACH TRAINING CLUSTER INTO MODEL TRAINING ALONG WITH AUXILIARY FILES, PLOIDY CALLS AND PLOIDY MODEL.
    training_samples_contig_ploidy_ch // [ cluster, meta ] for each cluster
        | FitPredictGermlineCNV       // [ cluster, meta, model, calls ] for each cluster
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[1].contig_ploidy_model,
                contig_ploidy_calls      : it[1].contig_ploidy_calls,
                gcnv_model               : it[2],
                gcnv_calls               : it[3]
            ]

            [ cluster, meta ] // Added gCNV model and calls to meta.
        }
        | set { training_samples_gcnv_ch } // [ cluster, meta ]

    // PUBLISH THE TRAINED MODELS TO THE RESULTS DIRECTORY. 
    training_samples_gcnv_ch
        | map { it ->
            def tempDir = file(it[1].gcnv_model)
            tempDir.copyTo( params.results_dir + "/models/" + tempDir.getBaseName() )
        }




    // ===== STEP 2.3: FEED INFERENCES FROM TRAINING SAMPLES INTO POSTPROCESSING ===== //

    // Context: Viterbi algorithm to infer the most likely copy number by using ploidy calls, gCNV calls and gCNV model. 
    // Assigns Phred-scaled probabilities and generates intervals.vcf and segments.vcf. 

    // GATK `GermlineCNVCaller` CREATES A FOLDER FOR EACH PROCESSED SAMPLE TITLED SOMETHING LIKE "SAMPLE_25"
    // WITH A `sample_name.txt` INSIDE, WHICH CONTAINS THE ORIGINAL SAMPLE ID. `CollateTrainCalls` USES SOME
    // BASH SCRIPTING TO PARSE THE SAMPLE IDS OUT, TO HELP WITH PARALLELISING `PostprocessTrainCalls`. 

    training_samples_gcnv_ch                  // [ cluster, meta ] for each cluster
        | map { it -> it[1] }
        | CollateTrainCalls                   // [ path("sample_names_and_indexes.txt") ] for each cluster
        | splitCsv( sep: "\t", header: true ) // { cluster, sample_name, sample_index } for each sample
        | map { it ->
            def meta = [
                sample_id           : it.SAMPLE_ID,
                sample_index        : it.SAMPLE_INDEX,
                gcnv_calls          : it.GCNV_CALLS,
                gcnv_model          : it.GCNV_MODEL,
                contig_ploidy_calls : it.CONTIG_PLOIDY_CALLS
            ]

            [ meta ]
        }
        | set { training_samples_for_postprocessing_ch }

    training_samples_for_postprocessing_ch
        | PostprocessTrainCalls




    // ===== STEP 3.1: USE TRAINED CONTIG PLOIDY MODELS FROM STEP 2.1 TO MAKE PLOIDY CALLS FOR TEST SAMPLES ===== //

    // Context: Predict the ploidy of each contig in each sample, i.e., the number of copies/chromosomes in the contig. 

    // MATCH UP THE CONTIG PLOIDY MODEL WITH THE CORRESPONDING CLUSTER OF TEST SAMPLES AND MAKE CONTIG PLOIDY CALLS.
    test_sample_meta_per_cluster_ch                            // [ cluster, partition_of_sample_metas ] for each partition per cluster
        | combine ( training_samples_contig_ploidy_ch, by: 0 ) // [ cluster, partition_of_sample_metas, meta ]
        | map { it ->
            def cluster = it[0]
            def meta = [
                partition_of_sample_metas : it[1],
                contig_ploidy_model       : it[2].contig_ploidy_model,
            ]

            [ cluster, meta ]
        }
        | PredictContigPloidy // [ cluster, meta, contig_ploidy_calls ]
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[1].contig_ploidy_model,
                contig_ploidy_calls      : it[2]
            ]

            [ cluster, meta ] // Added contig ploidy calls to meta.
        }
        | set { test_samples_contig_ploidy_ch } // [ cluster, meta ]




    // ===== STEP 3.2 USE TRAINED GCNV MODELS FROM STEP 2.2 AND MAKE GCNV CALLS FOR TEST SAMPLES ===== //

    // Context: Use the ploidy prediction to normalise the coverage and get a copy ratio. 

    // MATCH UP THE CONTIG PLOIDY CALLS WITH THE CORRESPONDING CLUSTER OF TEST SAMPLES AND
    // THE CORRESPONDING GCNV MODEL TO MAKE GCNV CALLS
    test_samples_contig_ploidy_ch                     // [ cluster, meta ]
        | combine ( training_samples_gcnv_ch, by: 0 ) // [ cluster, meta, model_meta ]
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[1].contig_ploidy_model,
                contig_ploidy_calls      : it[1].contig_ploidy_calls,
                gcnv_model               : it[2].gcnv_model
            ]

            [ cluster, meta ] // Added gCNV model to meta.
        }
        | PredictGermlineCNV // [ cluster, meta, gcnv_calls ]
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[1].contig_ploidy_model,
                contig_ploidy_calls      : it[1].contig_ploidy_calls,
                gcnv_model               : it[1].gcnv_model,
                gcnv_calls               : it[2]
            ]

            [ cluster, meta ] // Added gCNV calls to meta.
        }
        | set { test_samples_gcnv_ch } // [ cluster, meta ]




    // ===== STEP 3.3: FEED INFERENCES FROM TEST SAMPLES INTO POSTPROCESSING ===== //

    // Context: Viterbi algorithm to infer the most likely copy number by using ploidy calls, gCNV calls and gCNV model. 
    // Assigns Phred-scaled probabilities and generates intervals.vcf and segments.vcf. 

    // SIMILAR TO STEP 2.3, BUT FOR TEST SAMPLES. 
    test_samples_gcnv_ch                      // [ cluster, meta ]
        | map { it -> it[1] }
        | CollateTestCalls                    // [ path("sample_names_and_indexes.txt") ]
        | splitCsv( sep: "\t", header: true ) // { cluster, sample_name, sample_index } for each sample
        | map { it ->
            def meta = [
                sample_id           : it.SAMPLE_ID,
                sample_index        : it.SAMPLE_INDEX,
                gcnv_calls          : it.GCNV_CALLS,
                gcnv_model          : it.GCNV_MODEL,
                contig_ploidy_calls : it.CONTIG_PLOIDY_CALLS
            ]

            [ meta ]
        }
        | set { test_samples_for_postprocessing_ch }

    test_samples_for_postprocessing_ch
        | PostprocessTestCalls




    // ===== STEP 4.1: USE PRETRAINED CONTIG PLOIDY MODELS FROM MANIFEST TO MAKE PLOIDY CALLS FOR SAMPLES WITH PRETRAINED MODELS PROVIDED ===== //

    // Context: Predict the ploidy of each contig in each sample, i.e., the number of copies/chromosomes in the contig. 

    pretrained_model_sample_meta_per_cluster_ch // [ cluster, partition_of_sample_metas ] for each partition per cluster
        | map { it ->
            def cluster = it[0]
            def meta = [
                partition_of_sample_metas : it[1],
                contig_ploidy_model       : it[1][0].PATH_TO_PRETRAINED_CONTIG_PLOIDY_MODEL,
            ]

            [ cluster, meta ]
        }
        | PredictContigPloidyFromPretrainedModel
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[1].contig_ploidy_model,
                contig_ploidy_calls      : it[2]
            ]

            [ cluster, meta ] // Added contig ploidy calls to meta.
        }
        | set { pretrained_model_samples_contig_ploidy_ch } // [ cluster, meta ]




    // ===== STEP 4.2 USE PRETRAINED GCNV MODELS FROM MANIFEST AND MAKE GCNV CALLS FOR SAMPLE WITH PRETRAINED MODELS PROVIDED ===== //

    // Context: Use the ploidy prediction to normalise the coverage and get a copy ratio. 

    pretrained_model_samples_contig_ploidy_ch // [ cluster, meta ]
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[1].contig_ploidy_model,
                contig_ploidy_calls      : it[1].contig_ploidy_calls,
                gcnv_model               : it[1].partition_of_sample_metas[0].PATH_TO_PRETRAINED_GCNV_MODEL,
            ]

            [ cluster, meta ] // Added gCNV model to meta.
        }
        | PredictGermlineCNVFromPretrainedModel // [ cluster, meta, gcnv_calls ]
        | map { it ->
            def cluster = it[0]
            def meta    = [
                partition_of_sample_metas: it[1].partition_of_sample_metas,
                filtered_intervals       : it[1].filtered_intervals,
                annotated_intervals      : it[1].annotated_intervals,
                blacklist_intervals      : it[1].blacklist_intervals,
                contig_ploidy_priors     : it[1].contig_ploidy_priors,
                contig_ploidy_model      : it[1].contig_ploidy_model,
                contig_ploidy_calls      : it[1].contig_ploidy_calls,
                gcnv_model               : it[1].gcnv_model,
                gcnv_calls               : it[2]
            ]

            [ cluster, meta ] // Added gCNV calls to meta.
        }
        | set { pretrained_model_samples_gcnv_ch } // [ cluster, meta ]




    // ===== STEP 4.3: FEED INFERENCES FROM SAMPLES WITH PRETRAINED MODELS INTO POSTPROCESSING ===== //

    // Context: Viterbi algorithm to infer the most likely copy number by using ploidy calls, gCNV calls and gCNV model. 
    // Assigns Phred-scaled probabilities and generates intervals.vcf and segments.vcf. 

    // SIMILAR TO STEP 2.3 and 3.3, BUT FOR SAMPLES WITH PRETRAINED MODELS. 
    pretrained_model_samples_gcnv_ch          // [ cluster, meta ]
        | map { it -> it[1] }
        | CollatePretrainedModelCalls         // [ path("sample_names_and_indexes.txt") ]
        | splitCsv( sep: "\t", header: true ) // { cluster, sample_name, sample_index } for each sample
        | map { it ->
            def meta = [
                sample_id           : it.SAMPLE_ID,
                sample_index        : it.SAMPLE_INDEX,
                gcnv_calls          : it.GCNV_CALLS,
                gcnv_model          : it.GCNV_MODEL,
                contig_ploidy_calls : it.CONTIG_PLOIDY_CALLS
            ]

            [ meta ]
        }
        | set { pretrained_model_samples_for_postprocessing_ch }

    pretrained_model_samples_for_postprocessing_ch
        | PostprocessPretrainedModelCalls    




    // ===== STEP 5.1: CHECK CONTIG PLOIDY OF TRAINING SAMPLES - CREATE `failed_training_ploidy.tsv` ===== //

    // READ THE CONTIG PLOIDY FILE AND FLAG SAMPLES AND CONTIGS WHOSE CONTIG PLOIDY IS NOT EQUAL TO 1,
    // THEN COLLATE THEM INTO A FILE. NOTE, THIS MEANS THAT SAMPLES WERE USED FOR TRAINING EVEN IF THEY
    // ENDED UP FAILING CONTIG PLOIDY. 
    training_samples_contig_ploidy_ch
        | CheckTrainPloidy
        | splitCsv( sep: "\t", header: true )
        | filter { it.PLOIDY.toString() != "1" } // Collect contigs where the ploidy was not 1
        | map { it ->
            "${it.SAMPLE}\t${it.CONTIG}\t${it.PLOIDY}\t${it.PLOIDY_GQ}"
        }
        | collectFile( name: "${params.results_dir}/03_execution/failed_training_ploidy.tsv", newLine: true, sort: true )




    // ===== STEP 5.2: CHECK CONTIG PLOIDY OF TEST SAMPLES - CREATE `failed_test_ploidy.tsv` ===== //

    // READ THE CONTIG PLOIDY FILE AND FLAG SAMPLES AND CONTIGS WHOSE CONTIG PLOIDY IS NOT EQUAL TO 1,
    // THEN COLLATE THEM INTO A FILE. 
    test_samples_contig_ploidy_ch
        | CheckTestPloidy
        | splitCsv( sep: "\t", header: true )
        | filter { it.PLOIDY.toString() != "1" } // Collect contigs where the ploidy was not 1
        | map { it ->
            "${it.SAMPLE}\t${it.CONTIG}\t${it.PLOIDY}\t${it.PLOIDY_GQ}"
        }
        | collectFile( name: "${params.results_dir}/03_execution/failed_test_ploidy.tsv", newLine: true, sort: true )




    // ===== STEP 5.3: CHECK CONTIG PLOIDY OF SAMPLES WITH PRETRAINED MODELS - CREATE `failed_pretrained_model_ploidy.tsv` ===== //

    // READ THE CONTIG PLOIDY FILE AND FLAG SAMPLES AND CONTIGS WHOSE CONTIG PLOIDY IS NOT EQUAL TO 1,
    // THEN COLLATE THEM INTO A FILE. 
    pretrained_model_samples_contig_ploidy_ch
        | CheckPretrainedModelPloidy
        | splitCsv( sep: "\t", header: true )
        | filter { it.PLOIDY.toString() != "1" } // Collect contigs where the ploidy was not 1
        | map { it ->
            "${it.SAMPLE}\t${it.CONTIG}\t${it.PLOIDY}\t${it.PLOIDY_GQ}"
        }
        | collectFile( name: "${params.results_dir}/03_execution/failed_pretrained_model_ploidy.tsv", newLine: true, sort: true )




}