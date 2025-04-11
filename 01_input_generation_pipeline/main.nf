#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    
    // ===== CREATE CHANNELS FROM INPUT FILES ===== //
    // make channel containing fasta, fasta.fai, dict reference genome
    Channel.fromPath([ params.reference,
                    "${params.reference}.fai",
                    "${params.reference}".replaceAll(/fasta/, 'dict')])
        | collect
        | map {it -> tuple( it[0], it[1], it[2] )}
        | set { reference_genome_ch }

    // make channels for path to BAM files (i.e., manifest) and user-defined intervals of interest
    intervals_of_interest_ch = Channel.fromPath(params.intervals_of_interest)
    manifest_ch              = Channel.fromPath(params.manifest)

    // ===== FEED DATA FROM CHANNELS TO PROCESSES ===== //
    // First, segment the entire reference genome into intervals,
    reference_genome_ch
        | GenerateWholeGenomeIntervalList
        | set { interval_list_of_genome_ch }

    // then restrict the segmented genome to only include user-defined regions of interest,
    intervals_of_interest_ch
        | concat( interval_list_of_genome_ch ) // TODO: check shape here to remove collect?
        | collect
        | map {it -> tuple( it[0], it[1] )} // just feeding both channels as a tuple
        | RestrictIntervalList
        | set { restricted_list_of_intervals_ch }

    // and finally use the restricted list of intervals to count the number of reads per interval for each sample.
    manifest_ch
        | splitCsv(sep: "\t")
        | map { row -> tuple(row[0], row[1]) } // supply sample_id+path_to_bam pairs to process for parallel computing
        | combine( restricted_list_of_intervals_ch )
        | GenerateReadCountFile
        | map { it -> it.getName().tokenize(".")[0] + "\t${params.results_dir}/01_execution/readcounts/${it.getName()}"} // get pairs of sample_id+path_to_readcounts
        | collectFile( name: "${params.results_dir}/02_paths_to_readcounts.tsv", newLine: true, sort: true) // and write them to 02_paths_to_readcounts.tsv
    
    // Also annotate this intervals list of interest with GC-content in preparation for step 03 pipeline.
    restricted_list_of_intervals_ch
        | concat( reference_genome_ch )
        | collect
        | map {it -> tuple( it[0], it[1], it[2], it[3] )}
        | AnnotateIntervalList // concat the three reference channels to the restricted_list_of_intervals_ch and feed into process


}

process GenerateWholeGenomeIntervalList {
    // Takes the reference genome and segments it into intervals (or "bins") of equal length
    input:
        tuple path(reference_genome_fasta),
              path(reference_genome_fai),
              path(reference_genome_dict)

    output:
        path "whole_genome.interval_list"

    script:
        def java_Xmx = task.memory.mega-200

        """
        gatk PreprocessIntervals \
            --java-options "-Xmx${java_Xmx}m" \
            --reference             ${reference_genome_fasta} \
            --bin-length            ${params.bin_length} \
            --padding               0 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output                whole_genome.interval_list
        """
}

process RestrictIntervalList {
    // Takes the whole-genome list of intervals and restricts it
    // to the user-defined list of intervals of interest

    publishDir "${params.results_dir}", mode: "copy", overwrite: true, failOnError: true

    input:
        tuple path(intervals_of_interest),
              path(interval_list_of_genome) 

    output:
        path "03_restricted.interval_list"

    script:
        def java_Xmx = task.memory.mega-200

        """
        gatk IntervalListTools \
            --java-options "-Xmx${java_Xmx}m" \
            --INPUT        ${interval_list_of_genome} \
            --SECOND_INPUT ${intervals_of_interest} \
            --ACTION       OVERLAPS \
            --OUTPUT       03_restricted.interval_list.TEMP

        mv 03_restricted.interval_list.TEMP 03_restricted.interval_list
        """
}


process AnnotateIntervalList {
    // Annotates restricted.interval_list with GC content which allows
    // the model to avoid training on AT-rich regions of the genome
    // in 03_coverage_based_pipeline

    publishDir "${params.results_dir}", mode: "copy", overwrite: true, failOnError: true
    
    input:
        tuple path(restricted_list_of_intervals),
              path(reference_genome_fasta),
              path(reference_genome_fai),
              path(reference_genome_dict)

    output:
        path "03_annotated.interval_list"

    script:
        def java_Xmx = task.memory.mega-200

        """
        gatk AnnotateIntervals \
            --java-options "-Xmx${java_Xmx}m" \
            --reference             ${reference_genome_fasta} \
            --intervals             ${restricted_list_of_intervals} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output                03_annotated.interval_list.TEMP
        
        mv 03_annotated.interval_list.TEMP 03_annotated.interval_list # add comments
        """
}


process GenerateReadCountFile {
    // Takes the BAM files for each sample and calculates read counts
    // within the regions specified by the restricted_list_of_intervals
    
    publishDir "${params.results_dir}/01_execution/readcounts/", mode: "copy", overwrite: true, failOnError: true

    input: 
        tuple val (sample_id),
              path(path_to_bam),
              path(restricted_list_of_intervals)
    
    output:
        path "${sample_id}.counts.tsv"

    script:
        def java_Xmx = task.memory.mega-200
        
        """
        gatk CollectReadCounts \
            --java-options "-Xmx${java_Xmx}m" \
            --intervals             ${restricted_list_of_intervals} \
            --input                 ${path_to_bam} \
            --format                TSV \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output                ${sample_id}.counts.tsv.TEMP
        
        mv ${sample_id}.counts.tsv.TEMP ${sample_id}.counts.tsv
        """
}