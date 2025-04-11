process FitPredictContigPloidy {
    // Use training samples to fit the contig ploidy model, then make contig ploidy calls.
    // This is executed PER CLUSTER. 
    
    input:
        tuple val(cluster),
              val(meta)

    output:
        tuple val (cluster),
              val (meta),
              path("training-${cluster}-contig-ploidy-model"),
              path("training-${cluster}-contig-ploidy-calls")

    script:
        def input_flag_concatenated_readcount_paths = meta.partition_of_sample_metas.collect{ it.PATH_TO_READCOUNTS }.join(" --input ")
        def java_Xmx = task.memory.mega-200

        """
        mkdir -p tmp
        
        THEANO_FLAGS="base_compiledir=tmp" gatk DetermineGermlineContigPloidy \
            --java-options          "-Xmx${java_Xmx}m" \
            --tmp-dir               tmp \
            --input                 ${input_flag_concatenated_readcount_paths} \
            --intervals             ${meta.filtered_intervals} \
            --exclude-intervals     ${meta.blacklist_intervals} \
            --contig-ploidy-priors  ${meta.contig_ploidy_priors} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --verbosity             DEBUG \
            --output                . \
            --output-prefix         "training-${cluster}-contig-ploidy" \
            --mean-bias-standard-deviation ${params.contig_ploidy_model_mean_bias_standard_deviation} \
            --mapping-error-rate           ${params.contig_ploidy_model_mapping_error_rate} \
            --global-psi-scale             ${params.contig_ploidy_model_global_psi_scale} \
            --sample-psi-scale             ${params.contig_ploidy_model_sample_psi_scale}
        
        rm -rf tmp
        """
}