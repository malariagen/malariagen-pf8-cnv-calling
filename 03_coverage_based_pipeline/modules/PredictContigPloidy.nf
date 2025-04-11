process PredictContigPloidy {
    // Use the previously trained contig ploidy model to make contig ploidy calls. 
    
    input:
        tuple val(cluster),
              val(meta)

    output:
        tuple val (cluster),
              val (meta),
              path("test-${cluster}-contig-ploidy-calls")

    script:
        def input_flag_concatenated_readcount_paths = meta.partition_of_sample_metas.collect{ it.PATH_TO_READCOUNTS }.join(" --input ")
        def java_Xmx = task.memory.mega-200

        """
        mkdir -p tmp
        
        THEANO_FLAGS="base_compiledir=tmp" gatk DetermineGermlineContigPloidy \
            --java-options       "-Xmx${java_Xmx}m" \
            --tmp-dir            tmp \
            --input              ${input_flag_concatenated_readcount_paths} \
            --model              ${meta.contig_ploidy_model} \
            --verbosity          DEBUG \
            --output             . \
            --output-prefix      "test-${cluster}-contig-ploidy" \
            --mapping-error-rate ${params.contig_ploidy_model_mapping_error_rate} \
            --sample-psi-scale   ${params.contig_ploidy_model_sample_psi_scale}
        
        rm -rf tmp
        """
}