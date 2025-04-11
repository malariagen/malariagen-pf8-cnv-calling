process FilterIntervalsList {
    // Filter interval list using samples from training set to dispose of bad intervals.
    
    publishDir "${params.results_dir}", pattern: "03_filtered.interval_list", mode: "copy", overwrite: true, failOnError: true

    input:
        tuple val(partition_of_sample_metas),
              val(aux_files_meta)

    output:
        tuple path("03_filtered.interval_list"),
              val(aux_files_meta)

    script:
        // Concatenate paths to read counts with the --input flag so that each path has its own flag. 
        def input_flag_concatenated_readcount_paths = partition_of_sample_metas.collect{ it.PATH_TO_READCOUNTS }.join(" --input ")
        def java_Xmx = task.memory.mega-200

        """
        mkdir -p tmp

        gatk FilterIntervals \
            --java-options "-Xmx${java_Xmx}m" \
            --tmp-dir             tmp \
            --intervals           ${aux_files_meta.intervals} \
            --output              tmp/03_filtered.interval_list \
            --input               ${input_flag_concatenated_readcount_paths} \
            --annotated-intervals ${aux_files_meta.annotated_intervals} \
            --exclude-intervals   ${aux_files_meta.blacklist_intervals} \
            --minimum-gc-content  ${params.filter_intervals_minimum_gc_content} \
            --maximum-gc-content  ${params.filter_intervals_maximum_gc_content} \
            --minimum-mappability ${params.filter_intervals_minimum_mappability} \
            --maximum-mappability ${params.filter_intervals_maximum_mappability} \
            --interval-merging-rule                      OVERLAPPING_ONLY \
            --minimum-segmental-duplication-content      ${params.filter_intervals_minimum_segmental_duplication_content} \
            --maximum-segmental-duplication-content      ${params.filter_intervals_maximum_segmental_duplication_content} \
            --low-count-filter-count-threshold           ${params.filter_intervals_low_count_filter_count_threshold} \
            --low-count-filter-percentage-of-samples     ${params.filter_intervals_low_count_filter_percentage_of_samples} \
            --extreme-count-filter-minimum-percentile    ${params.filter_intervals_extreme_count_filter_minimum_percentile} \
            --extreme-count-filter-maximum-percentile    ${params.filter_intervals_extreme_count_filter_maximum_percentile} \
            --extreme-count-filter-percentage-of-samples ${params.filter_intervals_extreme_count_filter_percentage_of_samples}
        
        mv tmp/03_filtered.interval_list 03_filtered.interval_list
        """
}