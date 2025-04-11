process FitPredictGermlineCNV {
    // Use contig ploidy calls to train the gCNV model, then make gCNV calls.
    
    input:
        tuple val(cluster),
              val(meta)
    
    output:
        tuple val (cluster),
              val (meta),
              path("training-${cluster}-gcnv-model"),
              path("training-${cluster}-gcnv-calls")

    script:
        def input_flag_concatenated_readcount_paths = meta.partition_of_sample_metas.collect{ it.PATH_TO_READCOUNTS }.join(" --input ")
        def java_Xmx = task.memory.mega-200

        """
        mkdir -p tmp
        
        THEANO_FLAGS="base_compiledir=tmp" gatk GermlineCNVCaller \
            --java-options          "-Xmx${java_Xmx}m" \
            --tmp-dir               tmp \
            --run-mode              COHORT \
            --contig-ploidy-calls   ${meta.contig_ploidy_calls} \
            --intervals             ${meta.filtered_intervals} \
            --annotated-intervals   ${meta.annotated_intervals} \
            --exclude-intervals     ${meta.blacklist_intervals} \
            --input                 ${input_flag_concatenated_readcount_paths} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --verbosity             DEBUG \
            --output                . \
            --output-prefix         "training-${cluster}-gcnv" \
            --p-alt                                  ${params.gcnv_model_p_alt} \
            --p-active                               ${params.gcnv_model_p_active}\
            --cnv-coherence-length                   ${params.gcnv_model_cnv_coherence_length} \
            --class-coherence-length                 ${params.gcnv_model_class_coherence_length} \
            --max-copy-number                        ${params.gcnv_model_max_copy_number} \
            --max-bias-factors                       ${params.gcnv_model_max_bias_factors} \
            --mapping-error-rate                     ${params.gcnv_model_mapping_error_rate} \
            --interval-psi-scale                     ${params.gcnv_model_interval_psi_scale} \
            --sample-psi-scale                       ${params.gcnv_model_sample_psi_scale} \
            --depth-correction-tau                   ${params.gcnv_model_depth_correction_tau} \
            --log-mean-bias-standard-deviation       ${params.gcnv_model_log_mean_bias_standard_deviation} \
            --init-ard-rel-unexplained-variance      ${params.gcnv_model_init_ard_rel_unexplained_variance} \
            --num-gc-bins                            ${params.gcnv_model_num_gc_bins} \
            --gc-curve-standard-deviation            ${params.gcnv_model_gc_curve_standard_deviation} \
            --copy-number-posterior-expectation-mode ${params.gcnv_model_copy_number_posterior_expectation_mode} \
            --enable-bias-factors                    ${params.gcnv_model_enable_bias_factors} \
            --active-class-padding-hybrid-mode       ${params.gcnv_model_active_class_padding_hybrid_mode} \
            --learning-rate                          ${params.gcnv_inference_learning_rate} \
            --adamax-beta-1                          ${params.gcnv_inference_adamax_beta_1} \
            --adamax-beta-2                          ${params.gcnv_inference_adamax_beta_2} \
            --log-emission-samples-per-round         ${params.gcnv_inference_log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ${params.gcnv_inference_log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds           ${params.gcnv_inference_log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch              ${params.gcnv_inference_max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs        ${params.gcnv_inference_max_advi_iter_subsequent_epochs} \
            --min-training-epochs                    ${params.gcnv_inference_min_training_epochs} \
            --max-training-epochs                    ${params.gcnv_inference_max_training_epochs} \
            --initial-temperature                    ${params.gcnv_inference_initial_temperature} \
            --num-thermal-advi-iters                 ${params.gcnv_inference_num_thermal_advi_iters} \
            --convergence-snr-averaging-window       ${params.gcnv_inference_convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold      ${params.gcnv_inference_convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window       ${params.gcnv_inference_convergence_snr_countdown_window} \
            --max-calling-iters                      ${params.gcnv_inference_max_calling_iters} \
            --caller-update-convergence-threshold    ${params.gcnv_inference_caller_update_convergence_threshold} \
            --caller-internal-admixing-rate          ${params.gcnv_inference_caller_internal_admixing_rate} \
            --caller-external-admixing-rate          ${params.gcnv_inference_caller_external_admixing_rate} \
            --disable-annealing                      ${params.gcnv_inference_disable_annealing}
        
        rm -rf tmp
        """
}