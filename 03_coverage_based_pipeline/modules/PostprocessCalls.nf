process PostprocessCalls {
    // Takes the gCNV model, gCNV calls and contig ploidy calls to perform postprocessing,
    // which generates three output files: intervals, segments, and dCR (denoised copy ratio).

    // There is some strange-looking Python code down below. This code immediately parses the files generated
    // by the PostprocessGermlineCNVCalls, because we found Nextflow/LSF to prematurely publish incomplete files
    // to the results directory. This Python code ensures that the process fails and restarts when this occurs. 
    
    publishDir "${params.results_dir}/03_execution/postprocessing/${meta.sample_id}", mode: "move", overwrite: true, failOnError: true

    input:
        val meta
    
    output:
        path("${meta.sample_id}.*")

    script:
        def java_Xmx = task.memory.mega-200

        """
        mkdir -p tmp
        
        THEANO_FLAGS="base_compiledir=tmp" gatk PostprocessGermlineCNVCalls \
            --java-options                "-Xmx${java_Xmx}m" \
            --tmp-dir                     tmp \
            --calls-shard-path            ${meta.gcnv_calls} \
            --model-shard-path            ${meta.gcnv_model} \
            --contig-ploidy-calls         ${meta.contig_ploidy_calls} \
            --sample-index                ${meta.sample_index} \
            --autosomal-ref-copy-number   1 \
            --output-denoised-copy-ratios tmp/${meta.sample_id}.dCR.tsv \
            --output-genotyped-intervals  tmp/${meta.sample_id}.intervals.vcf.gz \
            --output-genotyped-segments   tmp/${meta.sample_id}.segments.vcf.gz
        
        python -c '
import pandas as pd
pd.read_csv("tmp/${meta.sample_id}.dCR.tsv", sep = "\\t", skiprows = 18)
pd.read_csv("tmp/${meta.sample_id}.intervals.vcf.gz", sep = "\\t", skiprows = 24)
pd.read_csv("tmp/${meta.sample_id}.segments.vcf.gz", sep = "\\t", skiprows = 35)
        '

        mv tmp/${meta.sample_id}.intervals.vcf.gz     ${meta.sample_id}.intervals.vcf.gz
        mv tmp/${meta.sample_id}.intervals.vcf.gz.tbi ${meta.sample_id}.intervals.vcf.gz.tbi
        mv tmp/${meta.sample_id}.segments.vcf.gz      ${meta.sample_id}.segments.vcf.gz
        mv tmp/${meta.sample_id}.segments.vcf.gz.tbi  ${meta.sample_id}.segments.vcf.gz.tbi
        mv tmp/${meta.sample_id}.dCR.tsv              ${meta.sample_id}.dCR.tsv
        
        rm -rf tmp
        """
}