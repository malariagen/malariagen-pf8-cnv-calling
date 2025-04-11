process CheckPloidy {
    // Look within each sample subdirectory within `meta.contig_ploidy_calls` and read each sample's
    // contig_ploidy.tsv file. Collect the name of the sample and the contig where the contig ploidy
    // was not equal to 1, as the HMM outputs for these cases should not be trusted. 

    input:
        tuple val(cluster),
              val(meta)
    
    output:
        path "combined_contig_ploidy.tsv"
    
    script:
        """
        list_of_sample_directories=\$(find "${meta.contig_ploidy_calls}/" -type d -name "SAMPLE_*")
        echo -e "SAMPLE\\tCONTIG\\tPLOIDY\\tPLOIDY_GQ" > combined_contig_ploidy.tsv

        for sample_dir in \$list_of_sample_directories; do
            contig_ploidy_file=\$(find "\$sample_dir" -name "contig_ploidy.tsv")
            sample_name=\$(awk '/@RG/ {print \$3}' "\$contig_ploidy_file" | sed 's/SM://')
            awk -v sample="\$sample_name" 'NR>2 {print sample"\t"\$0}' "\$contig_ploidy_file" >> combined_contig_ploidy.tsv
        done
        """
}