process CollateCalls {
    // Helper bash function for parsing GATK's nested output file structure. GCNV caller
    // creates a folder for each processed sample title "SAMPLE_25" with a `sample_name.txt`
    // inside, which contains the original sample ID. This function creates a tab-separated
    // .txt file containing the cluster identity, sample ID and sample index (e.g., "25"
    // if the sample was put into folder "SAMPLE_25"). 
    
    input:
        val(meta)
    
    output:
        path("sample_names_and_indexes.txt")
    
    script:
        """
        list_of_sample_directories=\$(find "${meta.gcnv_calls}/" -type d -name "SAMPLE_*")
        list_of_sample_indexes=(\$(echo "\$list_of_sample_directories" | sed 's|.*/SAMPLE_||'))

        sample_names=()
        for dir in \$list_of_sample_directories; do
            sample_name=\$(cat "\$dir/sample_name.txt")
            sample_names+=("\$sample_name")
        done

        output_file="sample_names_and_indexes.txt"

        echo -e "GCNV_CALLS\tGCNV_MODEL\tCONTIG_PLOIDY_CALLS\tSAMPLE_ID\tSAMPLE_INDEX" > "\$output_file"

        for ((i=0; i<\${#sample_names[@]}; i++)); do
            echo -e "${meta.gcnv_calls}\t${meta.gcnv_model}\t${meta.contig_ploidy_calls}\t\${sample_names[\$i]}\t\${list_of_sample_indexes[\$i]}" >> "\$output_file"
        done
        """
}