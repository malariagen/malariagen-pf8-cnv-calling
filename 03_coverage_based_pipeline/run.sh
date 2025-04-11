#!/bin/bash
# Run using ./run.sh (you may need to chmod)

module load nextflow-23.10.0

bsub \
    -J cnv-03 \
    -o cnv-03.o \
    -e cnv-03.e \
    -n 1 \
    -M 2000 \
    -R 'select[mem>=2000] rusage[mem=2000] span[hosts=1]' \
    -q basement \
    "
    nextflow run main.nf \
        -resume \
        -profile                 sanger_lsf \
        --results_dir            <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8 \
        --manifest               <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/03_sample_cluster_assignment.tsv \
        --intervals              <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/03_restricted.interval_list \
        --annotated_intervals    <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/03_annotated.interval_list \
        --blacklist_intervals    <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/03_blacklist.interval_list \
        --contig_ploidy_priors   <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/03_contig_ploidy_priors.tsv
    "