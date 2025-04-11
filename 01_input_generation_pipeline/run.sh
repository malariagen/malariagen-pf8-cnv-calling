#!/bin/bash
# Run using ./run.sh (you may need to chmod)

module load nextflow-23.10.0

bsub \
    -J cnv-01 \
    -o cnv-01.o \
    -e cnv-01.e \
    -n 1 \
    -M 4000 \
    -R 'select[mem>=4000] rusage[mem=4000] span[hosts=1]' \
    -q normal \
    "
    nextflow run main.nf \
        -profile                sanger_lsf \
        --results_dir           <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/ \
        --manifest              <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/01_paths_to_bams.tsv \
        --reference             <INSERT PATH HERE>/pf7/Pfalciparum.genome.fasta \
        --intervals_of_interest <INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/01_intervals_of_interest.interval_list \
    "