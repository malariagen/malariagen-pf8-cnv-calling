#!/bin/bash
# Run using ./run.sh (you may need to chmod)

CHUNK_SIZE=100

PATH_TO_POSTPROCESSING="<INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/03_execution/postprocessing/"
TOTAL_SAMPLES=$(ls $PATH_TO_POSTPROCESSING | wc -l)
N_JOBS=$(expr \( $TOTAL_SAMPLES + $CHUNK_SIZE - 1 \) / $CHUNK_SIZE)

bsub \
    -J "diagnostics[1-$N_JOBS]" \
    -o generate_plots%J%I.output \
    -e generate_plots%J%I.error \
    -M 500 \
    -R "select[mem>500] rusage[mem=500] span[hosts=1]" \
    -q normal \
    "python3 generate_plots.py \${LSB_JOBINDEX} $CHUNK_SIZE"
