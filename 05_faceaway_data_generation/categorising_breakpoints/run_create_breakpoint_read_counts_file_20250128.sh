
MANIFEST_FN=$1
JOB=$LSB_JOBINDEX

IN=$(sed "$JOB q;d" $MANIFEST_FN)
read -a LINE <<< "$IN"
SAMPLE=${LINE[0]}
BAM_FN=${LINE[1]}

echo $SAMPLE
echo $BAM_FN



python3 create_breakpoint_read_counts_file_20250128.py     --sample $SAMPLE     --bam_fn $BAM_FN     --output_dir samples     --tandem_dup_breakpoints_fn tandem_dup_breakpoints_20250128.txt     --dup_trpinv_dup_breakpoints_fn dup_trpinv_dup_breakpoints_20250128.txt
