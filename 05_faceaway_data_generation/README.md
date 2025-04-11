# 05_faceaway_data_generation

This codebase was originally developed by Richard. Chiyun has since integrated it with the coverage-based CNV calling workflow.

---

## Overview

This subdirectory generates `breakpoint_read_counts_20250128.tsv`, a file processed later in `../06_final_calls_and_pf8_figures/2025_01_29_-_create_combined_amplification_calls.ipynb` to produce the final faceaway calls. Note that the faceaway calls themselves are not generated here. This step consists of:

1. **Searching for new breakpoints**: Using notebooks in `breakpoint_searching_notebooks/`.
2. **Calculating read counts for each breakpoint**: Using scripts in `categorising_breakpoints/`.

---

## Detailed Instructions

### 1. Search for New Breakpoints

We used the following notebooks to manually refine the list of known breakpoints for tandem duplications and duplication-triplication/inversion-duplications (DTDs):

- `breakpoint_searching_notebooks/2024_12_03_-_breakpoint_searching_crt.ipynb`: Identifies CRT breakpoints.
- `breakpoint_searching_notebooks/2025_01_27_-_breakpoint_searching_gch1_mdr1.ipynb`: Identifies breakpoints in GCH1 and MDR1.
- Output files:
  - `breakpoint_searching_notebooks/tandem_dup_breakpoints_20250128.txt`: Stores tandem duplication breakpoints.
  - `breakpoint_searching_notebooks/dup_trpinv_dup_breakpoints_20250128.txt`: Stores DTD breakpoints.

### 2. Calculate Faceaway Read Counts

All computations were performed in `categorising_breakpoints/2025_01_28_-_determine_breakpoint_read_evidence.ipynb`. This notebook:

- Develops and verifies Python functions, which are then manually incorporated into `categorising_breakpoints/create_breakpoint_read_counts_file_20250128.py`.
- Generates auxiliary files for submission of an LSF array job:
  - `categorising_breakpoints/run_create_breakpoint_read_counts_file_20250128.sh`
  - `categorising_breakpoints/pf8_qc_pass_bams_20250128.tsv` (list of BAM file paths).
- Submits an array job to process all samples, producing per-sample `.tsv` files in `categorising_breakpoints/all_samples/`.
- Concatenates results into:
  - `categorising_breakpoints/breakpoint_read_counts_noheader_20250128.tsv` (data without column headers).
  - `categorising_breakpoints/breakpoint_read_counts_header_20250128.tsv` (column headers).
  - `breakpoint_read_counts_20250128.tsv` (final read count proportions file).

---

## Directory Contents

### **Search for New Breakpoints**
- `breakpoint_searching_notebooks/`
  - `2024_12_03_-_breakpoint_searching_crt.ipynb`: Identifies CRT breakpoints.
  - `2025_01_27_-_breakpoint_searching_gch1_mdr1.ipynb`: Identifies GCH1 and MDR1 breakpoints.
  - `dup_trpinv_dup_breakpoints_20250128.txt`: Stores DTD breakpoints.
  - `tandem_dup_breakpoints_20250128.txt`: Stores tandem duplication breakpoints.
  - `pf-haploatlas-PF3D7_0709000_sample_summary.csv`: Pf-HaplAtlas data aiding CRT breakpoint identification.
  - `full_fws_20241125.tsv`: Pf8 fws data used to identify clonal samples for CRT breakpoints.

### **Calculate Faceaway Read Counts**
- `categorising_breakpoints/`
  - `2025_01_28_-_determine_breakpoint_read_evidence.ipynb`: Primary computational notebook.
  - `create_breakpoint_read_counts_file_20250128.py`: Script computing faceaway proportions.
  - `pf8_qc_pass_bams_20250128.tsv`: List of BAM file paths for processing.
  - `run_create_breakpoint_read_counts_file_20250128.sh`: Bash script for LSF job submission.
  - `all_samples/`: Contains per-sample `.tsv` outputs from LSF jobs.
  - `breakpoint_read_counts_noheader_20250128.tsv`: Intermediate concatenated file without headers.
  - `breakpoint_read_counts_header_20250128.tsv`: Column headers for the read counts file.

- `breakpoint_read_counts_20250128.tsv`: Final aggregated read count proportion file.