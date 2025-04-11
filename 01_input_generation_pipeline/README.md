# 01_input_generation_pipeline

## Overview
This pipeline processes input data for coverage-based CNV calling. It prepares segmented reference intervals and extracts read counts from BAM files, generating key intermediate files used in later pipeline steps.

## Prerequisites
Ensure the following files are correctly configured before running the pipeline:

- **`../base.config`** – Verify that Nextflow profiles and resource constraints are properly set for your computing cluster.
- **`nextflow.config`** – Ensure this file is correctly set up.
- **`./run.sh`** – Check that the file paths are correct for required files:
  - `-profile` – Confirm the correct Nextflow profile from `../base.config` is used.
  - `--results_dir` – Set this to the desired output directory.
  - `--manifest` – Provide a properly formatted manifest file:
    - Must be a tab-separated file containing sample names and BAM file paths.
    - Should not have a header.
    - See `assets_pf8/01_paths_to_bams.tsv` for an example.
  - `--reference` – Reference the correct *Pf3D7* genome `.fasta` file:
    - Ensure `.fasta.fai` and `.dict` files are in the same directory.
  - `--intervals_of_interest` – Provide a user-defined list of genomic regions of interest:
    - Example format available in `assets_pf8/01_intervals_of_interest.interval_list`.

## Execution
Run the pipeline with the following commands:
```sh
chmod 744 run.sh
./run.sh
```

## Pipeline Workflow
1. **Segmentation**: The *Pf3D7* reference genome is divided into 500bp bins.
2. **Filtering**: Only bins within the specified `01_intervals_of_interest.interval_list` are retained, producing:
   - **`03_restricted.interval_list`** – Defines the bins for read count extraction.
3. **Annotation**: GC content is computed and added to the filtered bins, generating:
   - **`03_annotated.interval_list`** – Used in later steps to prevent training on AT-rich regions.

## Output
The pipeline creates a subdirectory within `assets_pf8` (or the specified `results_dir`) called `01_execution`, containing:

- **`timeline.html`** – Nextflow execution report.
- **`readcounts/`** – Contains per-sample read count files.
- **`02_paths_to_readcounts.tsv`** – Maps sample IDs to their respective `.counts.tsv` files.

These outputs are essential for subsequent steps in the CNV calling pipeline.

