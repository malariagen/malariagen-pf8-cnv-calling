# 03_coverage_based_pipeline

## Overview

This pipeline performs copy number variation (CNV) calling based on read depth coverage across genome intervals. It is designed for Plasmodium falciparum samples and is a core component of the Pf8 data release. The pipeline leverages Nextflow, GATK's gCNV toolkit, and optionally pretrained models to efficiently process large batches of sequencing data.

Key features:
- Organizes samples into training and test sets using cluster assignments.
- Trains ploidy and CNV models on high-quality training samples.
- Applies trained models to test samples in parallel batches.
- Supports pretrained model inference for small sample sets.
- Filters low-quality or biologically implausible CNV calls using GC- and AT-content heuristics.

Outputs include gCNV denoised copy ratio files, interval- and segment-level VCFs, trained models, and intermediate files useful for QC and troubleshooting. 

## Requirements
Ensure the following files and configurations are properly set up:

- `../base.config` – Verify that the correct Nextflow profiles and resource constraints are specified for your computing cluster.
- `nextflow.config` – Ensure it is correctly set up.
- `./run.sh` – Check that file paths are correctly referenced:
  - `-profile` – Ensure the correct Nextflow profile from `../base.config` is being used.
  - `--results_dir` – Specify the correct output directory.
  - `--manifest` – Ensure a properly formatted manifest file is referenced (tab-separated, no header). See `assets_pf8/03_sample_cluster_assignment.tsv` for an example.
  - `--intervals` – Output from `01_input_generation_pipeline`: `assets_pf8/03_restricted.interval_list`.
  - `--annotated_intervals` – Output from `01_input_generation_pipeline`: `assets_pf8/03_annotated.interval_list`.
  - `--blacklist_intervals` – User-defined list of ignored intervals. See `assets_pf8/03_blacklist.interval_list` for an example.
  - `--contig_ploidy_priors` – User-defined priors. See `assets_pf8/03_contig_ploidy_priors.tsv` for an example.

## Pipeline Outputs
The pipeline generates output folders for each sample in `assets_pf8/03_execution/postprocessing/{SAMPLE_NAME}`, containing:
- `.dCR.tsv`
- `.intervals.vcf.gz` and `.intervals.vcf.gz.tbi`
- `.segments.vcf.gz` and `.segments.vcf.gz.tbi`

Additional outputs:
- `assets_pf8/03_execution/failed_test_ploidy.tsv` and `failed_train_ploidy.tsv` – Lists of sample IDs with contig ploidy ≠ 1, making CNV calling unreliable.
- `03_filtered.interval_list` – The actual list of intervals used in training, refined from `03_restricted.interval_list` by removing AT-rich bins.

For in-depth details, refer to `main.nf`, where pipeline steps are numerically annotated (e.g., `STEP 1.2`).

### Pipeline Steps Overview
1. **STEP 1.x** – Organises auxiliary files and prepares sample channels for parallel processing.
2. **STEP 2.x** – Trains models using samples marked `IN_TRAINING_SET` in `03_sample_cluster_assignment.tsv`.
3. **STEP 3.x** – Applies trained models to make CNV calls for test samples.
4. **STEP 4.x** – Evaluates the pipeline’s outputs.

### Training vs. Test Sample Processing
- Training samples within each `CLUSTER` are processed as a batch for both model fitting and inference (`FitPredict...`).
- Test samples are grouped in batches of up to 200 for inference using the trained models from `STEP 2.x`:
  ```
  | groupTuple ( size : 200, remainder : true )
  ```
- For training:
  - Contig ploidy models are fitted and used for ploidy calls.
  - These calls inform the training of the gCNV model, which then generates CNV calls.
  - Outputs pass through `gatk PostProcessGermlineCNVCalls`.
- For test samples:
  - Pre-trained contig ploidy models generate ploidy calls.
  - These calls, along with the trained gCNV model, produce CNV calls.

### Providing pretrained models
Samples can be provided paths to pretrained contig ploidy and gCNV models in `03_sample_cluster_assignment.tsv` so that CNV calls can be made without the need for training these models from scratch (which would require hundreds of samples). This feature also enables the user to perform experiments on pretrained models. Provide paths to models (found in `../assets/models`) in the `PATH_TO_PRETRAINED_CONTIG_PLOIDY_MODEL` and `PATH_TO_PRETRAINED_GCNV_MODEL` fields of `03_sample_cluster_assignment.tsv` to make CNV calls for a small number of samples. Note, that these pretrained models were pretrained using samples from Pf8 -- if the coverage of your samples is too unique when compared to the gDNA_MDA or sWGA samples used in pretraining, the calls may be unreliable and thus called as missing in stage `04_gcnv_calls_validation`. 

## Codebase Notes

### THEANO_FLAGS and Singularity Containers
Some GATK functions use:
```bash
THEANO_FLAGS="base_compiledir=tmp"
```
This workaround prevents parallel processes from clashing when accessing compiled Theano modules. Experiment with these flags to optimise performance. You may wish to consider changing it to `THEANO_FLAGS="device=cpu,base_compiledir=tmp"` or even removing it. 

### Mixing `module load` and Singularity
Both methods are currently in use due to unresolved Singularity container issues. While not ideal, this setup works for now.

### Nextflow’s `-rerun` Parameter
The `-rerun` option sometimes behaves unpredictably with LSF. Monitor its interactions on your HPC system.

### Managing `publishDir` Issues
Nextflow’s `publishDir` may terminate early, leaving partially written output files. A workaround (as seen in `FilterIntervalsList.nf`) involves renaming files with `mv` at the end of processing to trigger `publishDir` only after file writing completes.

### High File Counts
The pipeline generates hundreds of thousands of files. Consider adding `rm -rf` commands within Nextflow processes to manage disk usage. Adjust computational parameters as needed to optimise performance.

