# 02_sample_clustering

## Overview
This step assigns samples to appropriate cohorts for training the HMM model in step 03. Refer to `02_sample_clustering/Pf8_manifest_generation.ipynb` for details on how best to pick samples for training. The notebook in this directory annotates `02_paths_to_readcounts.tsv` with training inclusion details and cohort assignments. The output file is `../assets_pf8/03_sample_cluster_assignment.tsv`.

## Output Format
The resultant file must contain the following columns:

- **SAMPLE_ID** – Provided by `02_paths_to_readcounts.tsv`.
- **PATH_TO_READCOUNTS** – Provided by `02_paths_to_readcounts.tsv`.
- **IN_TRAINING_SET** – Boolean (True/False) indicating if the sample is used for training the HMM.
  - Ensure a balanced number of training samples. Samuel Lee suggested avoiding under 300 samples to prevent underfitting. We train on 500 samples.
  - Overfitting risks are not well-defined, but it is advisable to train on clonal samples (From Pf8 Supplementary Materials: "_a maximum of 10 heterozygous, biallelic, coding, QC-passed SNPs with allele frequency greater than 5%_").
- **CLUSTER** – Labels differentiating sample cohorts: `gDNA` or `sWGA_MDA`.
- **PATH_TO_PRETRAINED_CONTIG_PLOIDY_MODEL** – Either "." or a path to a pretrained contig ploidy model. 
- **PATH_TO_PRETRAINED_GCNV_MODEL** – Either "." or a path to a pretrained gCNV model. 

### Additional Column (Optional)
- **IN_TRAINING_BLACKLIST** – A helper column used for logging to prevent non-clonal samples from being included in training. Not required but retained for reference.

---

### Providing pretrained models
Samples can be provided paths to pretrained contig ploidy and gCNV models in `03_sample_cluster_assignment.tsv` so that CNV calls can be made without the need for training these models from scratch (which would require hundreds of samples). This feature also enables the user to perform experiments on pretrained models. Provide paths to models (found in `../assets/models`) in the `PATH_TO_PRETRAINED_CONTIG_PLOIDY_MODEL` and `PATH_TO_PRETRAINED_GCNV_MODEL` fields of `03_sample_cluster_assignment.tsv` to make CNV calls for a small number of samples. Note, that these pretrained models were pretrained using samples from Pf8 -- if the coverage of your samples is too unique when compared to the gDNA_MDA or sWGA samples used in pretraining, the calls may be unreliable and thus called as missing in stage `04_gcnv_calls_validation`. 