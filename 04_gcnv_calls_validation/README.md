# 04_gcnv_calls_validation

---

## Overview

This subdirectory contains scripts and workflows for processing raw pipeline outputs to obtain final CNV genotype calls as part of the Pf8 CNV calling process. The workflow consists of three main steps:

1. **Genotyping CNV Calls**: Using `cnv_coverage_genotype_caller.ipynb`, raw coverage pipeline outputs (`.dCR.tsv`, `intervals.vcf`, and `segments.vcf`) are converted into genotype (GT) calls {-1, 0, 1} (where 1 denotes an amplification/deletion, 0 denotes lack of amplification/deletion, and -1 denotes uncallable). A breakdown of these file types is provided in the appendix of this README.

2. **Generating Diagnostic Plots**: `generate_plots.py` and `run.sh` are used to create diagnostic plots for each gene and sample.

3. **Identifying Discordance and Manual Curation**: Using `app.py`, discrepancies with previous CNV callsets (Pf7 final calls, Pf7 faceaway calls, or Broad pipeline coverage calls) are identified. A subset of potentially low-confidence calls is manually reviewed and modified using diagnostic plots.

---

## Detailed instructions

#### 1. Running the genotyper

`cnv_coverage_genotype_caller.ipynb` contains heuristics to determine CNV genotypes, following rules outlined in the Pf8 Supplementary Materials:
> The raw coverage pipeline outputs were then processed using a custom genotyping algorithm so that a single genotype value was created for each call region in a sample. For gene amplification calls, the genotype values could take one of three values: missing/uncallable (-1), not amplified (0), or amplified (1). In mixed samples, the copy ratio must have exceeded 1.5 for it to have been called as 1 (i.e., at least 50% of the parasites would be expected to contain an amplification). For gene deletion calls, the genotype values could be one of three values: missing/uncallable (-1), not deleted (0), or deleted (1). A call region’s genotype value automatically designated -1 if the number of segments in the corresponding contig exceeded 5, if the contig’s ploidy value was not 1, or if less than 20% of the contig was occupied by segments containing copy number values not equal to 1. The genotyping algorithm then proceeded differently, depending on whether the call region required a deletion call or an amplification call.
> 
>  For deletion calls, the genotyping algorithm only used data 5’ of, or within the call region for both the previous and subsequent genotyping steps, in order to reduce noise introduced by the subtelomeric regions 3’ of the call region. If the mean absolute residuals between the denoised copy ratio values and the corresponding segment’s copy number (CN) value was no greater than 0.3, segments which occupied at least part of the call region and contained a QS value (further details available in Supplementary figure 1) of at least 400 were used to determine the genotype (GT) value, otherwise the following rescue condition was triggered: if at least 80% of the contig was occupied by segments which had both a CN value of 1 and a QS value of 1000, and at least 20% of the call region was occupied by a segment with a CN value not equal to 1, then the segment was used to determine the GT value. This was done to save samples with true deletions from the tendency for GATK PostprocessGermlineCNVCalls to strongly penalise the QS value of segments with CN values not equal to 1. The genotypes for the raw coverage pipeline outputs for the gene deletion calls can be found in the columns with suffix “_uncurated_coverage_only”.
> 
> For amplification calls, if the mean absolute residuals between the denoised copy ratio values and the corresponding segment’s CN value was no greater than 0.3, segments which occupied at least 75% of the call region, contained a QS value of at least 400 and a QA value greater than 0 were used to determine the GT value, otherwise the following rescue condition was triggered: if at least 80% of the contig was occupied by segments which had both a CN value of 1 and a QS value of 1000, and at least 75% of the call region was occupied by a segment with a CN value not equal to 1, then the segment was used to determine the GT value. This was done to save samples with true amplifications from the tendency for GATK PostprocessGermlineCNVCalls to strongly penalise the QS value of segments with CN values not equal to 1. In all other cases, the GT was -1. The genotypes for the raw coverage pipeline outputs for the gene amplification calls can be found in the columns with suffix “_uncurated_coverage_only”. 

The output of the genotyper is the file `app_files/draft_coverage_calls.tsv`. This notebook also performs a series of analyses to group samples based on their Pf7 CNV call and their new draft Pf8 call, and prints out a series of Python lists of sample IDs, which can easily be copied into the Streamlit app for inspecting just those samples (see further below for details). This way of working allowed for a highly dynamic, flexible and rapidly iterating way of interrogating different sets of samples. 

#### 2. Generating diagnostic plots

To generate plots, run `./run.sh`, which submits an LSF array job invoking `generate_plots.py`. The resulting plots are saved in `./results/` and also available at the [Pf8 data resource S3 bucket](https://pf8-release.cog.sanger.ac.uk/cnv-diagnostic-plots/index.html). Adjust `CHUNK_SIZE` to control batch processing. 

#### 3. Identifying discordance and manually curating the calls

Install dependencies from `requirements.txt`, then `cd` to this directory and run `streamlit run app.py` or `python3 -m streamlit run app.py` or equivalent. The terminal should list several localhost URLs which you can click to access the app. You may wish to add the flag `--server.fileWatcherType "poll"` if the app is taking too long to load (this was the case on Sanger farm22 but not gen3). 

The app has five "browser modes":
- **"General browsing"** allows the user to browse all the diagnostic plots for all the genes. There is a toggle for restricting to Pf7's QC-passed samples only. Each page additionally displays calls from the Broad's Pf7 pipelines (named "Interim1" and "Interim2") as well as Pf7 calls and Pf8's draft calls. There is a textbox available where the user can paste Python lists of sample IDs in to. This subsets the samples displayed, allowing for rapid iteration between investigating samples by their meta data on Jupyter notebooks to then browsing them in the Streamlit app to find abnormalities or biases. The user can then navigate through the collection of samples using the "Next" and "Previous" buttons. Where available, a plot is shown containing the dCR outputs from the Broad's pipeline below the diagnostic plot. 
- **"Modify draft calls"** allows the user to see the distribution of samples possessing different combinations of genotype calls from the Pf7 CNV callset and the draft Pf8 callset. Of greatest interest are the samples that changed from GT=0 to GT=1 or vice versa (as these indicate a change from no amplification/deletion to amplification/deletion). The user can click on a row of the dataframe to, for example, show all samples that changed from GT=0 to GT=1. The user is then provided the option to manually change the call if deemed appropriate. Any modifications are stored **in memory** as a json file. Clicking `Save changes_checkpoint.json` allows the user to export a checkpoint save file storing all the changes made so far, so that the user can pick up from where they left off using `Load changes_checkpoint.json`. Finally, once the user is content with all of the manual changes, the manual modifications can be integrated with `draft_coverage_calls.tsv` to create `final_coverage_calls.tsv` using the "Convert manual changes to final_coverage_calls.tsv" button. This resulting output file is essentially the final coverage-based calls. 
  - Notebooks used to assist with manual curation can be found in `analysis_notebooks/`
  - From Pf8's Supplementary Materials:
      > For both deletion and amplification calls, samples with GT values different from Pf7 were manually curated to determine whether a manual modification to the GT value was needed using diagnostic plots, if the total number of samples containing that combination of old and new GT values was no greater than 100. For crt and gch1 amplification calls, samples were manually curated for potential manual GT modification if the total number of samples containing that GT value was no greater than 100. [...] We emphasise that the manual curation step was designed to target only the most difficult of coverage-based calls. We additionally inspected thousands of diagnostic plots of samples that were not curated to help confirm that the final GT values were of satisfactory quality.
      > 
      > [...]
      >
      > When curating gene deletion calls, sequencing read breakpoints were looked for in new samples identified to have potential gene deletions using IGV genome browser. Samples with both evidence from the coverage pipeline and from breakpoint evidence (either from Pf7 or from IGV) were called as a deletion in the final deletion calls columns suffixed “_final_deletion_call” with an additional “_breakpoint” (showing the genomic coordinates of the breakpoint) and “_deletion_type” column (describing the molecular nature of the deletion, e.g., recombination with chromosome 11, telomere healing, etc.). Samples with evidence from the coverage pipeline but without breakpoint evidence were called as -1 (missing). 
- **"Modification stats"** allows the user to see an interactive pie chart of how the manual curation stage went. You can hover over and click on many of the elements. The Streamlit app also automatically saves the chart as `final_coverage_results.html`. 
- **"Inspect modifications"** allows the user to flit through all of the samples that had their genotype changed. 
- **"Inspect final calls"** allows the user to inspect samples in groups of unique genotype combinations using the final calls, similar to as was done in "Modify draft calls", but this time with the final manually curated calls. 

---

## Full contents of this subdirectory

## Directory Contents

### Diagnostic Plot Generation
- `generate_plots.py` - Generates diagnostic plots per sample.
- `run.sh` - Submits array jobs for batch processing.
- `results/` - Stores output plots (not included in repo).

### CNV Call Processing
- `cnv_genotype_caller.ipynb` - Generates `draft_coverage_calls.tsv`.
- `app.py` - Streamlit app for manual curation.
- `.streamlit/config.toml` - App UI configuration.
- `app_files/`
  - `draft_coverage_calls.tsv` - Initial genotype calls before manual review.
  - `Pf7_duplication_calls_*.tsv` - Previous Broad Institute CNV calls.
  - `Pf7_inferred_resistance_status_classification.tsv` - Drug resistance classification.
  - `final_coverage_calls.tsv` - Final curated CNV calls.
  - `final_coverage_results.html` - Summary of manual modifications.

### Manual Curation Notebooks
- `analysis_notebooks/`
  - `2024_12_01_-_breakpoint_hrp2.ipynb` - Investigates HRP2 deletions.
  - `2024_11_30_-_breakpoint_hrp3.ipynb` - Inspects HRP3 deletions and generates `hrp_calls_pf8_20241220.tsv`.
  - `hrp_calls_pf7.tsv` : Pf7's HRP breakpoints file
  - `hrp_calls_pf8.tsv` : Pf8's HRP breakpoints file, generated using the notebook `2024_11_30_-_breakpoint_hrp3.iypnb`
  - `hrp_calls_pf8_20241220.tsv` : Same as `hrp_calls_pf8.tsv` but renamed for sharing
  - `2024_12_02_pf8_coverage_calls.tsv` : Same as `final_call.tsv` but renamed for sharing

- Miscellaneous
  - `example_data/` : contains examples of .dCR.tsv, .intervals.vcf, segments.vcf

---

## History

Initially, coverage-based calls were to be used as final CNV calls. However, late-stage analysis revealed underreporting of plasmepsin2/3 duplications in AS-SE-E samples. This led to the reintroduction of the Pf7 faceaway method, requiring two additional weeks of work.

Chiyun currently believes that the prevalence of plasmepsin duplications are so high in the AS-SE-E-looking coverage profiles of samples that the model mistook the duplicated amount of DNA for the CN=1 amount. To be continued as to whether this is correct. In future, more sophisticated training sample selection techniques may be needed to overcome this. 

---

# Appendix

To better understand how `.dCR.tsv`, `intervals.vcf` and `segments.vcf` relate to each other, please take a look at some example files in the `example_data/` directory. 

---

### Example of dCR.tsv - see `/example_data/PD1033-C.dCR.tsv` for full file

| CONTIG       | START   | END     | LINEAR_COPY_RATIO      |
|--------------|---------|---------|------------------------|
| Pf3D7_05_v3  | 635001  | 635500  | 1.0654295637909725     |
| Pf3D7_05_v3  | 635501  | 636000  | 1.141918548448794      |
| Pf3D7_05_v3  | 636001  | 636500  | 1.0602725935939856     |
| Pf3D7_05_v3  | 636501  | 637000  | 0.8658393171405365     |

There is one floating-point copy-ratio value for each 500-bp bin. For Pf8, there were 4850 rows per sample, because there were 1000 500-bp bins for each of the 6 genes (so 6000 in total), before FilterIntervals removed some bins due to noise. 

---

### Example of intervals.vcf - see `/example_data/PD1033-C.intervals.vcf` for full file

| CHROM       |	POS    | ... | INFO	      | FORMAT         | PD1033-C                      |
|-------------|--------|-----|------------|----------------|-------------------------------|
| Pf3D7_05_v3 |	635001 | ... | END=635500	| GT:CN:CNLP:CNQ | 0:1:368,0,123,216,299,367:123 |
| Pf3D7_05_v3 | 635501 | ... | END=636000 | GT:CN:CNLP:CNQ | 0:1:368,0,92,120,146,169:92   |
| Pf3D7_05_v3 | 636001 | ... | END=636500 | GT:CN:CNLP:CNQ | 0:1:368,0,109,158,202,240:109 |
| Pf3D7_05_v3 | 636501 | ... | END=637000 | GT:CN:CNLP:CNQ | 0:1:368,0,116,159,196,227:116 |

There is one integer copy-number value (CN) for each 500-bp bin. For Pf8, there were 4850 rows per sample, because there were 1000 500-bp bins for each of the 6 genes (so 6000 in total), before FilterIntervals removed some bins due to noise. 

---

### Example of segments.vcf - see `/example_data/PD1033-C.segments.vcf` for full file
| #CHROM      | POS    | INFO        | FORMAT                 | PD1033-C               |
|-------------|--------|-------------|------------------------|------------------------|
| Pf3D7_05_v3 | 635001 | END=947500  | GT:CN:NP:QA:QS:QSE:QSS | 0:1:548:58:3077:55:124 |
| Pf3D7_05_v3 | 947501 | END=954000  | GT:CN:NP:QA:QS:QSE:QSS | 1:2:13:49:618:45:51    |
| Pf3D7_05_v3 | 954001 | END=962000  | GT:CN:NP:QA:QS:QSE:QSS | 1:3:11:45:347:38:45    |
| Pf3D7_05_v3 | 962001 | END=970000  | GT:CN:NP:QA:QS:QSE:QSS | 1:2:14:38:657:51:38    |
| Pf3D7_05_v3 | 970001 | END=1135000 | GT:CN:NP:QA:QS:QSE:QSS | 0:1:277:47:3077:94:51  |

Consecutive regions of intervals which suggest the same CN are concatenated together to form long segments along the contig. Note that the first segment ranges from position 635,001 to 947,500 and is CN=1. Also, note this is a DUP-TRP/INV-DUP variant of MDR1 so the CN changes from 1 to 2 to 3 to 2 to 1. 