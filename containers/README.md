# Creating GATK containers for Nextflow pipelines

Needed for `01_input_generation_pipeline` and `03_coverage_based_pipeline`

Nextflow base config file, `base.config` requires there to be a container at this location `container = "$projectDir/../containers/gatk_4.5.0.0.sif"`. Run `run.sh` to pull the Singularity container. 