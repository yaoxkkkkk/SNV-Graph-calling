# SNV calling based on genome graph

## Dependent softwares

- [fastp](https://github.com/OpenGene/fastp)
- [vg](https://github.com/vgteam/vg)
- [KMC](https://github.com/refresh-bio/KMC)
- [DeepVariant](https://github.com/google/deepvariant)
- [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- [SAMtools](https://github.com/samtools/samtools)
- [BCFtools](https://github.com/samtools/bcftools)

## What to input
- Genome graph file in GFA format
- Reference haplotype genome which the graph build on
- WGS fastq files

## What to output
- SNV set (missingrate < 0.2, MAF > 0.05, bi-allelic sites)

## Usage

### 1. Prepare your working directory

```shell
├── raw_data
├── genome_index
└── logs
```

Please storage your resequence data in `raw_data/` folder, graph file and genome file in `genome_index/` folder. Script files, pipeline files and configuration files can be stored the way you like.

### 2. Prepare the config file

The config file needs to be at the same folder of snakefile.

#### 2.1 Move the graph file and genome file to `genome_index/` folder, add the file absolute path like:

```shell
# Absolute path to the pangenome gfa file
graph_gfa: "/genome_index/path/to/gfa"

# Reference haplotype code that pangenome graph built on
haplotype_code: ""

# Absolute path to reference haplotype genome file
haplotype_genome: "/genome_index/path/to/fa"
```

#### 2.2 Sometimes the fastq files may be ended with `.fastq.gz` or `.fq.gz`, specify the suffix of the fastq files if it's necessary.

```shell
# Fastq file suffix
fastq_suffix: ".fq.gz" # Default value is ".fq.gz"
```

#### 2.3 Software configuration

In this workflow:

- `vg` and `GLnexus` are installed via mamba.
- `DeepVariant` is executed through a Singularity container.

Please make sure to:

- Update the conda environment paths for `vg` and `GLnexus`.
- Set the absolute path to the Singularity image for DeepVariant.

If you use different installation methods, modify the Snakemake file accordingly to match your environment.

#### 2.4 Fill in the name of the samples. The samples name need to be filled with specific format like:

```shell
# Sample list, samples' name should start with letters.
sample:
    - "sample1"
    - "sample2"
    - "sample3"
    - "sample4"
    - ...
    - "samplen"
```

You can use following command to add sample list to the config file if you have a sample list txt file (for example `sample.list`):

```shell
# sample.list
sample1
sample2
sample3
sample4

# Add samples to the config file:
awk '{print "    - \"" $0 "\""}' sample.list >> ${working_dir}/SNPcalling_config.yaml
```

### 3. Submit the pipeline to HPC cluster

Put snakefile and configuration file in the same directory, then running it.

For example:

```bash
snakemake \
	--snakefile ${snakefile} \
    --configfile ${configfile} \
	-d ${working_dir} \
	--cores ${cores_num} \
	--use-conda \
	--use-singularity \
	--rerun-incomplete \
	--nolock
```