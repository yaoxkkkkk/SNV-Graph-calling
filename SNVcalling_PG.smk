import os
import gzip

# # # 提取文件名的基部分（去除路径和扩展名）
graph_basename=os.path.splitext(os.path.basename(config["graph_gfa"]))[0]
ref_basename=os.path.splitext(os.path.basename(config["haplotype_genome"]))[0]
fastq_suffix=config.get("fastq_suffix")

rule all:
    input:
        f"genome_index/{graph_basename}.dist",
        f"genome_index/{graph_basename}.shortread.withzip.min",
        f"genome_index/{graph_basename}.shortread.zipcodes",
        f"genome_index/{graph_basename}.gbz",
        f"genome_index/{graph_basename}.hapl",
        f"genome_index/{graph_basename}.path_list.txt",
        f"genome_index/{ref_basename}.fai",
        expand("kff/{sample}.fq.kmer.kff", sample=config["sample"]),
        expand("mapping/{sample}.sorted.bam.bai", sample=config["sample"]),
        expand("vcf/gvcf/{sample}.g.vcf.gz.csi", sample=config["sample"])
        "vcf/cohort.vcf.gz",
        "vcf/snv.core.vcf.gz"

rule vg_autoindex_giraffe:
    input:
        gfa_file=config["graph_gfa"]
    output:
        dist=f"genome_index/{graph_basename}.dist",
        min=f"genome_index/{graph_basename}.shortread.withzip.min",
        zipcodes=f"genome_index/{graph_basename}.shortread.zipcodes",
        gbz=f"genome_index/{graph_basename}.giraffe.gbz"
    threads: 10
    params:
        prefix=f"genome_index/{graph_basename}"
    conda:
        config["conda_env"]["vg_env"]
    log:
        f"logs/autoindex_{graph_basename}.log"
    shell:
        """
        vg autoindex \
        --workflow giraffe \
        -t {threads} \
        -g {input.gfa_file} \
        -p {params.prefix} \
        2> {log}
        """

rule gbz_rename:
    input:
        gbz=f"genome_index/{graph_basename}.giraffe.gbz"
    output:
        rename_gbz=f"genome_index/{graph_basename}.gbz"
    shell:
        """
        mv {input.gbz} {output.rename_gbz}
        """

rule RIndex:
    input:
        gbz=f"genome_index/{graph_basename}.gbz"
    output:
        ri=f"genome_index/{graph_basename}.ri"
    threads: 10
    conda:
        config["conda_env"]["vg_env"]
    log:
        "logs/rindex.log"
    shell:
        """
        vg gbwt \
        -p \
        --num-threads {threads} \
        -r {output.ri} \
        -Z {input.gbz} \
        2> {log}
        """

rule HaplotypeIndex:
    input:
        f"genome_index/{graph_basename}.gbz",
        f"genome_index/{graph_basename}.ri"
    output:
        f"genome_index/{graph_basename}.hapl"
    threads: 10
    conda:
        config["conda_env"]["vg_env"]
    log:
        "logs/haplotype_index.log"
    shell:
        """
        vg haplotypes \
        -H {output} \
        -t {threads} \
        {input[0]} \
        2> {log}
        """

rule VgPath:
    input:
        gbz=f"genome_index/{graph_basename}.gbz"
    output:
        path_list=f"genome_index/{graph_basename}.path_list.txt"
    params:
        config["haplotype_code"]
    threads: 8
    conda:
        config["conda_env"]["vg_env"]
    log:
        "logs/vg_path.log"
    shell:
        """
        vg paths \
        -x {input.gbz} \
        -L \
        -Q {params} \
        1> {output.path_list} \
        2> {log}
        """

rule samtools_fai_index:
    input:
        reference_genome=config["haplotype_genome"]
    output:
        "genome_index/{ref_basename}.fai"
    log:
        "logs/index/samtools_index_{ref_basename}.log"
    shell:
        """
        samtools faidx {input.reference_genome} 2> {log}
        cp {input.reference_genome}.fai {output}
        """

rule QualityControlfastp:
    input:
        f"raw_data/{{sample}}_1{fastq_suffix}",
        f"raw_data/{{sample}}_2{fastq_suffix}"
    output:
        "clean_data/{sample}_1_clean.fq.gz",
        "clean_data/{sample}_2_clean.fq.gz",
        "logs/fastp/fastp_report/{sample}.fastp.html"
    threads: 2
    params:
        qualified_quality_phred=config["qualified_quality_phred"],
        unqualified_percent_limit=config["unqualified_percent_limit"],
        trim_front=config["trim_front"]
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp \
        --thread {threads} \
        -i {input[0]} \
        -I {input[1]} \
        -o {output[0]} \
        -O {output[1]} \
        -h {output[2]} \
        -j /dev/null \
        -q {params.qualified_quality_phred} \
        -u {params.unqualified_percent_limit} \
        -f {params.trim_front} \
        2> {log}
        """

rule FastqFileList:
    input:
        r1="clean_data/{sample}_1_clean.fq.gz",
        r2="clean_data/{sample}_2_clean.fq.gz"
    output:
        "list/{sample}.fq.list"
    shell:
        """
        cat > {output} <<- EOM
$(readlink -f {input.r1})
$(readlink -f {input.r2})
EOM
        """

rule FastqKmerCount:
    input:
        fq_list="list/{sample}.fq.list"
    output:
        kff_file="kff/{sample}.fq.kmer.kff",
        tmp_dir=temp(directory("kff/{sample}_kmer_temp/"))
    threads: 8
    params:
        k=config["kmer_length"],
        prefix="kff/{sample}.fq.kmer"
    log:
        "logs/kmer_count/{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.kff_file})
        mkdir -p {output.tmp_dir}
        kmc \
        -k{params.k} \
        -okff \
        -t{threads} \
        @{input.fq_list} \
        {params.prefix} \
        {output.tmp_dir} \
        &> {log}
        """

rule Giraffe_map:
    input:
        path_list=f"genome_index/{graph_basename}.path_list.txt",
        r1="clean_data/{sample}_1_clean.fq.gz",
        r2="clean_data/{sample}_2_clean.fq.gz",
        gbz=f"genome_index/{graph_basename}.gbz",
        kff="kff/{sample}.fq.kmer.kff",
        hapl=f"genome_index/{graph_basename}.hapl"
    output:
        bam_file="mapping/{sample}.sorted.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA",
        haplotype_code=config["haplotype_code"]
    threads: 4
    conda:
        config["conda_env"]["vg_env"]
    log:
        "logs/giraffe/giraffe_map_{sample}.log"
    shell:
        """
        vg giraffe \
        --progress \
        --read-group {params.rg} \
        --sample {wildcards.sample} \
        -o SAM \
        --ref-path {input.path_list} \
        -P \
        -L 3000 \
        -f {input.r1} \
        -f {input.r2} \
        -Z {input.gbz} \
        --kff-name {input.kff} \
        --haplotype-name {input.hapl} \
        -t {threads} \
        | sed -e "s/{params.haplotype_code}#0#//g" \
        | samtools sort --threads {threads} -m 2G -O BAM \
        -o {output.bam_file} \
        2> {log}
        rm genome_index/*.{wildcards.sample}.*
        """

rule BamIndex:
    input:
        "mapping/{sample}.sorted.bam"
    output:
        "mapping/{sample}.sorted.bam.bai"
    threads: 2
    shell:
        """
        samtools index \
        -@ {threads} \
        {input} \
        2> /dev/null
        """

rule GVCFCalling:
    input:
        haplotype_file=config["haplotype_genome"],
        bam="mapping/{sample}.sorted.bam",
        bam_index="mapping/{sample}.sorted.bam.bai",
        ref_index=f"genome_index/{ref_basename}.fai"
    output:
        gvcf="vcf/gvcf/{sample}.g.vcf.gz"
    params:
        meea=r"min_mapping_quality=0,keep_legacy_allele_counter_behavior=true,normalize_reads=true",
        tmpdir=lambda wildcards: f"tmpdir/{wildcards.sample}"
    container:
        config["singularity"]["deepvariant"]
    threads: 4
    log:
        "logs/vcf/gvcf/{sample}.gvcf.log"
    shell:
        """
        run_deepvariant \
        --model_type=WGS \
        --ref {input.haplotype_file} \
        --reads={input.bam} \
        --output_gvcf={output.gvcf} \
        --output_vcf=/dev/null \
        --make_examples_extra_args={params.meea} \
        --num_shards={threads} \
        --intermediate_results_dir {params.tmpdir} \
        &> {log}
        """

rule index_gvcf:
    input:
        "vcf/gvcf/{sample}.g.vcf.gz"
    output:
        "vcf/gvcf/{sample}.g.vcf.gz.csi"
    log:
        "logs/vcf/gvcf/{sample}.g.vcf.index.log"
    params:
        extra="--csi"
    shell:
        """
        bcftools index {params.extra} {input} 2> {log}
        """

rule GVCFJointCalling:
    input:
        gvcfs=expand("vcf/gvcf/{sample}.g.vcf.gz", sample=config["sample"]),
        gvcf_indices=expand("vcf/gvcf/{sample}.g.vcf.gz.csi", sample=config["sample"])
    output:
        vcfs="vcf/cohort.vcf.gz",
        idx="vcf/cohort.vcf.gz.tbi",
        dbdir=directory("vcf/cohort.DB")
    params:
        cfg="DeepVariantWGS"
    threads: 192
    conda:
        config["conda_env"]["glnexus_env"]
    log:
        "logs/vcf/joint/glnexus.log"
    shell:
        """
        glnexus_cli \
          --config {params.cfg} \
          --threads {threads} \
          --dir {output.dbdir} \
          {input.gvcfs} \
        2> {log} \
        | bcftools view -O z -o {output.vcfs}
        tabix -f {output.vcfs}
        """

rule FilterCoreVariants:
    input:
        vcf="vcf/cohort.vcf.gz",
        idx="vcf/cohort.vcf.gz.tbi"
    output:
        vcf="vcf/snv.core.vcf.gz",
        idx="vcf/snv.core.vcf.gz.tbi"
    log:
        "logs/core_filter_variants.log"
    threads: 20
    params:
        expr=lambda wildcards: (
            f"QUAL > 20.0 & MAF > {config['core']['maf']} "
            f"& F_MISSING < {config['core']['missingrate']} & N_ALT = 1"
        )
    shell:
        """
        bcftools view \
          -i '{params.expr}' \
          --threads {threads} \
          {input.vcf} \
        | bcftools view -e 'GT~"\./[0-9]" || GT~"[0-9]/\."' \
        | bcftools annotate -x INFO,FORMAT -O z \
          -o {output.vcf} \
        2> {log}
        tabix -f {output.vcf}
        """