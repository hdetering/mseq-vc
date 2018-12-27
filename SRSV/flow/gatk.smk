# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule ref_index:
    input:
        "{ref}.fa"
    output:
        fai="{ref}.fa.fai",
        dict="{ref}.dict"
    log:    "log/{ref}.ref_index.log"
    threads: 1
    group:  "mutect1"
    shell:
        """
        time (
        module load gcc/6.4.0 
        module load samtools/1.8
        module load jdk/8u181
        #module load picard/2.17.11

        samtools faidx {input} 2> {log}
        {config[tools][picard]} CreateSequenceDictionary \
            REFERENCE={input} \
            OUTPUT={output.dict}
        ) >{log} 2>&1
        """

rule mutect2:
    input:
        ref="ref.fa",
        ref_fai="ref.fa.fai",
        ref_dict="ref.dict",
        tumor="{sam}.bam",
        normal="RN.bam"
    output:
        "mutect2/{sam}.raw.vcf"
    params:
        #fastq="data/{rep}/{sam}.fq.gz",
        #sam="data/{rep}/{sam.sam"
    log:
        "log/{sam}.mutect2.log"
    threads: 1
    group:  "mutect2"
    shell:
        """
        time (
        module load gatk/4.0.10.0
        #module load jdk/8u181
        if [[ ! -d mutect2 ]]; then mkdir mutect2; fi

        {config[tools][gatk4]} Mutect2 \
            -R {input.ref} \
            -I {input.tumor} \
            -tumor {wildcards.sam} \
            -I {input.normal} \
            -normal RN \
            -O {output} \
        ) >{log} 2>&1
        """

rule mutect2_filter:
    input:
        "mutect2/{sam}.raw.vcf"
    output:
        "mutect2/{sam}.vcf"
    log:
        "log/{sam}.mutect2_filter.log"
    threads: 1
    group:  "mutect2"
    shell:
        """
        time (
        module load gatk/4.0.10.0
        #module load jdk/8u181
        {config[tools][gatk4]} FilterMutectCalls \
            -V {input} \
            -O {output}
        ) >{log} 2>&1
        """

rule mutect1:
    input:
        ref="ref.fa",
        ref_fai="ref.fa.fai",
        ref_dict="ref.dict",
        tumor="{sam}.bam",
        normal="RN.bam"
    output:
        "mutect1/{sam}.vcf"
    params:
        workdir="mutect1"
    log:
        "log/{sam}.mutect1.log"
    threads: 1
    group:  "mutect1"
    shell:
        """
        time(
        #module load jdk/1.7.0
        if [[ ! -d "{params.workdir}" ]]; then mkdir -p {params.workdir}; fi

        {config[tools][mutect1]} \
            --analysis_type MuTect \
            --reference_sequence {input.ref} \
            --input_file:tumor {input.tumor} \
            --input_file:normal {input.normal} \
            -vcf {output} \
        ) 1>{log} 2>&1
        """

# Merge VCF files from multiple samples into single VCF
rule merge_vcf:
    input:
       ref="ref.fa",
       vcf1="{tool}/R1.vcf",
       vcf2="{tool}/R2.vcf",
       vcf3="{tool}/R3.vcf",
       vcf4="{tool}/R4.vcf",
       vcf5="{tool}/R5.vcf"
    output: 
       "{tool,[^\./]+}.vcf"
    log:
       "log/{tool}.GATK_CombineVariants.log"
    shell:
        """
        #module load jdk/8u171
        module load jdk/8u181
        time (
           {config[tools][gatk3]} -T CombineVariants \
             -R {input.ref} \
             --variant {input.vcf1} \
             --variant {input.vcf2} \
             --variant {input.vcf3} \
             --variant {input.vcf4} \
             --variant {input.vcf5} \
             -o {output} \
             -genotypeMergeOptions UNIQUIFY
        ) 1>{log} 2>&1
        """

# Add read group to BAM files (required by GATK tools)
rule add_rg:
    input:  "data/{rep}/{sam}.bam"
    output: "data/{rep}/{sam}.RG.bam"
    log:    "data/{rep}/log/{sam}.picard_AddRG.log"
    group:  "mutect1"
    shell:
        """
        module load jdk/9.0.4
        java -jar {config[tools][picard]} AddOrReplaceReadGroups \
            I={input} \
            OUTPUT={output} \
            RGID={wildcards.sam} \
            RGSM={wildcards.sam} \
            RGPL="ILLUMINA" \
            RGLB="ART_500_pe" \
            RGPU="simulation" \
            CREATE_INDEX=true
        """

# Remove PCR duplicates
rule rm_dup:
    input:  "data/{rep}/{sam}.bam"
    output: "data/{rep}/{sam}.rmdup.bam"
    log:    "data/{rep}/log/{sam}.picard_MarkDup.log"
    params:
        tmp_dir="/tmp/{rep}",
        dupfile="data/{rep}/{sam}.dup.txt"
    shell:
        """
        module load jdk/9.0.4
        #module load picard/2.2.1

        time(
        java -jar {config[tools][picard]} MarkDuplicates \
            INPUT={input} \
            OUTPUT={output} \
            CREATE_INDEX=true \
            REMOVE_DUPLICATES=true \
            TMP_DIR={params.tmp_dir} \
            M={params.dupfile} \
            VALIDATION_STRINGENCY=LENIENT \
        > {log} 2>&1
        )
        """
