# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

VCFS = [
  "sim/somatic.vcf",
  "bcftools.vcf",
  "caveman.vcf",
  "haplotypecaller.filt.vcf",
#  "muclone.vcf",
  "multisnv.vcf",
  "mutect1.vcf",
  "mutect2.vcf",
  "mutect2m.vcf",
  "neusomatic.vcf",
  "shimmer.vcf",
  "snooper.fix.vcf",
  "snv-ppilp.vcf",
  "somaticsniper.vcf",
  "strelka1.fix.vcf",
  "strelka2.fix.vcf",
  "vardict.vcf",
  "varscan.vcf"
]

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

rule gatk2_collect_allelic_counts:
    input:
        ref="ref.fa",
        bam="{sample}.bam",
    output:
        "{sample}.allelicCounts.tsv"
    log:
        "log/{sample}.gatk4_alleliccounts.log"
    threads: 1
    shell:
        """
        time (
        module load gatk/4.0.10.0

        cat {VCFS} | grep -v "^#" | cut -f 1,2 \
        | while read chr pos; do printf "%s:%s\n" $chr $pos; done \
        | sort -u \
        > vars.pos.list

        {config[tools][gatk4]} CollectAllelicCounts \
            -R {input.ref} \
            -I {input.bam} \
            -L vars.pos.list \
            -O {output} \
        ) >{log} 2>&1
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
