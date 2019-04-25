# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

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
