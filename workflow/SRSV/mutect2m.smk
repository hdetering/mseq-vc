# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule mutect2m:
    input:
        ref="ref.fa",
        ref_fai="ref.fa.fai",
        ref_dict="ref.dict",
        tum1="R1.bam",
        tum2="R2.bam",
        tum3="R3.bam",
        tum4="R4.bam",
        tum5="R5.bam",
        nrm="RN.bam"
    output:
        "mutect2m.raw.vcf"
    params:
        #fastq="data/{rep}/{sam}.fq.gz",
        #sam="data/{rep}/{sam.sam"
    log:
        "log/mutect2m.log"
    threads: 1
    group:  "mutect2"
    shell:
        """
        time (
        #module load gatk/4.0.10.0
        module load jdk/8u181

        {config[tools][gatk4.1]} Mutect2 \
            -R {input.ref} \
            -I {input.tum1} \
            -I {input.tum2} \
            -I {input.tum3} \
            -I {input.tum4} \
            -I {input.tum5} \
            -I {input.nrm} \
            -normal RN \
            --output {output} \
        ) >{log} 2>&1
        """

rule mutect2m_filter:
    input:
        vcf="mutect2m.raw.vcf",
        ref="ref.fa"
    output:
        "mutect2m.vcf"
    log:
        "log/mutect2m_filter.log"
    threads: 1
    group:  "mutect2"
    shell:
        """
        time (
        #module load gatk/4.0.10.0
        module load jdk/8u181
        {config[tools][gatk4.1]} FilterMutectCalls \
            --variant {input.vcf} \
            --reference {input.ref} \
            --filtering-stats {input.vcf}.stats \
            --output {output}
        ) >{log} 2>&1
        """
