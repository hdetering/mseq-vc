# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule vardict:
    input:
        ref="ref.fa",
        ref_fai="ref.fa.fai",
        bed="ref.bed",
        tumor="{sam}.bam",
        normal="RN.bam"
    output:
        vcf="vardict/{sam}.vcf"
    params:
        af_min=0.01
    log:
        out="log/{sam}.vardict.out",
        err="log/{sam}.vardict.err"
    threads: 1
    group:  "vardict"
    shell:
        """
        #module load gcc/6.4.0 vardict/20180502
        #module load intel/2016 vardict/1.0

        time(
            module load jdk/8u181
            {config[tools][vardictjava]} \
              -G {input.ref} \
              -f {params.af_min} \
              -N "{wildcards.sam}" \
              -b "{input.tumor}|{input.normal}" \
              -c 1 -S 2 -E 3 {input.bed} \
            > {wildcards.sam}.vardict.tsv

            module load gcc/6.4.0 R/3.5.1
            cat {wildcards.sam}.vardict.tsv \
            | {config[tools][vardictscripts]}/testsomatic.R \
            > {wildcards.sam}.testsom.tsv

            cat {wildcards.sam}.testsom.tsv \
            | {config[tools][vardictscripts]}/var2vcf_paired.pl \
              -N "{wildcards.sam}|RN" \
              -f {params.af_min} \
            > {output.vcf}
        ) 1>{log.out} 2>{log.err}
        """
