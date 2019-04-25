# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

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
