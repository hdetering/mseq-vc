# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule muse_call:
    input:
        ref="ref.fa",
        ref_idx="ref.fa.fai",
        tumor="{sample}.bam",
        normal="RN.bam"
    output:
        "muse/{sample}.MuSE.txt"
    log:
        "log/{sample}.muse_call.log"
    params:
        workdir="muse",
        out_pfx="muse/{sample}"
    threads: 1
    shell:
        """
        time (

        # create output dir if necessary
        if [[ ! -d "{params.workdir}" ]]; then
          mkdir -p {params.workdir}
        fi

        {config[tools][muse]} call \
          -f {input.ref} \
          -O {params.out_pfx} \
          {input.tumor} \
          {input.normal}

       ) 1>{log} 2>&1
       """

rule muse_sump:
    input:
        calls="muse/{sample}.MuSE.txt",
    output:
        "muse/{sample}.vcf"
    log:
        "log/{sample}.muse_sump.log"
    params:
        raw_vcf="muse/{sample}.raw.vcf"
    threads: 1
    shell:
        """
        time (
        module load bcftools/1.9

        {config[tools][muse]} sump -G \
          -I {input.calls} \
          -O {params.raw_vcf}
        ) 1>{log} 2>&1

        # rename samples in VCF header (reported as TUMOR; NORMAL)
        #bcftools view -s NORMAL,TUMOR {params.raw_vcf}.gz 
        bcftools reheader \
          -s <(printf "%s\n%s\n" {wildcards.sample} RN) \
          {params.raw_vcf} \
        > {output}
        """
