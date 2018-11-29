# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule somaticsniper:
    input:
        ref="ref.fa",
        ref_idx="ref.fa.fai",
        tumour="{sample}.bam",
        normal="RN.bam"
    output:
        "somaticsniper/{sample}.vcf"
    log:
        "log/{sample}.somaticsniper.log"
    params:
    shell:
        """
        time (
        {config[tools][somaticsniper]} \
          -n RN -t {wildcards.sample} -F vcf \
          -f {input.ref} \
          {input.tumour} \
          {input.normal} \
          {output}
        ) 1>{log} 2>&1
        """
