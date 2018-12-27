# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule snv_ppilp:
    input:
        "haplotypecaller.filt.vcf"
    output:
        "snv-ppilp.out.txt"
    log:
        "log/snv-ppilp.log"
    params:
    threads: 1
    shell:
        """
        time (
        python2 {config[tools][snv-ppilp]} \
          -i {input} \
          -o {output}
        ) 1>{log} 2>&1
        """
