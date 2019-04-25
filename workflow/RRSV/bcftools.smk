# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#SAMPLES_ALL="RN R1 R2 R3 R4 R5".split()
SAMPLES_TUM=   "T1 T2 T3 T4 T5".split()
REF = "/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa"

rule bcftools_call:
  input:
    ref=REF,
    nrm="ChineseSon.H.bam",
    tum=["ChineseSon.%s.{scenario}.{replicate}.bam" % s for s in SAMPLES_TUM]
  output:
    "bcftools.ChineseSon.{scenario}.{replicate}.raw.vcf"
  log:
    "log/bcftools.ChineseSon.{scenario}.{replicate}.call.log"
  threads:  1
  shell:
    """
    time (
    module load bcftools/1.9

    bcftools mpileup \
      --output-type u \
      --fasta-ref {input.ref} \
      --annotate DP,AD \
      {input.nrm} {input.tum} \
    | bcftools call \
      --variants-only \
      --multiallelic-caller \
      --output-type v \
      --output {output}
    ) >{log} 2>&1
    """

rule bcftools_view:
  input:
    "bcftools.ChineseSon.{scenario}.{replicate}.raw.vcf"
  output:
    "bcftools.ChineseSon.{scenario}.{replicate}.vcf"
  log:
    "log/bcftools.ChineseSon.{scenario}.{replicate}.view.log"
  params:
    sfx="{scenario}.{replicate}"
  threads:  1
  shell:
    """
    time (
    module load bcftools/1.9

    bcftools view \
      --samples T1.{params.sfx},T2.{params.sfx},T3.{params.sfx},T4.{params.sfx},T5.{params.sfx} \
      --include 'FMT/GT[0]="0/0"' \
    {input} > {output}
    ) >{log} 2>&1
    """
