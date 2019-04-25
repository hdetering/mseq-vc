# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
SAMPLES_ALL="RN R1 R2 R3 R4 R5".split()
SAMPLES_TUM=   "R1 R2 R3 R4 R5".split()

rule bcftools_call:
  input:
    ref="ref.fa",
    bam=["%s.bam" % s for s in SAMPLES_ALL]
  output:
    "bcftools.raw.vcf"
  log:
    "log/bcftools.call.log"
  threads:  1
  shell:
    """
    time (
    module load bcftools/1.9

    bcftools mpileup \
      --output-type u \
      --fasta-ref {input.ref} \
      --annotate DP,AD \
      {input.bam} \
    | bcftools call \
      --variants-only \
      --multiallelic-caller \
      --output-type v \
      --output {output}
    ) >{log} 2>&1
    """

rule bcftools_view:
  input:
    "bcftools.raw.vcf"
  output:
    "bcftools.vcf"
  log:
    "log/bcftools.view.log"
  threads:  1
  shell:
    """
    time (
    module load bcftools/1.9

    bcftools view \
      --samples R1,R2,R3,R4,R5 \
      --include 'FMT/GT[0]="0/0"' \
    {input} > {output}
    ) >{log} 2>&1
    """
