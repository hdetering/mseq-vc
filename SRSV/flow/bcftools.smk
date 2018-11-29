# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
SAMPLE="R1 R2 R3 R4 R5 RN".split()

rule bcftools_call:
  input:
    ref="ref.fa",
    bam=["%s.bam" % s for s in SAMPLE]
  output:
    "bcftools.vcf"
  log:
    "log/bcftools.log"
  threads:  1
  shell:
    """
    time (
    module purge
    module load gcc/5.3.0 bcftools/1.4.1

    bcftools mpileup \
      --output-type u \
      --fasta-ref {input.ref} \
      --annotate DP,AD \
      {input.bam} \
    | bcftools call \
      --variants-only \
      --consensus-caller \
      --output-type v \
      --output {output}
    ) >{log} 2>&1
    """

rule bcftools_mfilt:
  input:
    "{caller}.vcf"
  output:
    "{caller}.mfilt.vcf"
  log:
    "log/{caller}.mfilt.log"
  threads:  1
  shell:
    """
    time (
    module purge
    module load gcc/5.3.0 bcftools/1.4.1

    bcftools view \
      --samples R1.variant,R2.variant2,R3.variant3,R4.variant4,R5.variant5 \
      --include "INFO/AC > 1 || FMT/AF > 0.05" \
    {input} > {output}
    ) >{log} 2>&1
    """
