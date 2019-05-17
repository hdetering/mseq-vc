# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule samtools_mpileup:
  input:
    ref="ref.fa",
    bam="{sample}.bam"
  output:
    pup="varscan/{sample}.pu"
  log:"log/{sample}.varscan_mpileup.log"
  threads: 1
  shell:
    """
    module load gcc/6.4.0 samtools/1.9

    time (
    samtools mpileup \
      --fasta-ref {input.ref} \
      --min-MQ 40 \
      {input.bam} \
      > {output.pup}
    ) >{log} 2>&1
    """ 

rule varscan:
  input:
    ref="ref.fa",
    tum="{sample}.bam",
    nrm="varscan/RN.pu"
  output:
    vcf="varscan/{sample}.raw.vcf",
    snp="varscan/{sample}.snp.vcf",
    ind="varscan/{sample}.indel.vcf"
  params:
    pfx="varscan/{sample}"
  log: "log/varscan.{sample}.log"
  threads: 2
  shell:
    """
    time(

    module load jdk/8u181
    module load gcc/6.4.0 samtools/1.9
    module load bcftools/1.9

    {config[tools][varscan]} somatic \
      {input.nrm} \
      <(samtools mpileup -q 1 -f {input.ref} {input.tum} ) \
      {params.pfx} \
      --output-vcf

    # rename samples in VCF header (VarScan reports as NORMAL, TUMOUR)
    bcftools reheader -s <(printf "%s\n%s\n" RN {wildcards.sample}) \
      {output.snp} > {output.vcf}
    #rm {output.vcf}.idx
    ) 1>{log} 2>&1
    """

rule varscan_post:
  input:  "varscan/{sample}.raw.vcf"
  output: "varscan/{sample}.vcf"
  threads: 1
  shell:
    """
    awk -f {config[scripts]}/varscan_vcf_fix.awk {input} > {output}
    """
