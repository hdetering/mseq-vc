# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
NORMAL_PU = "/mnt/netapp1/posadalab/Tama_Laura/TestingVC/ChineseSon/workdir/ChineseSon.H.chr21.bam.q10.pileup"
REF = "/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa"

rule samtools_mpileup:
  input:
    ref=REF,
    bam="{sample}.bam"
  output:
    pup="varscan/{sample}.pu"
  log: "log/{sample}.mpileup.log"
  threads: 1
  shell:
    """
    module load gcc/6.4.0 samtools/1.9

    time (
    samtools mpileup \
      --fasta-ref {input.ref} \
      --min-MQ 10 \
      {input.bam} \
      > {output.pup}
    ) >{log} 2>&1
    """ 

rule varscan:
  input:
    ref=REF,
    tum="varscan/{sample}.pu",
    nrm=NORMAL_PU
  output:
    vcf="varscan/{sample}.raw.vcf",
    snp="varscan/{sample}.snp.vcf",
    ind="varscan/{sample}.indel.vcf"
  params:
    pfx="varscan/{sample}"
  log: "log/{sample}.varscan.log"
  threads: 1
  shell:
    """
    time(

    module load jdk/8u181
    module load gcc/6.4.0 samtools/1.9
    module load bcftools/1.9

    {config[tools][varscan]} somatic \
      {input.nrm} \
      {input.tum} \
      {params.pfx} \
      --output-vcf

    # rename samples in VCF header (VarScan reports as NORMAL, TUMOUR)
    bcftools reheader -s <(printf "%s\n%s\n" RN {wildcards.sample}) \
      {output.snp} > {output.vcf}
    ) 1>{log} 2>&1
    """

rule varscan_post:
  input:  "varscan/{sample}.raw.vcf"
  output: "{sample}.VarScan.vcf"
  threads: 1
  shell:
    """
    awk -f {config[scripts]}/varscan_vcf_fix.awk {input} > {output}
    """
