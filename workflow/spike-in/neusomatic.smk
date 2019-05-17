# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
NORMAL_PU = "/mnt/netapp1/posadalab/Tama_Laura/TestingVC/ChineseSon/workdir/ChineseSon.H.chr21.bam.q10.pileup"
REF_FA = "/mnt/netapp1/posadalab/phylocancer/RESOURCES/hs37d5.fa"
REF_BED = "/mnt/netapp1/posadalab/Tama_Laura/TestingVC/ChineseSon/workdir/chr21.forbamreadcount.bed"

rule neusomatic_prep:
  input:
    ref=REF_FA,
    bed=REF_BED,
    tum="{sample}.bam",
    nrm="ChineseSon.H.chr21.bam"
  output: "neusomatic/{sample}/call/dataset"
  log: "log/neusomatic_prep.{sample}.log"
  params:
    workdir="neusomatic/{sample}/call"
  threads: 1
  shell:
    """
    time (

    module load gcc/6.4.0 samtools/1.9
    source activate {config[tools][neusomatic]}/env

    python {config[tools][neusomatic]}/neusomatic/python/preprocess.py \
      --mode call \
      --reference {input.ref} \
      --region_bed {input.bed}  \
      --tumor_bam {input.tum} \
      --normal_bam {input.nrm} \
      --work {params.workdir} \
      --min_mapq 10 \
      --num_threads {threads} \
      --scan_alignments_binary {config[tools][neusomatic]}/neusomatic/bin/scan_alignments
    ) >{log} 2>&1
    """

rule neusomatic_call:
  input:
    ref=REF_FA,
    dat="neusomatic/{sample}/call/dataset"
  output:
    vcf="neusomatic/{sample}/call/pred.vcf"
  params:
    workdir="neusomatic/{sample}/call"
  log: "log/neusomatic_call.{sample}.log"
  threads: 1
  shell:
    """
    time(

    module load gcc/6.4.0 samtools/1.9
    source activate {config[tools][neusomatic]}/env

    CUDA_VISIBLE_DEVICES= python {config[tools][neusomatic]}/neusomatic/python/call.py \
      --candidates_tsv {input.dat}/*/candidates*.tsv \
      --reference {input.ref} \
      --out {params.workdir} \
      --checkpoint {config[tools][neusomatic]}/neusomatic/models/NeuSomatic_v0.1.3_standalone_Dream3.pth \
      --num_threads {threads} \
      --batch_size 100

    ) 1>{log} 2>&1
    """

rule neusomatic_post:
  input:
    ref=REF_FA,
    tum="{sample}.bam",
    vcf="neusomatic/{sample}/call/pred.vcf"
  output: "neusomatic/{sample}.vcf"
  log: "log/neusomatic_post.{sample}.log"
  params:
    workdir="neusomatic/{sample}/call"
  threads: 1
  shell:
    """
    time (

    module load gcc/6.4.0 samtools/1.9
    source activate {config[tools][neusomatic]}/env

    python $POSADALAB/APPS/neusomatic/neusomatic/python/postprocess.py \
      --reference {input.ref} \
      --tumor_bam {input.tum} \
      --pred_vcf {params.workdir}/pred.vcf \
      --candidates_vcf {params.workdir}/work_tumor/filtered_candidates.vcf \
      --output_vcf {output} \
      --work {params.workdir}
    ) >{log} 2>&1
    """
