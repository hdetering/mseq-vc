# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#SAMPLES = ['R1', 'R2', 'R3', 'R4', 'R5']
#MODEL="/mnt/lustre/scratch/home/uvi/be/hde/m-seq_varcall/src/SNooPer_v0.03_trained_models/SNooPer_v0.03_trained_model_45-55X/model/"

def get_snooper_model(wildcards):
  model_30x  = "/mnt/netapp2/posadalab2/uvibehde/mseq-vc/de-novo/src/SNooPer_v0.03_trained_models/SNooPer_v0.03_trained_model_25-35X/model/"
  model_50x  = "/mnt/netapp2/posadalab2/uvibehde/mseq-vc/de-novo/src/SNooPer_v0.03_trained_models/SNooPer_v0.03_trained_model_45-55X/model/"
  model_100x = "/mnt/netapp2/posadalab2/uvibehde/mseq-vc/de-novo/src/SNooPer_v0.03_trained_models/SNooPer_v0.03_trained_model_85-95X/model/"
  model_300x = "/mnt/netapp2/posadalab2/uvibehde/mseq-vc/de-novo/src/SNooPer_v0.03_trained_models/SNooPer_v0.03_trained_model_85-95X/model/"
  if ("depth" in config):
    if (config["depth"]=="30"):
      return model_30x
    elif (config["depth"]=="50"):
      return model_50x
    elif (config["depth"]=="100"):
      return model_100x
    elif (config["depth"]=="300"):
      return model_300x
  else:
    return model_50x

rule snooper:
    input:
        tumour="snooper/{sample}.pu",
        normal="snooper/RN.pu"
        #ref="data/{replicate}/ref.fa",
        #fai="data/{replicate}/ref.fa.fai",
        #nrm="data/{replicate}/RN.bam",
        #tum=["data/{replicate}/%s.bam" % x for x in SAMPLES]
    output: "snooper/{sample}.vcf"
    log:    "log/{sample}.snooper.log"
    params:
        workdir="snooper/{sample}",
        model=get_snooper_model
    threads: 1
    shell:
        """
        time (
        module load gcc/6.4.0 snooper/0.03
        module load bcftools/1.9

        # create workdir for sample
        if [[ -d {params.workdir} ]]; then
          rm -r {params.workdir}
        fi
        mkdir -p {params.workdir}
        cd {params.workdir}
        ln -s ../{wildcards.sample}.pu tset_T_{wildcards.sample}.pu
        ln -s ../RN.pu tset_N_{wildcards.sample}.pu 
        cd -

        SNooPer.pl \
            -i {params.workdir} \
            -o {params.workdir} \
            -a1 somatic \
            -a2 classify \
            -m {config[model]} \
            -w $EBROOTWEKA/weka.jar \
            -id 0

        # rename samples in VCF header (SNooPer reports as SOMATIC)
        bcftools reheader \
          -s <(printf "%s\n" {wildcards.sample}) \
          {params.workdir}/*/SNooPer_0_classification.snp \
        > {output}
        #| gunzip -c > {output}
 
        # cleanup
        rm -r {params.workdir}

        ) >{log} 2>&1
        """

rule snooper_pileup:
    input:
        bam="{sample}.bam",
        ref="ref.fa"
    output:
        "snooper/{sample}.pu"
    log:
        "log/{sample}.snooper_pileup.log"
    threads: 1
    shell:
        """
        time (
        module load gcc/6.4.0 samtools/1.8
        samtools mpileup \
          --fasta-ref {input.ref} \
          --output-MQ \
          --output {output} \
          {input.bam}
        ) >{log} 2>&1 
        """
