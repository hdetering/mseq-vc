# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule jsm_train:
    input:
        ref="ref.fa",
        ref_idx="ref.fa.fai",
        tumor="{sample}.bam",
        normal="RN.bam"
    output:
        "jointsnvmix/{sample}.jsm_train.params"
    log:
        "log/{sample}.jsm_train.log"
    params:
        workdir="jointsnvmix"
    threads: 10
    shell:
        """
        time (
        module load bcftools/1.9
        source activate $POSADALAB/APPS/JointSNVMix-0.8/env

        # create output dir if necessary
        if [[ ! -d "{params.workdir}" ]]; then
          mkdir -p {params.workdir}
        fi

        jsm.py train \
          --model beta_binomial \
          {input.ref} \
          {input.normal} \
          {input.tumor} \
          {output}

        ) 1>{log} 2>&1
        """

rule jsm_classify:
    input:
        ref="ref.fa",
        ref_idx="ref.fa.fai",
        tumor="{sample}.bam",
        normal="RN.bam",
        params="jointsnvmix/{sample}.jsm_train.params"
    output:
        "jointsnvmix/{sample}.jsm_classify.tsv"
    log:
        "log/{sample}.jsm_classify.log"
    params:
    threads: 10
    shell:
        """
        time (
        source activate $POSADALAB/APPS/JointSNVMix-0.8/env
        
        jsm.py classify \
          --parameters_file {input.params} \
          --model beta_binomial \
          --out_file {output} \
          {input.ref} \
          {input.normal} \
          {input.tumor}
        ) 1>{log} 2>&1
        """
 
rule jsm_tsv_to_vcf:
    input:
        tsv="jointsnvmix/{sample}.jsm_classify.tsv",
    output:
        "jointsnvmix/{sample}.vcf"
    log:
        "log/{sample}.jsm_tsv2vcf.log"
    threads: 1
    shell:
        """
        time (
        python {config[scripts]}/jointsnvmix2vcf.py \
          --min_p 0.5 \
          --id_sample {wildcards.sample} \
          {input.tsv} \
        > {output}
        ) 1>{log} 2>&1
        """
