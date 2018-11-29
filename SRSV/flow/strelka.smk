# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule strelka2:
    input:
        ref="ref.fa",
        ref_idx="ref.fa.fai",
        tumor="{sample}.bam",
        normal="RN.bam"
    output:
        "strelka2/{sample}.vcf"
    log:
        "log/{sample}.strelka2.log"
    params:
        workdir="strelka2/{sample}",
        result="strelka2/{sample}/results/variants/somatic.snvs.vcf.gz"
    threads: 10
    shell:
        """
        time (
        module load bcftools/1.9

        {config[tools][strelka2_dir]}/bin/configureStrelkaSomaticWorkflow.py \
          --normalBam {input.normal} \
          --tumorBam {input.tumor} \
          --ref {input.ref} \
          --runDir {params.workdir}

        {params.workdir}/runWorkflow.py -m local -j {threads}

        # rename samples in VCF header (Strelka reports as NORMAL, TUMOR)
        bcftools reheader \
          -s <(printf "%s\n%s\n" RN {wildcards.sample}) \
          {params.result} \
        | gunzip -c > {output}
        ) 1>{log} 2>&1
        """

rule strelka1:
    input:
        ref="ref.fa",
        ref_idx="ref.fa.fai",
        tumor="{sample}.bam",
        normal="RN.bam"
    output:
        "strelka1/{sample}.vcf"
    log:
        "log/{sample}.strelka1.log"
    params:
        rootdir="strelka1",
        workdir="strelka1/{sample}",
        result="strelka1/{sample}/results/passed.somatic.snvs.vcf"
    threads: 10
    shell:
        """
        time (
        #module load bcftools/1.9
        module load gcccore/6.4.0 bcftools/1.9

        # create working direktory if necessary
        if [[ ! -d {params.rootdir} ]]; then
            mkdir -p {params.rootdir}
        fi

        # configure
        {config[tools][strelka1_dir]}/bin/configureStrelkaWorkflow.pl \
            --normal=./{input.normal} \
            --tumor=./{input.tumor} \
            --ref=./{input.ref} \
            --config={config[tools][strelka1_conf]} \
            --output-dir={params.workdir}

        # run analysis
        make -j {threads} -C {params.workdir}

        # rename samples in VCF header (Strelka reports as NORMAL, TUMOR)
        bcftools reheader \
          -s <(printf "%s\n%s\n" RN {wildcards.sample}) \
          {params.result} \
        > {output}
 
        ) 1>{log} 2>&1
        """
