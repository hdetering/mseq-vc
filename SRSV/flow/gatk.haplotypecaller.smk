# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule haplotypecaller:
    input:
        ref="ref.fa",
        ref_fai="ref.fa.fai",
        ref_dict="ref.dict",
        tumor="{sam}.bam",
        normal="RN.bam"
    output:
        "haplotypecaller/{sam}.g.vcf"
    params:
        workdir="haplotypecaller"
    log:
        "log/{sam}.haplotypecaller.log"
    threads: 1
    shell:
        """
        time(
        module load gatk/4.0.10.0
        if [[ ! -d "{params.workdir}" ]]; then mkdir -p {params.workdir}; fi

        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.tumor} \
            -O {output} \
            -ERC GVCF
        ) 1>{log} 2>&1
        """

rule genomicsdbimport:
    input:
        "haplotypecaller/R1.g.vcf",
        "haplotypecaller/R2.g.vcf",
        "haplotypecaller/R3.g.vcf",
        "haplotypecaller/R4.g.vcf",
        "haplotypecaller/R5.g.vcf",
        "haplotypecaller/RN.g.vcf"
    output:
        "haplotypecaller/db"
    log:
        "log/haplotypecaller.genomicsdbimport.log"
    threads: 1
    shell:
        """
        time(
        module load gatk/4.0.10.0
        gatk GenomicsDBImport \
          -V {input[0]} \
          -V {input[1]} \
          -V {input[2]} \
          -V {input[3]} \
          -V {input[4]} \
          -V {input[5]} \
          --genomicsdb-workspace-path {output} \
          --intervals chr1
        ) >{log} 2>&1
        """

rule genotypegvcfs:
    input:
        ref="ref.fa",
        db="haplotypecaller/db"
    output:
        "haplotypecaller.raw.vcf"
    log:
        "log/haplotypecaller.genotypegvcfs.log"
    threads: 1
    shell:
        """
        time(
        module load gatk/4.0.10.0
        gatk GenotypeGVCFs \
          -R {input.ref} \
          -V gendb://{input.db} \
          -G StandardAnnotation \
          -O {output}
        ) >{log} 2>&1
        """

rule selectvariants:
    input:
        ref="ref.fa",
        raw="haplotypecaller.raw.vcf"
    output:
        "haplotypecaller.filt.vcf"
    log:
        "log/haplotypecaller.selectvariants.log"
    threads: 1
    shell:
        """
        time(
        module load gatk/4.0.10.0
        gatk SelectVariants \
          -R {input.ref} \
          -V {input.raw} \
          -select 'vc.getGenotype("RN").isHomRef()' \
          -O {output}
        ) >{log} 2>&1
        """
