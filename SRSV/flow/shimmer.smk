# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule shimmer:
    input:
        ref="ref.fa",
        ref_idx="ref.fa.fai",
        tumor="{sample}.bam",
        normal="RN.bam"
    output:
        "shimmer/{sample}.vcf"
    log:
        "log/{sample}.muse_call.log"
    params:
        workdir="shimmer",
        out_dir="shimmer/{sample}"
    threads: 1
    shell:
        """
        time (
        export PATH=/mnt/netapp1/posadalab/APPS/Shimmer/bin:$PATH
        module load gcc/6.4.0 samtools/1.9 R/3.5.1

        # create output dir if necessary
        if [[ ! -d "{params.workdir}" ]]; then
          mkdir -p {params.workdir}
        fi

        {config[tools][shimmer]} \
          {input.normal} \
          {input.tumor} \
          --ref {input.ref} \
          --outdir {params.out_dir}
        
        # rename samples in VCF header (reported as NORMAL, TUMOR)
        bcftools reheader \
          -s <(printf "%s\n%s\n" RN {wildcards.sample}) \
          {params.out_dir}/somatic_diffs.vcf \
        > {output}

       ) 1>{log} 2>&1
       """
