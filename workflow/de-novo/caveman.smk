# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule caveman_setup:
    input:
        ref="ref.fa",
        ref_fai="ref.fa.fai",
        ref_dict="ref.dict",
        tumor="{sample}.bam",
        normal="RN.bam"
    output:
        "caveman/{sample}/caveman.cfg.ini"
    params:
        dir="caveman/{sample}"
    log:
        "log/{sample}.caveman.out"
    threads: 1
    group:  "caveman"
    shell:
        """
        time (
        echo "*** CaVEMan setup ***"
        # create working directory if not exists
        if [[ -d {params.dir} ]]; then
          mkdir -p {params.dir}
          # create (empty) list of regions to exclude
          touch {params.dir}/ignore-regions.tsv
        fi

        module load caveman/1.13.8
        caveman setup \
            --tumour-bam {input.tumor} \
            --normal-bam {input.normal} \
            --reference-index {input.ref_fai} \
            --ignore-regions-file {params.dir}/ignore-regions.tsv \
            --config-file {output} \
            --split-file {params.dir}/splitList \
            --alg-bean-file {params.dir}/alg_bean \
            --results-folder {params.dir}/results
        ) 1>{log} 2>&1
       """

rule caveman_split:
    input:
        cfg="caveman/{sample}/caveman.cfg.ini"
    output:
        lst="caveman/{sample}/splitList"
    params:
        dir="caveman/{sample}"
    log:
        "log/{sample}.caveman.out"
    threads: 1
    group:  "caveman"
    shell:
        """
        time(
        echo "*** CaVEMan Split ***"
        module load caveman/1.13.8
        caveman split \
            --config-file {input.cfg} \
            --index 1

        # cleanup: merge split lists for chromosomes
        cat {params.dir}/splitList.* > {params.dir}/splitList

        # create sentinetl files for jobids (CaVEMan uses job IDs to address regions...)
        #for jid in $(seq $(wc -l {params.dir}/splitList | cut -d" " -f1)); do
        #  touch {params.dir}/segment.$jid.mstep.todo
        #  touch {params.dir}/segment.$jid.estep.todo
        #done
        ) 1>>{log} 2>&1
        """

rule caveman_mstep:
    input:
        lst="caveman/{sample}/splitList"
    output:
        "caveman/{sample}/caveman_mstep.done"
    log:
        "log/{sample}.caveman.out"
    params:
        cfg="caveman/{sample}/caveman.cfg.ini"
    threads: 10
    group:  "caveman"
    shell:
        """
        time(
        echo "*** CaVEMan MStep ***"
        module load caveman/1.13.8 parallel/20180922
        seq $(wc -l {input.lst} | cut -d" " -f1) | \
        parallel \
          caveman mstep \
            --config-file {params.cfg} \
            --index {{}}
        ) 1>>{log} 2>&1
        
        touch {output}
        """

rule caveman_merge:
    input:
        "caveman/{sample}/caveman_mstep.done" 
    output:
        covs="caveman/{sample}/covs_arr",
        probs="caveman/{sample}/probs_arr"
    log:
        "log/{sample}.caveman.log"
    params:
        cfg="caveman/{sample}/caveman.cfg.ini"
    threads: 1
    group:  "caveman"
    shell:
        """
        time(
        echo "*** CaVEMan Merge ***"
        module load caveman/1.13.8
        caveman merge \
            --config-file {params.cfg} \
            --covariate-file {output.covs} \
            --probabilities-file {output.probs}
        ) 1>>{log} 2>&1
        """

rule caveman_estep:
    input:
        covs="caveman/{sample}/covs_arr",
        probs="caveman/{sample}/probs_arr"
    output:
        "caveman/{sample}/caveman_estep.done"
    log:
        "log/{sample}.caveman.log"
    params:
        cfg="caveman/{sample}/caveman.cfg.ini",
        lst="caveman/{sample}/splitList"
    threads: 10
    group:  "caveman"
    shell:
        """
        # NOTE: default values for normal-contamination, tumour-copy-number have been chosen
        #       according to recommendations in the GitHub README:
        #       https://github.com/cancerit/CaVEMan#default-settings
 
        time(
        echo "*** CaVEMan EStep ***"
        module load caveman/1.13.8 parallel/20180922
        seq $(wc -l {params.lst} | cut -d" " -f1) | \
        parallel \
          caveman estep \
            --index {{}} \
            --config-file {params.cfg} \
            --cov-file {input.covs} \
            --prob-file {input.probs} \
            --normal-contamination 0.1 \
            --tumour-copy-number 5 \
            --species Homo_syntheticus \
            --species-assembly simulation
        ) 1>>{log} 2>&1
       
        touch {output}
        """
 
rule caveman:
    input:
        "caveman/{sample}/caveman_estep.done"
    output:
        muts="caveman/{sample}.muts.vcf",
        snps="caveman/{sample}.snps.vcf",
        vcf="caveman/{sample}.vcf"
    log:
        out="log/{sample}.caveman.out"
    params:
        res="caveman/{sample}/results"
    threads: 1
    group:  "caveman"
    shell:
        """
        time (
        echo "*** Finalizing results ***"
        module load bcftools/1.9
        #source $POSADALAB/APPS/bcftools-1.9/prep_env.sh

        # create sorted list of partial VCF files
        for x in muts snps; do
          ls {params.res}/*/*.$x.vcf.gz \
          | awk 'BEGIN{{OFS="\t"}}{{split($1,a,"/");chr=a[4];fn=a[5];split(fn,b,"_");print chr,b[1],$1}}' \
          | sort -k2n | cut -f3 \
          > {params.res}/$x.vcf.gz.fofn
        done

        # concatenate partial VCFs
        bcftools concat --file-list {params.res}/muts.vcf.gz.fofn --output-type v > {params.res}/muts.vcf
        bcftools concat --file-list {params.res}/snps.vcf.gz.fofn --output-type v > {params.res}/snps.vcf

        # rename samples in VCF header (caveman reports as NORMAL, TUMOUR)
        bcftools reheader -s <(printf "%s\n%s\n" RN {wildcards.sample}) {params.res}/muts.vcf > {output.muts}
        bcftools reheader -s <(printf "%s\n%s\n" RN {wildcards.sample}) {params.res}/snps.vcf > {output.snps}

        # copy somatic mutations VCF to comply with naming convention
        cp {output.muts} {output.vcf}
        ) 1>>{log} 2>&1

        # compress intermediate files (reduce number of files / disk space)
        tar -zcf caveman/{wildcards.sample}.tar.gz caveman/{wildcards.sample} && rm -r caveman/{wildcards.sample}
        """
 
