# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule bwa_index:
    input:
        "{ref}"
    output:
        "{ref}.amb",
        "{ref}.ann",
        "{ref}.bwt",
        "{ref}.pac",
        "{ref}.sa"
    log:    "log/{ref}.bwa_index.log"
    group:  "bwa"
    shell:
        """
        time (
        module load gcc/6.4.0 bwa/0.7.17
        bwa index {input}
        ) >{log} 2>&1
        """

rule bwa_mem:
    input:
        ref="ref.fa",
        idx="ref.fa.bwt",
        sam="sim/{sample}.sam"
    output:
        bam="{sample}.bam"
        #stat="data/{replicate}/{sample}.stat"
    params:
        fastq="{sample}.fq.gz",
        sam="{sample}.sam"
    log:    "log/{sample}.bwa_mem.log"
    group:  "bwa"
    threads: 10
    shell:
        """
        time (
        module purge
        module load gcc/6.4.0 
        module load samtools/1.8
        module load bwa/0.7.17
        samtools fastq -n {input.sam} | gzip > {params.fastq}

        echo "Running BWA MEM on {threads} threads..."
        bwa mem \
            -t {threads} -p \
            -R '@RG\\tID:{wildcards.sample}\\tLB:ART\\tPL:ILLUMINA\\tSM:{wildcards.sample}\\tPU:simulation' \
            {input.ref} {params.fastq} \
        > {params.sam}

        samtools sort -@ {threads} {params.sam} > {output.bam} && rm {params.sam}
        samtools index {output.bam}
        #samtools flagstat {output.bam} > {{output.stat}}

        ) >{log} 2>&1
        """
