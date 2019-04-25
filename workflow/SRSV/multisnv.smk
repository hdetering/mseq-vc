# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from glob import glob
#BAMS_TUM = glob_wildcards("data/%s/{sample,R\d+\.bam}" % REPS)
SAMPLES = ['R1', 'R2', 'R3', 'R4', 'R5']

rule multisnv:
    input:
        ref="ref.fa",
        fai="ref.fa.fai",
        nrm="RN.bam",
        tum=["%s.bam" % x for x in SAMPLES]
    output: "multisnv.vcf"
    log:    "log/multisnv.log"
    params:
        #dir="data/{replicate}"
    threads: 1
    shell:
        """
        time (
        module purge
        module load gcc/6.4.0 multisnv/2.3
        
        BAM=({input.tum})
        echo "${{BAM[@]}}"
        NBAM=$(( ${{#BAM[@]}}+1 ))
        echo "$NBAM tumor BAMs"

        # NOTE: adapt dmax, medianN, medianT according to replicate scenario
        echo "Running MultiSNV with median depth {config[depth]}"

        multiSNV \
            --number $NBAM \
            --fasta {input.ref} \
            --bam {input.nrm} ${{BAM[@]}}\
            --fout {output} \
            --mu 0.000001 \
            --minBase 20 \
            --minMapQual 30 \
            --dmin 5 \
            --low_depth 6 \
            --dmax 500 \
            --mva 1 \
            --weak_evidence 0.03 \
            --normal_contamination 0.03 \
            --minVariantReadsForTriallelicSite 2 \
            --flag-homopolymer 5 \
            --medianN {config[depth]} \
            --medianT {config[depth]} \
            --Rmdups 0 \
            --include-germline 1 \
            --include-LOH 1 \
            --print 1 \
            --conv 1
        ) >{log} 2>&1
        """
