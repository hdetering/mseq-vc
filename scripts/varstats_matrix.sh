#!/bin/bash
# vim: syntax=sh tabstop=2 expandtab
# coding: utf-8

SIMS=$LUSTRE/m-seq_varcall/sims
DATA=$LUSTRE/m-seq_varcall/data

module load bcftools/1.9
# required for consistent sorting
export LC_ALL=C

printf "replicate,sample,caller,chrom,pos,status\n"
for rep in $(ls $SIMS); do
  som="$DATA/$rep/var.som.sites.tsv"
  awk -v OFS="\t" '!/^#/{print $1,$2}' $SIMS/$rep/output/somatic.vcf > $som
  
  for caller in caveman mutect1 mutect2 snooper somaticsniper strelka2 vardict; do
    for region in R1 R2 R3 R4 R5; do
      case "$caller" in
        caveman) ;&
        somaticsniper)
          join -a 1 <(awk '!/^#/{print $1"_"$2}' $SIMS/$rep/output/somatic.vcf | sort) \
                    <(awk -v OFS="\t" '!/^#/{print $1"_"$2,1}' $DATA/$rep/$caller/${region}.vcf | sort) |
          while read chrompos status; do
            if [ -z "$status" ]; then status=0; fi
            printf "%s,%s,%s,%s,%d\n" $rep $region $caller ${chrompos/_/,} $status
          done
        ;;
        mutect1) ;&
        mutect2) ;&
        snooper) ;&
        strelka2)
          join -a 1 <(awk '!/^#/{print $1"_"$2}' $SIMS/$rep/output/somatic.vcf | sort) \
                    <(awk -v OFS="\t" '!/^#/{s=-1;if($7=="PASS")s=1;print $1"_"$2,s}' $DATA/$rep/$caller/${region}.vcf | sort) |
          while read chrompos status; do
            if [ -z "$status" ]; then status=0; fi
            printf "%s,%s,%s,%s,%d\n" $rep $region $caller ${chrompos/_/,} $status
          done
        ;;
        vardict)
          join -a 1 <(awk '!/^#/{print $1"_"$2}' $SIMS/$rep/output/somatic.vcf | sort) \
                    <(awk -v OFS="\t" '!/^#/{s=-1;if($7=="PASS" && $8~/STATUS=[^;]*Somatic/)s=1;print $1"_"$2,s}' $DATA/$rep/$caller/${region}.vcf | sort) |
          while read chrompos status; do
            if [ -z "$status" ]; then status=0; fi
            printf "%s,%s,%s,%s,%d\n" $rep $region $caller ${chrompos/_/,} $status
          done
        ;;
         *) continue ;;
      esac
    done
  done
done
