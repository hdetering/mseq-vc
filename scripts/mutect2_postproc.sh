#!/bin/bash
module load jdk/1.8.0

# Filter Mutect2 calls
for rep in $(ls data); do 
  if [[ ! -f data/$rep/mutect2.filt.vcf ]]; then 
    for s in R1 R2 R3 R4 R5; do 
      /mnt/netapp1/posadalab/APPS/gatk-4.0.7.0/gatk FilterMutectCalls \
        -V data/$rep/mutect2/$s.vcf \
        -O data/$rep/mutect2/$s.filt.vcf 
    done 
  fi 
done

# Combine Variants
module load gatk/3.7

for rep in $(ls data); do 
  if [[ ! -f data/$rep/mutect2.filt.vcf ]]; then 
    java -jar $GATK -T CombineVariants \
      -R data/$rep/ref.fa \
      -V data/$rep/mutect2/R1.filt.vcf \
      -V data/$rep/mutect2/R2.filt.vcf \
      -V data/$rep/mutect2/R3.filt.vcf \
      -V data/$rep/mutect2/R4.filt.vcf \
      -V data/$rep/mutect2/R5.filt.vcf \
      -genotypeMergeOptions UNIQUIFY \
      -o data/$rep/mutect2.filt.vcf
  fi
done
