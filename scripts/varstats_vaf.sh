#!/bin/bash
DATA=$1
SIM=$LUSTRE/m-seq_varcall/sims
BIN=$LUSTRE/m-seq_varcall/scripts

vstats="varstats_type.tsv"
if [[ ! -f $vstats ]]; then
  echo "[ERROR]: missing file '$vstats'. (run scripts/varstats.type.vcf)"
  exit 1
fi

module purge
module load bcftools/1.9
module load miniconda3

printf "replicate\tcaller\tchrom\tpos\tmax_vaf\n"
for rep in $(ls $DATA); do
  
  # check if there are variant calls to process
  if [[ ! -f $DATA/$rep/mutect2.vcf ]]; then
    continue
  fi

  bed="$DATA/$rep/var.sites.bed"
  vcf="$DATA/$rep/var.sites.vcf"
  
  # extract sites of all simulated and called variants
  grep -P "^${rep}" $vstats | awk '{printf "%s\t%d\t%d\n",$4,$5-1,$5}' | sort -u > $bed
  # calculate read counts for given sites
  bcftools mpileup \
    --output-type v \
    --fasta-ref $DATA/$rep/ref.fa \
    --regions-file $bed \
    --skip-indels \
    --annotate DP,AD \
    --output $vcf \
    $DATA/$rep/R{1,2,3,4,5}.bam

  if [[ -f $vcf ]]; then
    # print VAF
    python3 scripts/varstats_vaf.py $vcf $DATA/$rep/mutect2.vcf $SIM/$rep/output/somatic.vcf "$rep,Mutect2" 
  fi

  #rm $bed $vcf 
done
