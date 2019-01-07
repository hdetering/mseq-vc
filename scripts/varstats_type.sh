#!/bin/bash
DATA=$1
SIM=$LUSTRE/m-seq_varcall/sims

module load gcccore/6.4.0 bedtools/2.27.1

printf "replicate\tcaller\ttype\tchrom\tpos\n"
for rep in $(ls $DATA); do
  # export true somatic variants to BED file
  vcf_som="$SIM/$rep/output/somatic.vcf"
  bed_som="$DATA/$rep/somatic.true.bed"
  cat $vcf_som | awk '!/^#/{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed_som
  # true number of somatic variants
  nT=$(wc -l ${bed_som} | cut -d" " -f1)

  # MuTect1 stats
  vcf="$DATA/$rep/mutect1.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '$7=="PASS"{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tMuTect1\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tMuTect1\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tMuTect1\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # Mutect2 stats
  vcf="$DATA/$rep/mutect2.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '$7=="PASS"{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tMutect2\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tMutect2\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tMutect2\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # VarDict stats (LikelySomatic + StrongSomatic)
  vcf="$DATA/$rep/vardict.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '$7=="PASS" && /STATUS=[^;]*Somatic/{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tVarDict\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tVarDict\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tVarDict\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # VarDict_mfilt stats (filtered across regions)
  vcf="$DATA/$rep/vardict.mfilt.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '$7=="PASS" && /STATUS=[^;]*Somatic/{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tVarDict_mfilt\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tVarDict_mfilt\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tVarDict_filt\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # Strelka2 stats
  vcf="$DATA/$rep/strelka2.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '$7=="PASS"{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tStrelka2\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tStrelka2\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tStrelka2\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # CaVEMan stats
  vcf="$DATA/$rep/caveman.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '!/^#/{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tCaVEMan\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tCaVEMan\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tCaVEMan\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # SNooPer stats
  vcf="$DATA/$rep/snooper.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '$7=="PASS"{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tSNooPer\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tSNooPer\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tSNooPer\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # SomaticSniper stats
  vcf="$DATA/$rep/somaticsniper.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '!/^#/{printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tSomaticSniper\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tSomaticSniper\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tSomaticSniper\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

  # MultiSNV stats
  vcf="$DATA/$rep/multisnv.vcf"
  bed=${vcf%.vcf}.filt.bed
  if [[ -f $vcf ]]; then
    cat $vcf | awk '$7=="PASS"{s=0;for(i=11;i<NF;i++){split($i,a,":");if(a[4]==2)s=1}if(s)printf "%s\t%s\t%s\n",$1,$2,$2}' > $bed
    # True Positives
    bedtools intersect -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tMultiSNV\tTP\t%s\t%s\n",r,$1,$2}'
    # False Positives
    bedtools subtract -b ${bed_som} -a ${bed} \
    | awk -v r=$rep '{printf "%s\tMultiSNV\tFP\t%s\t%s\n",r,$1,$2}'
    # False Negatives
    bedtools subtract -a ${bed_som} -b ${bed} \
    | awk -v r=$rep '{printf "%s\tMultiSNV\tFN\t%s\t%s\n",r,$1,$2}' 
  fi

done
