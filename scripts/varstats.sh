#!/bin/bash
DATA=$1
SIM=$LUSTRE/m-seq_varcall/sims

printf "replicate\tcaller\tn_T\tn_P\tn_TP\n"
for rep in $(ls $DATA); do
  vcf_som="$SIM/$rep/output/somatic.vcf"
  # true number of somatic variants
  nT=$(grep -vcP "^#" ${vcf_som})

  # MuTect1 stats
  nP_mutect1=0
  nTP_mutect1=0
  vcf="$DATA/$rep/mutect1.vcf"
  if [[ -f $vcf ]]; then
    nP_mutect1=$(cat $vcf | awk '$7=="PASS"' | wc -l)
    nTP_mutect1=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '$7=="PASS"{print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "MuTect1" $nT ${nP_mutect1} ${nTP_mutect1}

  # MuTect2 stats
  nP_mutect2=0
  nTP_mutect2=0
  vcf="$DATA/$rep/mutect2.vcf"
  if [[ -f $vcf ]]; then
    nP_mutect2=$(cat $vcf | awk '$7=="PASS"' | wc -l)
    nTP_mutect2=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '$7=="PASS"{print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "Mutect2" $nT ${nP_mutect2} ${nTP_mutect2}

  # VarDict
  nP_vardict=0
  nTP_vardict=0
  vcf="$DATA/$rep/vardict.vcf"
  if [[ -f $vcf ]]; then
    nP_vardict=$(cat $vcf | awk '!/^#/ && $7=="PASS" && /STATUS=[^;]*Somatic/' | wc -l)
    nTP_vardict=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '!/^#/ && $7=="PASS" && /STATUS=[^;]*Somatic/{print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "VarDict" $nT ${nP_vardict} ${nTP_vardict}

  # VarDict filtered across regions
  nP_vardict=0
  nTP_vardict=0
  vcf="$DATA/$rep/vardict.mfilt.vcf"
  if [[ -f $vcf ]]; then
    nP_vardict=$(cat $vcf | awk '!/^#/ && $7=="PASS" && /STATUS=[^;]*Somatic/' | wc -l)
    nTP_vardict=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '!/^#/ && $7=="PASS" && /STATUS=[^;]*Somatic/{print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "VarDict_mfilt" $nT ${nP_vardict} ${nTP_vardict}

  # Strelka1 stats
  nP_strelka=0
  nTP_strelka=0
  vcf="$DATA/$rep/strelka1.vcf"
  if [[ -f $vcf ]]; then
    nP_strelka=$(awk '$7=="PASS"' $vcf | wc -l)
    nTP_strelka=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '$7=="PASS"{print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "Strelka1" $nT ${nP_strelka} ${nTP_strelka}

  # Strelka2 stats
  nP_strelka=0
  nTP_strelka=0
  vcf="$DATA/$rep/strelka2.vcf"
  if [[ -f $vcf ]]; then
    nP_strelka=$(awk '$7=="PASS"' $vcf | wc -l)
    nTP_strelka=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '$7=="PASS"{print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "Strelka2" $nT ${nP_strelka} ${nTP_strelka}

  # CaVEMan stats
  nP_caveman=0
  nTP_caveman=0
  vcf="$DATA/$rep/caveman.vcf"
  if [[ -f $vcf ]]; then
    nP_caveman=$(grep -cvP "^#" $vcf)
    nTP_caveman=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(grep -vP "^#" ${vcf} | cut -f2 | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "CaVEMan" $nT ${nP_caveman} ${nTP_caveman}

  # SNooPer stats
  nP_snooper=0
  nTP_snooper=0
  vcf="$DATA/$rep/snooper.vcf"
  if [[ -f $vcf ]]; then
    nP_snooper=$(awk '$7=="PASS"' $vcf | wc -l)
    nTP_snooper=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '$7=="PASS"{print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "SNooPer" $nT ${nP_snooper} ${nTP_snooper}

  # SomaticSniper stats
  nP_sniper=0
  nTP_sniper=0
  vcf="$DATA/$rep/somaticsniper.vcf"
  if [[ -f $vcf ]]; then
    nP_sniper=$(grep -cvP "^#" $vcf)
    nTP_sniper=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(grep -vP "^#" ${vcf} | cut -f2 | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "SomaticSniper" $nT ${nP_sniper} ${nTP_sniper}

  # MultiSNV stats
  nP_multisnv=0
  nTP_multisnv=0
  vcf="$DATA/$rep/multisnv.vcf"
  if [[ -f $vcf ]]; then
    nP_multisnv=$(cat $vcf | awk '$7=="PASS"{s=0;for(i=11;i<=NF;i++){split($i,a,":");if(a[4]==2)s=1};if(s)print}' | wc -l)
    nTP_multisnv=$(join \
      <(grep -vP "^#" ${vcf_som} | cut -f2 | sort) \
      <(cat $vcf | awk '$7=="PASS"{s=0;for(i=11;i<=NF;i++){split($i,a,":");if(a[4]==2)s=1};if(s)print $2}' | sort) \
      | wc -l
    )
  fi
  printf "%s\t%s\t%d\t%d\t%d\n" $rep "MultiSNV" $nT ${nP_multisnv} ${nTP_multisnv}

done
