#!/bin/awk
# vim: syntax=awk tabstop=2 expandtab
# coding: utf-8

# usage: 
#  awk -f mutect2vcf.awk somatic.vcf MuClone-sample-results.tsv > muclone.vcf

# global settings
BEGIN {
  cutoff = 0.5; # probability threshold to accept a mutation as present
  gq_max = 99; # Genotype Quality cap
}
# store first seven fields for each variant from input VCF
NR==FNR && !/^#/ {
  var[$3] = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;
  next;
}
# read region names from MuClone result header; print VCF header
NR>FNR && FNR==1 {
  print "##fileformat=VCFv4.1";
  print "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples with Data\">";
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
  print "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">";
  printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (i=1; $i!~/cluster/; i++) {
    nreg++;
    printf "\t"$i;
  }
  printf "\n";
  next;
}
# extract mutation status for each sample from MuClone results
NR>FNR {
  printf var[$1];
  printf "\tNS="nreg; # INFO
  printf "\tGT:GQ"; # FORMAT
  for (i=2; i<=nreg+1; i++) {
    if ($i >= 1) {
      gt = "0/1";
      gq = gq_max;
    }
    else if ($i >= cutoff) {
      gt = "0/1";
      gq = (-log(1-$i)/log(10)*10);
    }
    else {
      gt = "0/0"
      gq = (-log($i)/log(10)*10);
    }
    if (gq > gq_max) gq = gq_max;
    printf "\t%s:%.0f", gt, gq;
  }
  printf "\n";
}
