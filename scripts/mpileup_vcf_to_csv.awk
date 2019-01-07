#!/bin/awk
# vim: tabstop=2 expandtab
# coding: utf-8

BEGIN {
  OFS=","
  print "chrom","pos","id_sample","A","C","G","T"
}

# store headers
/^#CHROM/ {
  split($0, hdr, "\t")
}

# process variant lines
!/^#/ {
   # determine position of AD tag within FORMAT field
  AD = 1
  split($9, f, ":")
  while (f[AD] != "AD") AD++ # this will cause an error if there is no "AD" tag

  # split ALT alleles
  split($5, a, ",")
  n = length(a) # number of detected alleles (REF+ALT)

  # analyse read counts for each sample
  for (i=10; i<=NF; i++) {
    # initialize read counts for alleles
    c["A"] = 0
    c["C"] = 0
    c["G"] = 0
    c["T"] = 0
    # split format tags
    split($i, tags, ":")
    # split read counts
    split(tags[AD], ad, ",")
    # set read count for REF allele
    c[$4] = ad[1]
    # set read counts for ALT alleles
    for (j=2; j<=n; j++) {
      c[a[j-1]] = ad[j]
    }
    # print allele depths for sample
    print $1,$2,hdr[i],c["A"],c["C"],c["G"],c["T"] 
  }
}
