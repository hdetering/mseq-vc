#!/bin/awk
# INPUT:
#  1. VCF file, samples in columns 10-14 will be considered
# PARAMS:
#  r : replicate ID
!/^#/ { 
  n=0
  f=0
  # loop over samples
  for(i=10; i<15; i++) {
    split($i,a,":")
      if(a[1]=="0/1"||a[1]=="1/1") {
        n++
        split(a[4], b, ",")
        if(b[2]/a[3] > f) {
          f=b[2]/a[3]
      }
    }
  }
  printf "%s\t%s\t%s\t%d\t%0.2f\n", r, $1, $2, n, f
}
