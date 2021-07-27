#!/bin/bash
if [[ $# -lt 1 ]]; then
  echo
  echo "usage: $0 file.vcf [AWK condition]"
  echo
  exit 1
fi
fn=$1
cond=${2:-1}

echo "fn: $fn"
echo "cond: $cond"
echo "---"

cat $fn | awk -v FS="\t" '/^#C/{split($0,h,"\t")}!/^#/&&'${cond}'{for(i=1;i<=NF;i++)print h[i],$i}'
