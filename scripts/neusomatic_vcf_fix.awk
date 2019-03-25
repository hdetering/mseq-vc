#!/bin/awk
# vim: syntax=awk tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Reformat Neusomatic 0.2.0 VCF output for downstream analysis.
# Formatting steps performed:
#  1. Fix INFO fields:
#     - move 'SCORE' to FORMAT
#     - 'DP', 'RO', 'AO', 'AF' appear twice (INFO, FORMAT) -> remove from INFO
#     - Remove final ';' in INFO column
#  2. Rename 'SAMPLE' to actual sample name (user param 'sid')
# NOTE: Assuming a single-sample VCF (tumour).
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-03-18
#------------------------------------------------------------------------------

BEGIN {
  FS = "\t";
  OFS = "\t";
  if (!sid) sid = "TUMOR";
}
# replace SAMPLE id in heading
/^#CHROM/ {
  $10 = sid;
  print; next;
}
!/^#/ { # fix non-header lines
  # remove redundant INFO fields (which are also in FORMAT)
  gsub("DP=[0-9]+;", "", $8);
  gsub("RO=[0-9]+;", "", $8);
  gsub("AO=[0-9]+;", "", $8);
  gsub("AF=[^;]+;", "", $8);
  # remove trailing semicolon from INFO
  gsub(";$", "", $8);
  # add remaining INFO fields to FORMAT (should be only 'SCORE')
  $9 = $9":SCORE";
  $10 = $10":"$8;
  gsub("SCORE=", "", $10);
  # set INFO column to empty
  $8 = ".";
  print; next;
}
{ # print all other lines (headers etc.)
  print; 
}
