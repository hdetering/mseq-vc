#!/usr/bin/env Rscript
# vim: syntax=r tabstop=4 expandtab
# coding: utf-8

script_name = "run_shearwater.R"
usage <- function() {
  sprintf("usage: %s target_region bam_file1 [bam_file2 [...]]", script_name)
}

# COMMAND LINE PARAMS
#------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop(usage())
}
# parse region to analyse
chr_int = strsplit(args[1], ":")[[1]]
chr = chr_int[1]
start_end = strsplit(chr_int[2], "-")[[1]]
pos_start = as.numeric(start_end[1])
pos_end = as.numeric(start_end[2])

# parse BAM file names
files = args[2:length(args)]

# MAIN
#------------------------------------------------------------------------------
require(deepSNV)
regions = GRanges(seqnames = chr, ranges = IRanges(pos_start, pos_end))
counts = loadAllData(files, regions)
# calculate Bayes factors
bf = bbb(counts, model="OR")
# convert to VCF object
vcf = bf2Vcf(bf, counts, regions)
# write results to file
writeVcf(vcf, "shearwater.vcf")
