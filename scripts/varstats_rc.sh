#!/bin/bash
SIMS=$LUSTRE/m-seq_varcall/sims
DATA=$LUSTRE/m-seq_varcall/data

module load bcftools/1.9

printf "replicate\tsample\tcaller\tchrom\tpos\tstatus\n"
