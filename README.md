# Benchmarking somatic variant calling using multi-regional sequencing data

This repository contains workflow and analysis scripts pertaining to a simulation study that aims to assess the performance of several popular somatic variant calling methods.

The study was carried out using two different simulation approaches:

## De-novo - Simulated Reads, Simulated Variants (full simulation)

1. Reads were simulated from an artificial reference genome. 
2. Germline and somatic variants were simulated and spiked into the reads.

## Spike-in - Real Reads, Simulated Variants (partial simulation)

1. An empirical sequencing data set was split into multiple samples.
2. Reads were sorted by germline allele (maternal/paternal allele) using germline variants.
3. Somatic variants were simulated and spiked into the reads.

# Structure of this repository

- R code to generate the plots for the publictation are in the scripts [de-novo.R](de-novo.R) and [spike-in.R](spike-in.R)
- Scripts that show how variant callers where run are in the respective directories under [workflow](workflow).
- Technical details about the simulations can be gleaned from the content of the [sim](sim) directories.
