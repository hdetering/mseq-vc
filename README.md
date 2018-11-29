# Benchmarking somatic variant calling from multi-regional sequencing data

This repository contains workflow and analysis scripts pertaining to a simulation study that aims to assess the performance of several popular somatic variant calling methods.

The study was carried out using two different simulation approaches:

## SRSV - Simulated Reads, Simulated Variants (full simulation)

1. Reads were simulated from an artificial reference genome. 
2. Germline and somatic variants were simulated and spiked into the reads.

## RRSV - Real Reads, Simulated Variants (partial simulation)

1. An empirical sequencing data set was split into multiple samples.
2. Reads were sorted by germline allele (maternal/paternal allele) using germline variants.
3. Somatic variants were simulated and spiked into the reads.
