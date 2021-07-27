#!/bin/bash
# Calculate run times from job scripts
for log in slurm-*.out; do
  rep=$(perl -ne '/Running clonipy for replicate sims\/(.+$)/ && print "$1"' $log)
  t_sim=$($HOME/scripts/datediff.py --format="%a_%b_%d_%H:%M:%S_%Z_%Y" $(grep "CEST" $log | tr " " "_" | tr "\n" " "))
  t_bwa=$()
  printf "%s\t%s\t%s\n" $log $rep $t_sim
done
