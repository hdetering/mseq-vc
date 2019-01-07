# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
from __future__ import print_function, division
import os, sys
import argparse

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('tsv', help='Output of jsm.py classify.')
  parser.add_argument('--min_p', type=float, default=0.0, help='Minimum probability (p_AA_AB + p_AA_BB) filter.')
  parser.add_argument('--id_sample', default='TUMOR', help='Sample ID to include in VCF.')
  args = parser.parse_args()

  # sanity checks
  if not os.path.exists(args.tsv):
    print("[ERROR] File does not exist: %s" % args.tsv, file=sys.stderr)
    sys.exit(1)

  return args

def tsv2vcf(fn_tsv, min_p, id_sample):
  # write VCF header
  print('''##fileformat=VCFv4.0
##source=JointSNVMix
##reference=ref.fa
##INFO=<ID=NR,Number=1,Type=Integer,Description="Number of reads supporting reference in normal">
##INFO=<ID=NV,Number=1,Type=Integer,Description="Number of reads supporting variant in normal">
##INFO=<ID=TR,Number=1,Type=Integer,Description="Number of reads supporting reference in tumor">
##INFO=<ID=TV,Number=1,Type=Integer,Description="Number of reads supporting variant in tumor">
##INFO=<ID=RPS,Number=.,Type=Float,Description="Raw probability of somatic mutation P=p_AA_AB+p_AA_BB">
##INFO=<ID=PS,Number=.,Type=Float,Description="Post processed probability of somatic mutation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Maximum-likelihood somatic genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'''.format(id_sample))

  with open(fn_tsv, 'rt') as tsv:
    # skip header line
    line = tsv.readline()

    # parse data lines
    for line in tsv:
      fields = line.strip().split('\t')
      chr = fields[0]
      pos = fields[1]
      ref = fields[2]
      var = fields[3]
      normal_a = fields[4]
      normal_b = fields[5]
      tumor_a = fields[6]
      tumor_b = fields[7]
      p_AA_AA = float(fields[8])
      p_AA_AB = float(fields[9])
      p_AA_BB = float(fields[10])
      p_AB_AA = float(fields[11])
      p_AB_AB = float(fields[12])
      p_AB_BB = float(fields[13])
      p_BB_AA = float(fields[14])
      p_BB_AB = float(fields[15])
      p_BB_BB = float(fields[16])

      p_somatic = p_AA_AB + p_AA_BB

      # output variant only if minimum probability filter is a pass
      if (p_somatic < min_p):
        continue

      # choose max-likelihood somatic genotype
      gt = '0/1' if p_AA_AB >= p_AA_BB else '!/!'

      print("%s\t%s\t.\t%s\t%s\t%f\t.\tNR=%s;NV=%s;TR=%s;TV=%s;RPS=%f\tGT\t%s" % (
        chr, 
        pos, 
        ref, 
        var, 
        p_somatic,
        normal_a,
        normal_b,
        tumor_a,
        tumor_b,
        p_somatic,
        gt
      ))

if __name__ == '__main__':
  args = parse_args()
  tsv2vcf(args.tsv, args.min_p, args.id_sample)
