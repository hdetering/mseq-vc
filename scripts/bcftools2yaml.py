#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert Bcftools multi-sample VCF output to PyClone/MuClone YAML input.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-02-14
#------------------------------------------------------------------------------

import os, sys
import argparse
import vcf
import yaml

def parse_args():
  parser = argparse.ArgumentParser(description='Convert Bcftools multi-sample VCF output to PyClone YAML input.')
  parser.add_argument('infile', type=argparse.FileType(), help='Strelka2 VCF file.')
  parser.add_argument('outdir', help='Output directory for YAML files.')
  parser.add_argument('--samples', help='Which samples\' variants to extract (default: all samples).')
  parser.add_argument('--nofilt', action='store_true', help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

  args = parser.parse_args()
  if not os.path.isdir(args.outdir):
    print("[ERROR] Directory does not exist: '%s'." % args.outdir, file=sys.stderr)
    sys.exit(1)

  return args

def main(fh_input, dir_out, id_samples=None, nofilt=False):
  '''Extract filter-passing variants from VCF file.'''
  muts = {}
  rdr = vcf.Reader(fh_input)
  for rec in rdr:
    if nofilt or rec.FILTER is None or len(rec.FILTER) == 0:
      calls = []
      # filter by sample if parameter given
      for x in rec.samples:
        if (id_samples==None) or (x.sample in id_samples):
          if x.sample not in muts:
            muts[x.sample] = []
          calls.append(x)
      # parse filtered variant calls
      for call in calls:
        rc_ref, rc_alt = call.data.AD[0:2]
        muts[call.sample].append({
          'id': '%s:%d' % (rec.CHROM, rec.POS),
          'ref_counts': rc_ref,
          'var_counts': rc_alt,
          'states': [{'g_n':'AA', 'g_r':'AA', 'g_v':x, 'prior_weight':1} for x in ['AA', 'AB']]
      })
  for sample in muts:
    fn_out = os.path.join(dir_out, sample + ".yaml")
    with open(fn_out, 'wt') as fh_out:
      data_out = {'mutations': muts[sample]}
      fh_out.write(yaml.dump(data_out))

if __name__ == '__main__':
  args = parse_args()
  main(args.infile, args.outdir, args.samples, args.nofilt)
