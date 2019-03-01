#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert Strelka2 VCF output to PyClone/MuClone YAML input.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-02-15
#------------------------------------------------------------------------------

import os, sys
import argparse
import vcf
import yaml

def parse_args():
  parser = argparse.ArgumentParser(description='Convert Strelka2 VCF output to PyClone YAML input.')
  parser.add_argument('infile', type=argparse.FileType(), help='Strelka2 VCF file.')
  parser.add_argument('outfile', type=argparse.FileType('wt'), help='Output YAML file.')
  parser.add_argument('--sample', help='Which sample\'s variants to extract (default: last sample).')
  parser.add_argument('--nofilt', action='store_true', help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

  args = parser.parse_args()
  return args

def main(fh_input, fh_output, id_sample=None, nofilt=False):
  '''Extract filter-passing variants from VCF file.'''
  muts = []
  rdr = vcf.Reader(fh_input)
  for rec in rdr:
    if nofilt or rec.FILTER is None or len(rec.FILTER) == 0:
      call = None
      # filter by sample if parameter given
      if id_sample:
        for x in rec.samples:
          if x.sample == id_sample:
            call = x
            break
      else:
        call = rec.samples[-1]
      # extract read counts for all bases
      rc = {
        'A': call.data.AU[0],
        'C': call.data.CU[0],
        'G': call.data.GU[0],
        'T': call.data.TU[0]
      }
      rc_ref = rc[rec.REF]
      rc_alt = rc[str(rec.ALT[0])]
      muts.append({
        'id': '%s:%d' % (rec.CHROM, rec.POS),
        'ref_counts': rc_ref,
        'var_counts': rc_alt,
        'states': [{'g_n':'AA', 'g_r':'AA', 'g_v':x, 'prior_weight':1} for x in ['AA', 'AB']]
      })
  data_out = {'mutations': muts}
  fh_output.write(yaml.dump(data_out))

if __name__ == '__main__':
  args = parse_args()
  main(args.infile, args.outfile, args.sample, args.nofilt)
