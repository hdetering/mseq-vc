#!/usr/bin/env python3
# vim: syntax=python tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert Mutect2 VCF output to PyClone/MuClone YAML input.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-10-08
#------------------------------------------------------------------------------

from __future__ import print_function
import os, sys
import argparse
import vcf
import yaml

def parse_args():
  parser = argparse.ArgumentParser(description='Convert Mutect2 VCF output to PyClone YAML input.')
  parser.add_argument('infile', type=argparse.FileType(), help='Mutect2 VCF file.')
  parser.add_argument('outfile', type=argparse.FileType('wt'), help='Output YAML file.')
  parser.add_argument('--sample', help='Which sample\'s variants to extract (default: detect from VCF header).')
  parser.add_argument('--nofilt', action='store_true', help='Disable filtering; otherwise only "PASS" variants are output (default: off).')

  args = parser.parse_args()
  return args

def main(fh_input, fh_output, sample=None, nofilt=False):
  '''Extract filter-passing biallelic SNVs from VCF file.'''
  muts = []
  rdr = vcf.Reader(fh_input)
  
  # attempt to determine tumor sample from VCF header
  lbl_sample = ''
  idx_sample = 0
  if sample:
    if  sample in rdr.samples:
      lbl_sample = sample
    else:
      print('[WARN] sample {} not present in VCF file. Attempting to guess tumor sample from header.'.format(sample))
  if not sample and 'tumor_sample' in rdr.metadata:
    lbl_sample = rdr.metadata['tumor_sample'][0]
  if lbl_sample in rdr.samples:
    idx_sample = rdr.samples.index(lbl_sample)
  print('[INFO] exporting records for sample "{}".'.format(rdr.samples[idx_sample]))  

  for rec in rdr:
    if nofilt or rec.FILTER is None or len(rec.FILTER) == 0:
      call = None
      call = rec.samples[idx_sample]
      
      # make sure variant is biallelic SNV
      if len(rec.REF) > 1 or len(rec.ALT) > 1 or len(rec.ALT[0]) > 1:
        continue
      
      rc_ref, rc_alt = call.data.AD[0:2]
      muts.append({
        'id': '%s:%d' % (rec.CHROM, rec.POS),
        'ref_counts': rc_ref,
        'var_counts': rc_alt,
        #'states': [{'g_n':'AA', 'g_r':'AA', 'g_v':x, 'prior_weight':1} for x in ['AB', 'BB']]
        'states': [{'g_n':'AA', 'g_r':'AA', 'g_v':x, 'prior_weight':1} for x in ['AB']]
      })
  data_out = {'mutations': muts}
  fh_output.write(yaml.dump(data_out))

if __name__ == '__main__':
  args = parse_args()
  main(args.infile, args.outfile, args.sample, args.nofilt)
