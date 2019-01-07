#!/usr/bin/env python3
import os, sys
import argparse
import vcf

def parse_args():
  parser = argparse.ArgumentParser(description='Extract variant allele freqs from pileup VCF.')
  parser.add_argument('vcf_pileup', help='Pileup VCF (bcftools mpileup)')  
  parser.add_argument('vcf_calls', help='VCF with variant calls')
  parser.add_argument('vcf_true', help='VCF with true variants')
  parser.add_argument('prefix', help='Identifier field to prepend output rows')

  args = parser.parse_args()
  if not os.path.exists(args.vcf_calls):
    print('[ERROR] file does not exist: {}'.format(args.vcf_calls, file=sys.stderr))
    sys.exit(1)
  if not os.path.exists(args.vcf_pileup):
    print('[ERROR] file does not exist: {}'.format(args.vcf_pileup, file=sys.stderr))
    sys.exit(1)
   
  return args

def extract_vafs(vcf_pileup, vars_called, vars_true, pfx):
  vcf_reader = vcf.Reader(open(vcf_pileup), 'r')
  for rec in vcf_reader:
    alt = ''
    # check if pileup position occurs in variant calls
    if rec.POS in vars_called:
      alt = vars_called[rec.POS].ALT[0]
    elif rec.POS in vars_true:
      alt = vars_true[rec.POS].ALT[0]
    else:
      print('[ERROR] POS does not match calls or truth: {}'.format(rec.POS), file=sys.stderr)
      print('[ERROR] bailing out...')
      sys.exit(1)

    # determine index of called ALT allele
    allele_found = False
    idx = 0
    while idx < len(rec.ALT) and not allele_found:
      if rec.ALT[idx] == alt:
        allele_found = True
      else:
        idx += 1
    if not allele_found:
      print("[ERROR] allele '{}' not found in pileup (POS: {})".format(alt, rec.POS), file=sys.stderr)

    # calulate ALT allele VAF for each sample 
    vaf_max = 0.0
    if allele_found:     
      vafs = []
      for call in rec.samples:
        dp = call.data.DP
        ad = call.data.AD[idx+1]
        vafs.append(ad/dp)
      vaf_max = max(vafs)
    
    print(','.join([pfx, str(rec.CHROM), str(rec.POS), '{:.4f}'.format(vaf_max)]))


def scan_var_calls(vcf_calls):
  variants = {}
  vcf_reader = vcf.Reader(open(vcf_calls, 'rt'))

  for rec in vcf_reader:
    # multi-allelic variants are not supported at present
    if len(rec.ALT) != 1:
      print('[ERROR] multi-allelic call at POS {}'.format(rec.POS), file=sys.stderr)
      print('[ERROR] bailing out...')
      sys.exit(1)

    variants[rec.POS] = rec

  return variants

def main(args):

  print('[INFO] importing variant calls from VCF:\n{}'.format(args.vcf_calls), file=sys.stderr)
  vars_called = scan_var_calls(args.vcf_calls)

  print('[INFO] importing true variants from VCF:\n{}'.format(args.vcf_true), file=sys.stderr)
  vars_true = scan_var_calls(args.vcf_true)

  print('[INFO] scanning, sites in pileup VCF:\n{}'.format(args.vcf_pileup), file=sys.stderr)
  extract_vafs(args.vcf_pileup, vars_called, vars_true, args.prefix)

if __name__ == '__main__':
  args = parse_args()
  main(args)
