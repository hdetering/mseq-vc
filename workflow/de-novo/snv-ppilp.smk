# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

rule snv_ppilp:
  input:
    "haplotypecaller.filt.vcf"
  output:
    "snv-ppilp.out.txt"
  log:
    "log/snv-ppilp.log"
  params:
  threads: 1
  shell:
    """
    time (
    python2 {config[tools][snv-ppilp]} \
      -i {input} \
      -o {output}
    ) 1>{log} 2>&1
    """

rule snv_ppilp_to_vcf:
  input:
    "snv-ppilp.out.txt"
  output:
    "snv-ppilp.vcf"
  log:
    "log/snv-ppilp.to_vcf.log"
  threads: 1
  run:
    import re

    samples = []
    chr_pos_sample = {}

    f = open(input[0], 'rt')
    # discard header ("Sample_ID", "list of SNVs")
    f.readline()
    for line in f:
      # fields are expected to be separated by commas
      sample_snvs = line.strip().split(', ', 1)
      # skip lines with only one field
      if len(sample_snvs) == 1:
        samples.append(sample_snvs[0])
        continue
      sample, snvs = sample_snvs
      samples.append(sample)
      for m in re.finditer(r'\((?P<chr>\w+):(?P<pos>\d+)\)', snvs):
        chr = m.group('chr')
        pos = int(m.group('pos'))
        if chr not in chr_pos_sample:
          chr_pos_sample[chr] = {}
        if pos not in chr_pos_sample[chr]:
          chr_pos_sample[chr][pos] = []
        chr_pos_sample[chr][pos].append(sample)

    # output variants in multi-sample VCF format
    with open(output[0], 'wt') as vcf:
      # write headers
      vcf.write('##fileformat=VCFv4.1\n')
      vcf.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">\n')
      vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
      vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % '\t'.join(samples))
      # write variants
      for chr in chr_pos_sample:
        for pos in sorted(chr_pos_sample[chr].keys()):
          info = "NS=%d" % len(chr_pos_sample[chr][pos])
          parts = [chr, str(pos), '.', '.', '.', '.', 'PASS', info, 'GT']
          for s_id in samples:
            gt = '0/1' if (s_id in chr_pos_sample[chr][pos]) else '.'
            parts.append(gt)
          vcf.write('%s\n' % '\t'.join(parts))
