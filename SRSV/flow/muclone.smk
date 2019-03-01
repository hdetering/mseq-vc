# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

from glob import glob
#BAMS_TUM = glob_wildcards("data/%s/{sample,R\d+\.bam}" % REPS)
SAMPLES = ['R1', 'R2', 'R3', 'R4', 'R5']

rule muclone_tsv:
  input:
    rc=["sim/%s.vars.csv" % x for x in SAMPLES]
  output:
    tsv=["muclone/%s.tsv" % x for x in SAMPLES]
  threads: 1
  run:
    for i in range(len(input.rc)):
      with open(output.tsv[i], 'wt') as f_tsv:
        # write header
        hdr = ["mutation_id","ref_counts","var_counts","normal_cn","minor_cn","major_cn"]
        f_tsv.write('\t'.join(hdr) + '\n')
        for line in open(input.rc[i], 'rt'):
          if line.startswith('s'):
            mid, tot, alt = line.strip().split('\t')
            rc_tot = int(tot)
            rc_alt = int(alt)
            rc_ref = rc_tot - rc_alt
            vaf = float(rc_alt)/rc_tot
            f_tsv.write("%s\t%d\t%d\t2\t1\t1\t%f\n" % (mid,rc_ref,rc_alt,vaf))

#rule muclone_yaml:
#  input:
#    rc=["sim/%s.vars.csv" % x for x in SAMPLES]
#  output:
#    yml=["muclone/%s.yaml" % x for x in SAMPLES]
#  threads: 1
#  run:
#    import yaml
#    for i in range(len(input.rc)):
#      muts = []
#      for line in open(input.rc[i], 'rt'):
#        if line.startswith('s'):
#          mid, tot, alt = line.strip().split('\t')
#          rc_tot = int(tot)
#          rc_alt = int(alt)
#          rc_ref = rc_tot - rc_alt
#          m = {
#            'id': mid,
#            'ref_counts': rc_ref,
#            'var_counts': rc_alt,
#            'states': [{'g_n':'AA', 'g_r':'AA', 'g_v':x, 'prior_weight': 1} for x in ['AB', 'BB']]
#          }
#          muts.append(m)
#      out = { 'mutations' : muts }
#      with open(output.yml[i], 'wt') as f_yml:
#        f_yml.write(yaml.dump(out))
rule muclone_yaml:
  input:
    vcf="mutect2/{sample}.raw.vcf"
  output:
    yml="muclone/{sample}.yaml"
  threads: 1
  shell:
    """
    python {config[scripts]}/mutect2yaml.py --nofilt {input.vcf} {output.yml}
    """
#rule muclone_yaml:
#  input:
#    vcf="bcftools.vcf"
#  output:
#    yml=["muclone/%s.yaml" % x for x in SAMPLES]
#  params:
#    outdir="muclone/"
#  threads: 1
#  shell:
#    """
#    python {config[scripts]}/bcftools2yaml.py --nofilt {input.vcf} {params.outdir}
#    """

#rule muclone_mm2clust:
#  input:  mm="sim/mm.csv"
#  output: tsv="muclone/clusters.tsv"
#  threads: 1
#  run:
#    import pandas as pd
#
#    clusters = {}
#    # read mutation matrix (named rows, no header, cells are 0 or 1)
#    df = pd.read_csv(input.mm, header=None, index_col=0)
#    n = len(df.columns)
#    # loop over columns (each column represents a mutation)
#    for colname in df:
#      col = df[colname]
#      clones = col[col==1].index
#      id_clust = '_'.join(clones)
#      if id_clust in clusters:
#        clusters[id_clust] += 1
#      else:
#        clusters[id_clust] = 1
#    # output
#    with open(output.tsv, 'wt') as f_out: 
#      for k in sorted(clusters.keys()):
#        f_out.write("%s\t%f\t%d\n" % (k, clusters[k]/n, clusters[k]))

#rule muclone_cluster2phi:
#  input:
#    clus="muclone/clusters.tsv",
#    prev="sim/clone_prev.csv"
#  output:
#    ["muclone/tmp/%s_cluster2Phi.tsv" % x for x in SAMPLES]
#  threads: 1
#  run:
#    import pandas as pd
#    # read clusters
#    clusters = pd.read_csv(input.clus, sep='\t', header=None)
#    # read prevalence matrix (rows: samples, cols: clones)
#    prev = pd.read_csv(input.prev, index_col=0)
#    for id_sample, prow in prev.iterrows():
#      out = [['prior', 'phi']]
#      for idx, crow in clusters.iterrows():
#        id_clust = crow[0].strip().split('_')
#        prior = crow[1]
#        phi = prow[id_clust].sum()
#        out.append([prior, phi])
#      with open("muclone/tmp/%s_cluster2Phi.tsv" % id_sample, 'wt') as f_out:
#        for elements in out:
#          f_out.write('\t'.join([str(e) for e in elements]) + '\n')

rule muclone_pyclone2clust:
  input:
    clusters="pyclone/result.clusters.tsv"
  output:
    ["muclone/tmp/%s_cluster2Phi.tsv" % x for x in SAMPLES]
  threads: 1
  run:
    import pandas as pd
    # read clusters
    # sample_id  cluster_id  size  mean  std
    clusters = pd.read_csv(input.clusters, sep='\t')
    # total number of mutations across all clusters
    n_mut = clusters[['cluster_id', 'size']].drop_duplicates()['size'].sum()
    # order clusters by samples
    samples = {}
    for idx, row in clusters.iterrows():
      if row['sample_id'] not in samples:
        samples[row['sample_id']] = []
      samples[row['sample_id']].append((float(row['size'])/n_mut, row['mean']))
    # output cluster frequencies to files
    for sid in samples:
      with open('muclone/tmp/%s_cluster2Phi.tsv' % sid, 'wt') as f_out:
        f_out.write('prior\tphi\n')
        f_out.write('\n'.join(['%g\t%g' % (p, q) for p, q in samples[sid]]) + '\n')

rule muclone_conf:
  input:
    yaml=["muclone/%s.yaml" % x for x in SAMPLES],
    tsv=["muclone/tmp/%s_cluster2Phi.tsv" % x for x in SAMPLES]
  output:
    cfg="muclone/config.yaml"
  run:
    import os, yaml
    cfg = {
      'working_dir': os.path.join(os.getcwd(), 'muclone'),
      'trace_dir'  : 'tmp',
      'samples'    : { x: {
        'mutations_file': "%s.yaml" % x,
        'flat_cluster_files': "tmp/%s_cluster2Phi.tsv" % x,
        'tumour_content': { 'value': 1.0 },
        'error_rate': 0.01
        } for x in SAMPLES}
    }
    with open(output.cfg, 'wt') as f_out:
      f_out.write( yaml.dump(cfg) )

rule muclone:
  input:
    yml=["muclone/%s.yaml" % x for x in SAMPLES],
    cfg="muclone/config.yaml",
    vcf="mutect2.vcf"
  output: "muclone.vcf"
  log:    "log/muclone.log"
  params:
    #dir="data/{replicate}"
  threads: 1
  shell:
    """
    time (
    
    set +u
    source activate {config[tools][muclone_dir]}/env
    export PYTHONPATH="{config[tools][muclone_dir]}/src:$PYTHONPATH"
    set -u

    python {config[tools][muclone_dir]}/src/mutation_classifier/main.py analysis \
      --config_file {input.cfg}
    
    if [[ -d "muclone/tmp/trace" ]]; then
      for f in muclone/tmp/trace/*/MuClone-sample-results.tsv; do
        awk -f {config[scripts]}/mutect2muclone.awk {input.vcf} $f > {output}
      done
    fi

    ) >{log} 2>&1
    """
