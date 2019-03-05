# vim: syntax=python tabstop=2 expandtab
# coding: utf-8

from glob import glob
#BAMS_TUM = glob_wildcards("data/%s/{sample,R\d+\.bam}" % REPS)
SAMPLES = "R1 R2 R3 R4 R5".split()

rule pyclone_yaml:
  input:
    vcf="mutect2/{sample}.vcf"
  output:
    yml="pyclone/{sample}.yaml"
  threads: 1
  shell:
    '''
    python3 {config[scripts]}/mutect2yaml.py {input} {output}
    '''

rule pyclone_conf:
  input:
    yaml=["pyclone/%s.yaml" % x for x in SAMPLES]
  output:
    cfg="pyclone/config.yaml"
  run:
    import os, yaml
    cfg = {
      # Specifies working directory for analysis. All paths in the rest of the file are relative to this.
      'working_dir': os.path.join(os.getcwd(), 'pyclone'),
      # Where the trace (output) from the PyClone MCMC analysis will be written.
      'trace_dir'  : 'trace',
      # Specifies which density will be used to model read counts. Most people will want pyclone_beta_binomial or pyclone_binomial
      'density'    : 'pyclone_beta_binomial',
      # Number of iterations of the MCMC chain.
      'num_iters'  : 10000,
      
      # Specifies parameters in Beta base measure for DP. Most people will want the values below.
      'base_measure_params': {
        'alpha' : 1,
        'beta'  : 1
      },
      
      # Specifies initial values and prior parameters for the prior on the concentration (alpha) parameter in the DP. If the prior node is not set the concentration will not be estimated and the specified value will be used.
      'concentration': {
        # Initial value if prior is set, or fixed value otherwise for concentration parameter.
        'value' : 1.0,
        # Specifies the parameters in the Gamma prior over the concentration parameter.
        'prior' : {
          'shape' : 1.0,
          'rate'  : 0.001
        }
      },
      
      # Beta-Binomial precision (alpha + beta) prior
      'beta_binomial_precision_params': {
        # Starting value
        'value': 1000,
        # Parameters for Gamma prior distribution
        'prior': {
          'shape': 1.0,
          'rate': 0.0001
        },
        # Precision of Gamma proposal function for MH step
        'proposal': {
          'precision': 0.01
        }
      },
      
      # Contains one or more sub-entries which specify details about the samples used in the analysis.
      'samples': { x: {
        'mutations_file': "%s.yaml" % x,
        'tumour_content': { 'value': 1.0 },
        'error_rate': 0.01
        } for x in SAMPLES}
    }
    with open(output.cfg, 'wt') as f_out:
      f_out.write( yaml.dump(cfg) )

rule pyclone:
  input:
    cfg="pyclone/config.yaml"
  output: 
    loci="pyclone/result.loci.tsv",
    clust="pyclone/result.clusters.tsv"
  log:  
    "log/pyclone.log"
  params:
    #dir="data/{replicate}"
  threads: 1
  shell:
    """
    time (
   
    source activate pyclone

    PyClone run_analysis --config_file {input.cfg} 
    PyClone build_table --config_file {input.cfg} --out_file {output.loci} --table_type loci
    PyClone build_table --config_file {input.cfg} --out_file {output.clust} --table_type cluster

    ) >{log} 2>&1
    """
