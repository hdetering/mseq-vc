require(tidyverse)
require(dbplyr)
require(RSQLite)
require(cowplot)
require(ggpubr)
require(RColorBrewer)
require(gridExtra)
require(grid)

# source required scripts
source( file.path('analysis', 'performance.analysis.R') )
source( file.path('analysis', 'admixture.analysis.R') )
source( file.path('analysis', 'similarity.analysis.R') )
source( file.path('analysis', 'upset.analysis.R') )
source( file.path('plotting', 'performance.plotting.R') )
source( file.path('plotting', 'rc_alt.plotting.R') )
source( file.path('plotting', 'similarity.plotting.R') )
source( file.path('plotting', 'upset.plotting.R') )

# define input/output paths
data_dir <- file.path( 'data', 'RRSV' )
plot_dir <- file.path( 'plot', 'RRSV' )

# load RRSV data
#-------------------------------------------------------------------------------
df_rep <- readRDS( file.path(data_dir, 'RRSV.replicates.rds') )
df_mut <- readRDS( file.path(data_dir, 'RRSV.mutations.rds') )
df_mut_sample <- readRDS( file.path(data_dir, 'RRSV.mutations_samples.rds') )
df_mut_clone <- readRDS( file.path(data_dir, 'RRSV.clones_mutations.rds') )
df_prev <- readRDS( file.path(data_dir, 'RRSV.clones_prev.rds') )
df_caller <- readRDS( file.path(data_dir, 'RRSV.callers.rds') )
df_varcall <- readRDS( file.path(data_dir, 'RRSV.varcalls.rds') ) %>% 
  mutate(chrom = as.character(chrom), id_caller = as.integer(id_caller))
df_rc <- readRDS( file.path(data_dir, 'RRSV.readcounts.rds') )
# df_snp.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.snps.rds') ) # does not exist

# define variant calling tools and their properties
#-------------------------------------------------------------------------------
callers <- tibble(
  name_caller = c(
    'Bcftools', 
    'CaVEMan', 
    'MuTect1', 
    'Mutect2', 
    'NeuSomatic', 
    'Shimmer', 
    'SNooPer', 
    'SomaticSniper', 
    'Strelka1', 
    'Strelka2', 
    'VarDict', 
    'VarScan',
    'MuClone', 
    'SNV-PPILP',
    'HaplotypeCaller', 
    'MultiSNV', 
    'Mutect2_ms'
  ),
  class = c(rep('marginal', 12), rep('two-step', 2), rep('joint', 3))
)

# determine status of variant calls
#   TP: true positives
#   FP: false positives
#   FN: false negatives
#-------------------------------------------------------------------------------
df_vars <- classify_variants( df_varcall, df_mut, df_mut_sample, df_caller )
df_vars <- df_vars %>% 
  left_join( df_snp, by = c('id_rep', 'chrom', 'pos') ) %>% 
  mutate( germline = (!is.na(id_mut)) ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos, type, germline )
# store variant calls for later use
saveRDS( df_vars, file.path(data_dir, 'df_vars.rds') )

# calculate performance metrics
# ------------------------------------------------------------------------------
df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_perf <- calculate_performance( df_vars, df_caller, df_rep )
# write summary stats to file
saveRDS( df_perf, file.path(data_dir, 'df_perf.rds') )

# plot performance metrics
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
p_perf <- plot_perf_rrsv( df_perf )
ggsave( file.path( plot_dir, 'Fig4.RRSV.performance.cvg.pdf'), plot = p_perf, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig4.RRSV.performance.cvg.png'), plot = p_perf, width = 8, height = 10)


################################################################################
# Variant allele freqs for TP, FP, FN variants
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

# p_vaf <- plot_vaf_dens_srsv( df_vars, df_rc, df_rep )
# ggsave( file.path( plot_dir, 'Fig3.SRSV.vaf.dens.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 10 )
# ggsave( file.path( plot_dir, 'Fig3.SRSV.vaf.dens.png'), plot = p_vaf, device = png(), width = 8, height = 10 )
# p_vaf <- plot_rc_alt_ridges_srsv( df_vars, df_rc, df_rep )
# ggsave( file.path( plot_dir, 'Fig3.SRSV.vaf.ridges.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 10 )
# ggsave( file.path( plot_dir, 'Fig3.SRSV.vaf.ridges.png'), plot = p_vaf, device = png(), width = 8, height = 10 )
p_vaf <- plot_vaf_bar_srsv( df_vars, df_rc, df_rep )
ggsave( file.path( plot_dir, 'Fig5.spike-in.vaf.bar.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'Fig5.spike-in.vaf.bar.png'), plot = p_vaf, device = png(), width = 8, height = 8 )


################################################################################
# Similarity between call sets
################################################################################

df_var <- readRDS( file.path(data_dir, 'df_vars.rds') )

## true positives
df_pres_tp <- get_var_pres( df_var, df_caller, 'TP' )
df_jacc_tp <- Jaccard.df( df_pres_tp %>% select(-id_mut) )
p_jacc_tp <- plot_jacc_idx( df_jacc_tp )

## false positives
df_pres_fp <- get_var_pres( df_var, df_caller, 'FP' )
df_jacc_fp <- Jaccard.df( df_pres_fp %>% select(-id_mut) )
p_jacc_fp <- plot_jacc_idx( df_jacc_fp )

## false negatives
df_pres_fn <- get_var_pres( df_var, df_caller, 'FN' )
df_jacc_fn <- Jaccard.df( df_pres_fn %>% select(-id_mut) )
p_jacc_fn <- plot_jacc_idx( df_jacc_fn )

## multi-plot
p_jacc_multi <- plot_jacc_idx_multi( p_jacc_tp, p_jacc_fn, p_jacc_fp )
ggsave( file.path( plot_dir, 'Fig8.spike-in.jaccard.pdf'), plot = p_jacc_multi, device = pdf(), width = 10, height = 4 )
ggsave( file.path( plot_dir, 'Fig8.spike-in.jaccard.png'), plot = p_jacc_multi, device = png(), width = 10, height = 4 )
