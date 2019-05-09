#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Load benchmarking data, perform analyses and generate plots for publication.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-05-07
#------------------------------------------------------------------------------#

require(tidyverse)
require(dbplyr)
require(RSQLite)
require(cowplot)
require(ggpubr)
require(RColorBrewer)
require(gridExtra)
require(grid)

# define input/output paths
data_dir.srsv <- file.path( 'data', 'SRSV' )
plot_dir.srsv <- file.path( 'plot', 'SRSV' )
data_dir.rrsv <- file.path( 'data', 'RRSV' )
plot_dir.rrsv <- file.path( 'plot', 'RRSV' )

# source required scripts
source( file.path('analysis', 'performance.analysis.R') )
source( file.path('plotting', 'performance.plotting.R') )
source( file.path('analysis', 'upset.analysis.R') )
source( file.path('plotting', 'upset.plotting.R') )
source( file.path('plotting', 'rc_alt.plotting.R') )

# load SRSV data
#-------------------------------------------------------------------------------
db.srsv <- file.path( data_dir.srsv, 'analysis.db' )
con.srsv <- DBI::dbConnect(RSQLite::SQLite(), db.srsv)
df_rep.srsv <- tbl( con.srsv, 'replicates' ) %>% collect()
df_mut.srsv <- tbl( con.srsv, 'mutations' ) %>% collect()
df_mut_sample.srsv <- tbl( con.srsv, 'mutations_samples' ) %>% collect()
df_mut_clone.srsv <- tbl( con.srsv, 'clones_mutations' ) %>% collect()
df_prev.srsv <- tbl( con.srsv, 'clones_prev' ) %>% collect()
df_caller.srsv <- tbl( con.srsv, 'callers' ) %>% collect()
df_varcall.srsv <- tbl( con.srsv, 'varcalls' ) %>% collect()
df_rc.srsv <- tbl( con.srsv, 'readcounts' ) %>% collect()
df_snp.srsv <- tbl( con.srsv, 'snps' ) %>% collect()

# load RRSV data
#-------------------------------------------------------------------------------
df_rep.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.replicates.rds') )
df_mut.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.mutations.rds') )
df_mut_sample.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.mutations_samples.rds') )
df_mut_clone.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.clones_mutations.rds') )
df_prev.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.clones_prev.rds') )
df_caller.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.callers.rds') )
df_varcall.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.varcalls.rds') ) %>% 
  mutate(chrom = as.character(chrom), id_caller = as.integer(id_caller))
df_rc.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.readcounts.rds') )
# df_snp.rrsv <- readRDS( file.path(data_dir.rrsv, 'RRSV.snps.rds') ) # does not exist

# define variant calling tools and their properties
#-------------------------------------------------------------------------------
callers <- tibble(
  name_caller = c(
    'Bcftools', 
    'HaplotypeCaller', 
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
    'MultiSNV', 
    'SNV-PPILP',
    'Mutect2_mseq'
  ),
  class = c(rep('general', 2), rep('tumor-normal', 11), rep('multi-sample', 4))
)
df_caller <- df_caller %>%
  inner_join( callers, by = 'name_caller' )

# determine status of variant calls
#   TP: true positives
#   FP: false positives
#   FN: false negatives
#-------------------------------------------------------------------------------
df_vars.srsv <- classify_variants( df_varcall.srsv, df_mut.srsv, df_mut_sample.srsv, df_caller.srsv )
df_vars.rrsv <- classify_variants( df_varcall.rrsv, df_mut.rrsv, df_mut_sample.rrsv, df_caller.rrsv )
# store variant calls for later use
saveRDS( df_vars.srsv, file.path(data_dir.srsv, 'df_vars.rds') )
saveRDS( df_vars.rrsv, file.path(data_dir.rrsv, 'df_vars.rds') )


# calculate performance metrics
# ------------------------------------------------------------------------------
#df_vars.srsv <- readRDS( file.path(data_dir.srsv, 'df_vars.rds') )
df_perf.srsv <- calculate_performance( df_vars.srsv, df_caller.srsv, df_rep.srsv )
df_perf.rrsv <- calculate_performance( df_vars.rrsv, df_caller.rrsv, df_rep.rrsv )
# write summary stats to file
#saveRDS( df_perf, file.path(data_dir, 'df_perf.rds') )

# plot performance metrics for SRSV 30x and RRSV 25x
# ------------------------------------------------------------------------------
df.srsv <- df_perf.srsv %>% dplyr::filter( cvg == 50 ) %>%
  select(id_rep, FN, FP, TP, recall, precision, F1, caller = name_caller) %>%
  mutate( cvg = 'SRSV')
df.rrsv <- df_perf.rrsv %>% 
  select(id_rep, FN, FP, TP, recall, precision, F1, caller = name_caller) %>%
  mutate( cvg = 'RRSV' )
df.rrsv$caller[df.rrsv$caller == 'bcftools'] <- 'Bcftools'
df.rrsv$caller[df.rrsv$caller == 'Haplotype Caller'] <- 'HaplotypeCaller'
df.rrsv$caller[df.rrsv$caller == 'Somatic Sniper'] <- 'SomaticSniper'

df <- df.srsv %>% bind_rows( df.rrsv ) %>% inner_join( callers, by = c('caller'='name_caller') )
p_perf <- plot_perf( df )
ggsave( file.path( plot_dir, 'SRSV.performance.pdf'), plot = p_perf, width = 10, height = 6)
ggsave( file.path( plot_dir, 'SRSV.performance.png'), plot = p_perf, width = 10, height = 6)

# plot ALT read counts for TP, FP, FN variants
# ------------------------------------------------------------------------------

df_vars.rrsv <- readRDS( file.path(data_dir.rrsv, 'df_vars.rds') )
df_vars.rrsv <- df_caller.rrsv %>% inner_join( df_vars.rrsv, by = 'id_caller' )

p_rc_alt.rrsv <- plot_rc_alt_ridges_rrsv( df_vars.rrsv, df_rc.rrsv, df_rep.rrsv )
ggsave( file.path( plot_dir.rrsv, 'RRSV.rc_alt.ridges.pdf'), plot = p_rc_alt.rrsv, device = pdf(), width = 8, height = 10 )
ggsave( file.path( plot_dir.rrsv, 'RRSV.rc_alt.ridges.png'), plot = p_rc_alt.rrsv, device = png(), width = 8, height = 10 )
