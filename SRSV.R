#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Load benchmarking data, perform analyses and generate plots for publication.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-04-25
#------------------------------------------------------------------------------

require(tidyverse)
require(dbplyr)
require(RSQLite)
require(cowplot)
require(ggpubr)
require(RColorBrewer)
require(gridExtra)
require(grid)

# define input/output paths
data_dir <- file.path( 'data', 'SRSV' )
plot_dir <- file.path( 'plot', 'SRSV' )

# source required scripts
source( file.path('analysis', 'performance.analysis.R') )
source( file.path('plotting', 'performance.plotting.R') )
source( file.path('analysis', 'upset.analysis.R') )
source( file.path('plotting', 'upset.plotting.R') )

# connection to analysis database
#-------------------------------------------------------------------------------
db <- file.path( data_dir, 'analysis.db' )
con <- DBI::dbConnect(RSQLite::SQLite(), db)
df_rep <- tbl( con, 'replicates' ) %>% collect()
df_mut <- tbl( con, 'mutations' ) %>% collect()
df_mut_sample <- tbl( con, 'mutations_samples' ) %>% collect()
df_mut_clone <- tbl( con, 'clones_mutations' ) %>% collect()
df_prev <- tbl( con, 'clones_prev' ) %>% collect()
df_caller <- tbl( con, 'callers' ) %>% collect()
df_varcall <- tbl( con, 'varcalls' ) %>% collect()
df_rc <- tbl( con, 'readcounts' ) %>% collect()
df_snp <- tbl( con, 'snps' ) %>% collect()

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
df_vars <- classify_variants( df_varcall, df_mut, df_mut_sample )
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
p_perf <- plot_perf( df_perf )
ggsave( file.path( plot_dir, 'SRSV.performance.pdf'), plot = p_perf, width = 10, height = 6)
ggsave( file.path( plot_dir, 'SRSV.performance.png'), plot = p_perf, width = 10, height = 6)

################################################################################
# UpSet plots
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

df_pres <- get_upset_pres( df_vars, df_rc )
df_pres %>% write_csv( file.path(data_dir, 'SRSV.muts_callsets.csv') )

df_pres <- read.csv( file.path(data_dir, 'SRSV.muts_callsets.csv') )
# annoying, but necessary... only if loaded from file
names(df_pres)[names(df_pres)=='SNV.PPILP'] <- 'SNV-PPILP'

# plot variant calls in relation to TRUE somatic variants
df <- df_pres
n <- c(callers$name_caller, 'TRUE_somatic')
pdf( file.path( plot_dir, 'SRSV.upset.som.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'FP' )
n <- c(callers$name_caller, 'TRUE_germline')
pdf( file.path( plot_dir, 'SRSV.upset.FP.GL.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'TP' )
n <- c(callers$name_caller, 'TRUE_somatic')
pdf( file.path( plot_dir, 'SRSV.upset.TP.som.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
