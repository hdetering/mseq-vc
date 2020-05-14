#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Load benchmarking data, perform analyses and generate plots for publication.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-01-29
#------------------------------------------------------------------------------

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
# source( file.path('analysis', 'tumor-wise.analysis.R') )
source("/Users/laura/GoogleDrive/PROJECTS/PROYECTOS/TestingVC_withTama/05_SCRIPTS/new_version_paper_20200514/tumor-wise.analysis.R")
source( file.path('analysis', 'ITH.analysis.R') )
source( file.path('analysis', 'admixture.analysis.R') )
source( file.path('analysis', 'similarity.analysis.R') )
source( file.path('analysis', 'upset.analysis.R') )
source( file.path('plotting', 'performance.plotting.R') )
source( file.path('plotting', 'tumor-wise.plotting.R') )
source( file.path('plotting', 'ITH.plotting.R') )
source( file.path('plotting', 'rc_alt.plotting.R') )
source( file.path('plotting', 'similarity.plotting.R') )
source( file.path('plotting', 'upset.plotting.R') )

# define input/output paths
data_dir <- file.path( 'data', 'spike-in' )
plot_dir <- file.path( 'plot', 'spike-in' )

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
df_snp <- readRDS( file.path(data_dir, 'RRSV.snps.rds') ) 

# define variant calling tools and their properties
#-------------------------------------------------------------------------------
# df_caller <- tibble(
#   name_caller = c(
#     'Bcftools', 
#     'CaVEMan', 
#     'MuTect1', 
#     'Mutect2_single', 
#     'NeuSomatic', 
#     'Shimmer', 
#     'SNooPer', 
#     'SomaticSniper', 
#     'Strelka1', 
#     'Strelka2', 
#     'VarDict', 
#     'VarScan',
#     'MuClone', 
#     'SNV-PPILP',
#     'HaplotypeCaller', 
#     'MultiSNV', 
#     'Mutect2_multi_F'
#   ),
#   class = c(rep('marginal', 12), rep('two-step', 2), rep('joint', 3))
# )

# determine status of variant calls
#   TP: true positives
#   FP: false positives
#   FN: false negatives
#-------------------------------------------------------------------------------
df_vars <- classify_variants( df_varcall, df_mut, df_mut_sample, df_caller )
df_vars <- df_vars %>% 
  left_join( df_snp, by = c( 'chrom', 'pos') ) %>% 
  mutate( germline = (!is.na(id_mut)) ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos, type, germline )
# store variant calls for later use
saveRDS( df_vars, file.path(data_dir, 'df_vars.rds') )

# calculate performance metrics
# ------------------------------------------------------------------------------
df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_perf <- calculate_performance_sample( df_vars, df_caller, df_rep )
# write summary stats to file
saveRDS( df_perf, file.path(data_dir, 'df_perf.rds') )

# plot performance metrics
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
# to look up median performance scores manually
df_perf_agg <- df_perf %>% 
  group_by( name_caller ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

p_perf <- plot_perf_min( df_perf )
ggsave( file.path( plot_dir, 'spike-in.performance.pdf'), plot = p_perf, width = 12, height = 4)
ggsave( file.path( plot_dir, 'spike-in.performance.png'), plot = p_perf, width = 12, height = 4)

# p_perf <- plot_perf_rrsv( df_perf )
# ggsave( file.path( plot_dir, 'Fig4.spike-in.performance.cvg.pdf'), plot = p_perf, width = 8, height = 10)
# ggsave( file.path( plot_dir, 'Fig4.spike-in.performance.cvg.png'), plot = p_perf, width = 8, height = 10)

# performance by admixture regime
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
df <- df_perf %>% mutate( ttype = fct_recode(ttype, 'med'='medium') )
df$ttype <- factor(df$ttype, levels = c('low', 'med', 'high') )
p_perf_admix <- plot_perf_admix( df )
ggsave( file.path( plot_dir, 'Fig5.spike-in.performance.admix.pdf'), plot = p_perf_admix, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig5.spike-in.performance.admix.png'), plot = p_perf_admix, width = 8, height = 10)

# correlation between F1 score and recall, precision
# ------------------------------------------------------------------------------
cor.test( df$F1, df$recall )
#       Pearson's product-moment correlation
# 
# data:  df$F1 and df$recall
# t = 35.457, df = 2398, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.5595922 0.6121199
# sample estimates:
#   cor 
# 0.5864723 
cor.test( df$F1, df$precision )
#       Pearson's product-moment correlation
# 
# data:  df$F1 and df$precision
# t = 115.25, df = 2398, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9140135 0.9262676
# sample estimates:
#   cor 
# 0.9203662 
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Two-out-of-two performance
# ------------------------------------------------------------------------------

# enumerate all possible pairs of callers
df_caller_pairs <- df_caller  %>% 
  mutate(name_caller = as.character(name_caller)) %>%
  dplyr::filter( !(name_caller %in% c('MuTect1', 'Mutect2_single')) ) %>%
  select( id_caller, name_caller ) %>%
  map_df( function(x) { combn(x, 2, paste, collapse = '_') } )
# calculate intersections of paired call sets
df_varcall_pairs <- df_varcall %>%
  dplyr::filter( !(id_caller %in% c(3, 4)) ) %>%
  inner_join( df_varcall, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
  dplyr::filter( id_caller.x < id_caller.y ) %>%
  unite( id_caller, id_caller.x, id_caller.y)
# assign type (TP, FP, FN) to var calls
df_vars_pairs <- classify_variants( df_varcall_pairs, df_mut, df_mut_sample, df_caller_pairs )
# assign germline status to var calls
df_vars_pairs <- df_vars_pairs %>% 
  left_join( df_snp, by = c( 'chrom', 'pos') ) %>% 
  mutate( germline = (!is.na(id_mut)) ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos, type, germline )

df_perf_pairs <- calculate_performance_sample( df_vars_pairs, df_caller_pairs, df_rep )
df_perf_pairs_agg <- df_perf_pairs %>% group_by( id_caller ) %>%
  summarise( med_F1 = median(F1), med_prec = median(precision), med_rec = median(recall) ) %>%
  separate( id_caller, c('id1', 'id2'), sep = '_') %>%
  mutate( id_caller1 = as.numeric(id1), id_caller2 = as.numeric(id2) ) %>%
  inner_join( df_caller, by = c('id_caller1' = 'id_caller')) %>%
  inner_join( df_caller, by = c('id_caller2' = 'id_caller')) %>%
  select( id_caller1, id_caller2, caller1 = name_caller.x, caller2 = name_caller.y, med_F1, med_prec, med_rec ) %>%
  ungroup()

df_perf_pairs_top <- df_perf_pairs %>%
  inner_join(
    df_perf_pairs %>% 
    group_by( id_caller ) %>%
    summarise( med_F1 = median(F1) ) %>% 
    mutate( rank_mF1 = rank(desc(med_F1)) ) %>%
    dplyr::filter( rank_mF1 <= 10 ) %>% 
    select( id_caller, rank_mF1 ),
    by = c('id_caller')
  )

# plot performance metrics
# ------------------------------------------------------------------------------

p_perf <- plot_perf_pairs( df_perf_pairs_top )
ggsave( file.path( plot_dir, 'spike-in.performance.pairs.pdf'), plot = p_perf, width = 8, height = 10)
ggsave( file.path( plot_dir, 'spike-in.performance.pairs.png'), plot = p_perf, width = 8, height = 10)

# performance by admixture regime
# ------------------------------------------------------------------------------
df <- df_perf_pairs_top %>% mutate( ttype = fct_recode(ttype, 'med'='medium') )
df$ttype <- factor(df$ttype, levels = c('low', 'med', 'high') )
p_perf_admix <- plot_perf_pairs_admix( df )
ggsave( file.path( plot_dir, 'spike-in.performance.pairs.admix.pdf'), plot = p_perf_admix, width = 8, height = 10)
ggsave( file.path( plot_dir, 'spike-in.performance.pairs.admix.png'), plot = p_perf_admix, width = 8, height = 10)

# plot heat map
# ------------------------------------------------------------------------------

p_mF1_hm <- plot_perf_pairs_hm( df_perf_pairs_agg )
ggsave( file.path( plot_dir, 'spike-in.median_F1.pairs.hm.pdf'), plot = p_mF1_hm, width = 8, height = 8)
ggsave( file.path( plot_dir, 'spike-in.median_F1.pairs.hm.png'), plot = p_mF1_hm, width = 8, height = 8)

# ------------------------------------------------------------------------------
# Two-out-of-three performance
# ------------------------------------------------------------------------------

# enumerate all possible triplets of callers
# df_caller_trip <- df_caller %>%
#   select( id_caller, name_caller ) %>%
#   map_df( function(x) { combn(x, 3, paste, collapse = '_') } )
df_caller_trip <- df_caller %>%
  dplyr::filter( !(name_caller %in% c('MuTect1', 'Mutect2_single')) ) %>%
  select( id_caller ) %>%
  map_df( function(x) { combn(x, 3) } ) %>%
  mutate( id_trip = floor((row_number()-1)/3)+1 )
# create accessory table for variants
df_var <- df_varcall %>%
  dplyr::filter( !(id_caller %in% c(3, 4)) ) %>%
  select( id_rep, id_sample, chrom, pos ) %>%
  unique() %>%
  rowid_to_column( 'id_mut' )
saveRDS( df_var, file.path(data_dir, 'df_var.rds') )
# determine presence/absecce of var calls
df_varcall_trip <- df_varcall %>%
  dplyr::filter( !(id_caller %in% c(3, 4)) ) %>%
  inner_join( df_var, by = c('id_rep', 'id_sample', 'chrom', 'pos' ) ) %>%
  select( id_mut, id_caller ) %>%
  inner_join( df_caller_trip )
saveRDS( df_varcall_trip, file.path(data_dir, 'df_varcall_trip.s.rds') )

df_varcall_trip <- readRDS( file.path(data_dir, 'df_varcall_trip.s.rds') )
class(df_varcall_trip$id_mut)

library(data.table)
dt_varcall_trip <- copy( df_varcall_trip )
setDT( dt_varcall_trip )
dt <- dt_varcall_trip[, .(n = .N), by = .(id_mut, id_trip)]
saveRDS( dt, file.path(data_dir, 'dt_varcall_trip.grouped.rds') )

df_caller_trip <- df_caller_trip %>% 
  inner_join( df_caller ) %>%
  group_by( id_trip ) %>%
  summarise( id_caller = str_c(id_caller, collapse = '_'), name_caller = str_c(name_caller, collapse = '_') )
# calculate intersections of paired call sets
df_varcall_trip <- dt %>% dplyr::filter( n > 1 ) %>%
  inner_join( df_var, by = 'id_mut' ) %>% 
  inner_join( df_caller_trip, by = 'id_trip' ) %>% 
  select( -id_mut )
# assign type (TP, FP, FN) to var calls
df_vars_trip <- classify_variants( df_varcall_trip, df_mut, df_mut_sample, df_caller_trip )
saveRDS( df_vars_trip, file.path(data_dir, 'df_vars_trip.rds') )
# assign germline status to var calls
df_vars_trip <- df_vars_trip %>%
  left_join( df_snp, by = c( 'chrom', 'pos') ) %>% 
  mutate( germline = (!is.na(id_mut)) ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos, type, germline )
saveRDS( df_vars_trip, file.path(data_dir, 'df_vars_trip.rds') )

df_perf_trip <- calculate_performance_sample( df_vars_trip, df_caller_trip, df_rep )
saveRDS( df_perf_trip, file.path(data_dir, 'df_perf_trip.rds') )

df_perf_trip <- readRDS( file.path(data_dir, 'df_perf_trip.rds') )

df_perf_trip_agg <- df_perf_trip %>% group_by( id_caller ) %>%
  summarise( med_F1 = median(F1), med_prec = median(precision), med_rec = median(recall) ) %>%
  separate( id_caller, c('id1', 'id2', 'id3') ) %>%
  mutate( id_caller1 = as.numeric(id1), 
          id_caller2 = as.numeric(id2), 
          id_caller3 = as.numeric(id3) ) %>%
  inner_join( df_caller, by = c('id_caller1' = 'id_caller')) %>%
  inner_join( df_caller, by = c('id_caller2' = 'id_caller')) %>%
  inner_join( df_caller, by = c('id_caller3' = 'id_caller')) %>%
  select( id_caller1, id_caller2, id_caller3,
          caller1 = name_caller.x, caller2 = name_caller.y, caller3 = name_caller, 
          med_F1, med_prec, med_rec ) %>%
  ungroup()
df_perf_agg <- df_perf %>% group_by( id_caller ) %>%
  summarise( med_F1 = median(F1), med_prec = median(precision), med_rec = median(recall) )

df_perf_trip_top <- df_perf_trip %>%
  inner_join(
    df_perf_trip %>% 
      group_by( id_caller ) %>%
      summarise( med_F1 = median(F1) ) %>% 
      mutate( rank_mF1 = rank(desc(med_F1)) ) %>%
      dplyr::filter( rank_mF1 <= 10 ) %>% 
      select( id_caller, rank_mF1 ),
    by = c('id_caller')
  )

# plot performance metrics
# ------------------------------------------------------------------------------

p_perf <- plot_perf_pairs( df_perf_trip_top )
ggsave( file.path( plot_dir, 'spike-in.performance.trip.pdf'), plot = p_perf, width = 8, height = 10)
ggsave( file.path( plot_dir, 'spike-in.performance.trip.png'), plot = p_perf, width = 8, height = 10)

p_perf <- plot_perf_trip_hm( df_perf_trip_agg )
ggsave( file.path( plot_dir, 'spike-in.performance.trip.all.pdf'), plot = p_perf, width = 16, height = 14)
ggsave( file.path( plot_dir, 'spike-in.performance.trip.all.png'), plot = p_perf, width = 16, height = 14)


# performance by admixture regime
# ------------------------------------------------------------------------------
df <- df_perf_trip_top %>% mutate( ttype = fct_recode(ttype, 'med'='medium') )
df$ttype <- factor(df$ttype, levels = c('low', 'med', 'high') )
p_perf_admix <- plot_perf_pairs_admix( df )
ggsave( file.path( plot_dir, 'spike-in.performance.trip.admix.pdf'), plot = p_perf_admix, width = 8, height = 10)
ggsave( file.path( plot_dir, 'spike-in.performance.trip.admix.png'), plot = p_perf_admix, width = 8, height = 10)


################################################################################
# Variant allele freqs for TP, FP, FN variants
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

p_vaf <- plot_vaf_bar_srsv( df_vars, df_rc, df_rep )
ggsave( file.path( plot_dir, 'FigS9.spike-in.vaf.bar.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'FigS9.spike-in.vaf.bar.png'), plot = p_vaf, device = png(), width = 8, height = 8 )


################################################################################
# Similarity between call sets
################################################################################

df_var <- readRDS( file.path(data_dir, 'df_vars.rds') )
callerorder = c(
  'Bcftools', 
  'CaVEMan', 
  'MuTect1', 
  'Mutect2_single', 
  'NeuSomatic', 
  'Shimmer', 
  'SNooPer', 
  'SomaticSniper', 
  'Strelka2', 
  'VarDict', 
  'VarScan', 
  'MuClone',
  'SNV-PPILP',
  'HaplotypeCaller', 
  'MultiSNV', 
  'Mutect2_multi_F'
)
## true positives

df_pres_tp <- get_var_pres( df_var, df_caller, 'TP' )
df_jacc_tp <- Jaccard.df( df_pres_tp %>% select(-id_mut)  %>% select(callerorder))
p_jacc_tp <- plot_jacc_idx( df_jacc_tp %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "TP", size = 6)

## false positives
df_pres_fp <- get_var_pres( df_var, df_caller, 'FP' )
df_jacc_fp <- Jaccard.df( df_pres_fp %>% select(-id_mut)  %>% select(callerorder))
p_jacc_fp <- plot_jacc_idx( df_jacc_fp %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "FP", size = 6)

## false negatives
df_pres_fn <- get_var_pres( df_var, df_caller, 'FN' )
df_jacc_fn <- Jaccard.df( df_pres_fn %>% select(-id_mut) %>% select(callerorder))
p_jacc_fn <- plot_jacc_idx( df_jacc_fn %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "FN", size = 6)

## multi-plot
p_jacc_multi <- plot_jacc_idx_multi( p_jacc_tp, p_jacc_fn, p_jacc_fp )
ggsave( file.path( plot_dir, 'spike-in.jaccard.pdf'), plot = p_jacc_multi, device = pdf(), width = 10.5, height = 4 )
ggsave( file.path( plot_dir, 'spike-in.jaccard.png'), plot = p_jacc_multi, device = png(), width = 10.5, height = 4 )



################################################################################
# UpSet plots
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

df_pres <- get_upset_pres( df_vars, df_rc )
df_pres %>% write_csv( file.path(data_dir, 'spike-in.muts_callsets.csv') )

df_pres <- read.csv( file.path(data_dir, 'spike-in.muts_callsets.csv') )
# annoying, but necessary... only if loaded from file
names(df_pres)[names(df_pres)=='SNV.PPILP'] <- 'SNV-PPILP'

# plot variant calls in relation to TRUE somatic variants
df <- df_pres
lbl_callers <- setdiff(df_caller$name_caller, c('Strelka1'))
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS13.spike-in.upset.som')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
# plot_upset_empirical( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()

# plot variant calls intersections (agnostic to truth)

df <- df_pres %>% select(-c(TRUE_somatic, TRUE_germline))
n <- lbl_callers
plot_upset( df, n )  
  
# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'FP' )
lbl_callers <- setdiff(callers$name_caller, c('Strelka1'))
n <- c(lbl_callers, 'TRUE_germline')
fn_pfx <- file.path( plot_dir, 'FigS11.spike-in.upset.FP.GL')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()

# plot TP variant calls in relation to somatic vars
df <- df_pres %>% dplyr::filter( type == 'TP' )
lbl_callers <- setdiff(callers$name_caller, c('Strelka1'))
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS14.spike-in.upset.TP.som')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()

# hierarchical clustering based on varcalls
#-------------------------------------------------------------------------------
require(ade4) # dist.binary()
require(ggdendro) # ggdendrogram()

df_pres <- read.csv( file.path(data_dir, 'spike-in.muts_callsets.csv') )
# annoying, but necessary... only if loaded from file
names(df_pres)[names(df_pres)=='SNV.PPILP'] <- 'SNV-PPILP'

df_pres <- df_pres_tp %>% bind_rows( df_pres_fn ) %>% bind_rows( df_pres_fp )
df_pres <- df_pres %>% select (
  Bcftools, CaVEMan, HaplotypeCaller, MuClone, MultiSNV, MuTect1, Mutect2_multi_F, Mutect2_single, NeuSomatic, 
  Shimmer, SNooPer, `SNV-PPILP`, SomaticSniper, Strelka2, VarDict, VarScan
)
df_jacc <- Jaccard.df( df_pres )
df_jacc_idx <- df_jacc %>% spread( caller1, jaccard_idx ) %>% as.data.frame() 
df_jacc_idx <- df_jacc_idx %>% set_rownames( df_jacc_idx$caller2 ) %>% select( -caller2 ) %>% ungroup()
d <- as.dist( 1-df_jacc_idx )
hc <- hclust(d)

fn_pfx <- 'spike-in.hclust.dendro.tp_fn_fp'
fn_pfx <- 'spike-in.hclust.dendro.df_pres'
pdf( file.path(plot_dir, paste0(fn_pfx, '.pdf')), width = 8, height = 4 )
par( mar = c(2, 0, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()
png( file.path(plot_dir, paste0(fn_pfx, '.png')), width = 8, height = 4, units = 'in', res = 300 )
par( mar = c(2, 1, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()


# calculate tumor-wise performance
#-------------------------------------------------------------------------------

df_vars_tw <- classify_variants_tumorwise( df_varcall, df_mut, df_caller )
df_vars_tw <- df_vars %>% 
  left_join( df_snp, by = c( 'chrom', 'pos') ) %>% 
  mutate( germline = (!is.na(id_mut)) ) %>%
  select( id_caller, id_rep, chrom, pos, type, germline )


df_perf_tw <- calculate_performance_tumorwise( df_vars_tw, df_caller, df_rep )

# plot tumor-wise performance

df_perf_agg <- df_perf_tw %>% 
  group_by( name_caller ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

p_perf <- plot_perf_rrsv( df_perf_tw )
ggsave( file.path( plot_dir, 'spike-in.performance.pdf'), plot = p_perf, width = 12, height = 4)
ggsave( file.path( plot_dir, 'spike-in.performance.png'), plot = p_perf, width = 12, height = 4)
