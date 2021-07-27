data_dir  = "/Users/laura/GoogleDrive/LAB FOLDERS/M-seq Variant Calling Benchmarking/de-novo/data"#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Load benchmarking data, perform analyses and generate plots for publication.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-08-01
#------------------------------------------------------------------------------

require(tidyverse)
require(dbplyr)
require(RSQLite)
require(cowplot)
require(ggpubr)
require(RColorBrewer)
require(gridExtra)
require(grid)
require(broom) # tidy()
require(rstatix) # wilcox_test()
require(stringr)

# source required scripts
source( file.path('analysis', 'performance.analysis.R') )
source( file.path('analysis', 'admixture.analysis.R') )
source( file.path('analysis', 'similarity.analysis.R') )
source( file.path('analysis', 'upset.analysis.R') )
source( file.path('analysis', 'af.denovo.analysis.R') )
source( file.path('plotting', 'performance.plotting.R') )
source( file.path('plotting', 'rc_alt.plotting.R') )
source( file.path('plotting', 'similarity.plotting.R') )
source( file.path('plotting', 'upset.plotting.R') )

# define input/output paths
data_dir <- file.path( 'data', 'de-novo' )
plot_dir <- file.path( 'plot', 'de-novo' )

# connection to analysis database
#-------------------------------------------------------------------------------
# db <- file.path( data_dir, 'analysis.db' )
# con <- DBI::dbConnect(RSQLite::SQLite(), db)
# df_rep <- tbl( con, 'replicates' ) %>% collect()
# df_mut <- tbl( con, 'mutations' ) %>% collect()
# df_mut_sample <- tbl( con, 'mutations_samples' ) %>% collect()
# df_mut_clone <- tbl( con, 'clones_mutations' ) %>% collect()
# df_prev <- tbl( con, 'clones_prev' ) %>% collect()
# df_caller <- tbl( con, 'callers' ) %>% collect()
# df_varcall <- tbl( con, 'varcalls' ) %>% collect()
# df_rc <- tbl( con, 'readcounts' ) %>% collect()
# df_snp <- tbl( con, 'snps' ) %>% collect()


df_rep <- readRDS(paste0(data_dir, "/df_rep.rds"))
df_mut <- readRDS(paste0(data_dir, "/df_mut.rds"))
df_mut_sample <- readRDS(paste0(data_dir, "/df_mut_sample.rds"))
df_mut_clone <- readRDS(paste0(data_dir, "/df_mut_clone.rds"))
df_prev <- readRDS(paste0(data_dir, "/df_prev.rds"))
df_caller <- readRDS(paste0(data_dir, "/df_caller.rds"))
df_varcall <- readRDS(paste0(data_dir, "/df_varcall.rds"))
df_rc <- readRDS(paste0(data_dir, "/df_rc.rds"))
df_snp <- readRDS(paste0(data_dir, "/df_snp.rds"))

# rename Mutect2 sub-modes
# df_caller <- df_caller %>% 
#   mutate( name_caller = str_replace(name_caller, "Mutect2$", "Mutect2_single") ) %>%
#   mutate( name_caller = str_replace(name_caller, "Mutect2_mseq$", "Mutect2_multi") )

# define variant calling tools and their properties
#-------------------------------------------------------------------------------
callers <- tibble(
  name_caller = c(
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
    'MuClone_perf',
    'SNV-PPILP',
    'HaplotypeCaller', 
    'MultiSNV', 
    'Mutect2_multi_F'
  ),
  class = c(rep('marginal', 11), rep('two-step', 3), rep('joint', 3))
)
df_caller <- df_caller %>%
  inner_join( callers, by = 'name_caller' )
# do not show these callers in main plots
noshow <- c( 'MuClone_perf' )

# store data frames on file system
# df_caller %>% saveRDS( file.path(data_dir, 'df_caller.rds') )
# df_rep %>% saveRDS( file.path(data_dir, 'df_rep.rds') )
# df_mut %>% saveRDS( file.path(data_dir, 'df_mut.rds') )
# df_mut_sample %>% saveRDS( file.path(data_dir, 'df_mut_sample.rds') )
# df_mut_clone %>% saveRDS( file.path(data_dir, 'df_mut_clone.rds') )
# df_prev %>% saveRDS( file.path(data_dir, 'df_prev.rds') )
# df_varcall %>% saveRDS( file.path(data_dir, 'df_varcall.rds') )
# df_rc %>% saveRDS( file.path(data_dir, 'df_rc.rds') )
# df_snp %>% saveRDS( file.path(data_dir, 'df_snp.rds') )

################################################################################
# Performance metrics (recall, precision, F1 score)
################################################################################

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

df_perf_sample <- calculate_performance_sample( df_vars, df_caller, df_rep )
df_perf_freq <- calculate_performance_freq( df_vars, df_caller, df_rep, df_rc ) 

# write summary stats to file
saveRDS( df_perf_sample, file.path(data_dir, 'df_perf_sample.rds') )
saveRDS( df_perf_freq, file.path(data_dir, 'df_perf_freq.rds') )


# plot performance metrics
# ------------------------------------------------------------------------------
df_perf_sample <- readRDS( file.path(data_dir, 'df_perf_sample.rds') )
df_perf_freq <- readRDS( file.path(data_dir, 'df_perf_freq.rds') )

# to look up median performance scores manually
df_perf_sample_agg <- df_perf_sample %>% dplyr::filter(!(name_caller %in% noshow)) %>% 
  group_by( name_caller ) %>% 
  dplyr::summarise( 
    med_rec = median(recall, na.rm = T),
    med_pre = median(precision, na.rm = T),
    med_F1  = median(F1, na.rm = T),
    avg_rec = mean(recall, na.rm = T),
    avg_pre = mean(precision, na.rm = T),
    avg_F1  = mean(F1, na.rm = T))

df_perf_freq_agg <- df_perf_freq %>%
  group_by( name_caller ) %>%
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

p_perf_sample <- plot_perf_min( df_perf_sample %>% dplyr::filter(!(name_caller %in% noshow)) )
ggsave( file.path( plot_dir, 'de-novo.performance.sample.pdf'), plot = p_perf_sample, width = 12, height = 4)
ggsave( file.path( plot_dir, 'de-novo.performance.sample.png'), plot = p_perf_sample, width = 12, height = 4)

p_perf_freq <- plot_perf_freq( df_perf_freq )
ggsave( file.path( plot_dir, 'de-novo.performance.freq.pdf'), plot = p_perf_freq, width = 12, height = 12)
ggsave( file.path( plot_dir, 'de-novo.performance.freq.png'), plot = p_perf_freq, width = 12, height = 12)

# plot F1 score histograms to check if scores are normally-distributed
#df_perf %>% ggplot( aes(x = F1) ) + geom_histogram() + facet_wrap( ~name_caller, ncol = 1 )
#ggsave( filename = 'de-novo.F1.hist.pdf', width = 4, height = 20 )


# determine status of variant calls for the per-tumor performance 
#   TP: true positives
#   FP: false positives
#   FN: false negatives
#-------------------------------------------------------------------------------
df_vars_tumor <- classify_variants_pertumor( df_varcall, df_mut, df_mut_sample, df_caller )
df_vars_tumor <- df_vars_tumor %>% 
  left_join( df_snp ) %>% 
  mutate( germline = (!is.na(id_mut)) ) %>% 
  select( id_caller, id_rep, chrom, pos, type, germline )
# store variant calls for later use
saveRDS( df_vars_tumor, file.path(data_dir, 'df_vars_tumor.rds') )

# calculate per-tumor performance metrics
# ------------------------------------------------------------------------------
df_vars_tumor <- readRDS( file.path(data_dir, 'df_vars_tumor.rds') )
df_perf_tumor <- calculate_performance_tumor( df_vars_tumor, df_caller, df_rep )

# write summary stats to file
saveRDS( df_perf_tumor, file.path(data_dir, 'df_perf_tumor.rds') )

# plot per-tumor performance metrics
# ------------------------------------------------------------------------------
df_perf_tumor <- readRDS( file.path(data_dir, 'df_perf_tumor.rds') )
# to look up median performance scores manually
df_perf_tumor_agg <- df_perf_tumor %>% 
  group_by( name_caller ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

p_perf_tumor <- plot_perf_min_sig( df_perf_tumor %>% dplyr::filter(!(name_caller %in% noshow)) )
ggsave( file.path( plot_dir, 'de-novo.performance.tumor.pdf'), plot = p_perf_tumor, width = 12, height = 4)
ggsave( file.path( plot_dir, 'de-novo.performance.tumor.png'), plot = p_perf_tumor, width = 12, height = 4)


# performance by coverage
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
# to look up median performance scores manually
df_perf_cvg_agg <- df_perf %>% dplyr::filter(!(name_caller %in% noshow)) %>% 
  group_by( name_caller, cvg ) %>% 
  summarise( 
    med_rec = median(recall, na.rm = T),
    med_pre = median(precision, na.rm = T),
    med_F1  = median(F1, na.rm = T),
    avg_rec = mean(recall, na.rm = T),
    avg_pre = mean(precision, na.rm = T),
    avg_F1  = mean(F1, na.rm = T))

p_perf <- plot_perf_cvg_sig( df_perf %>% dplyr::filter(!(name_caller %in% noshow)) )
ggsave( file.path( plot_dir, 'de-novo.performance.cvg.mean.pdf'), plot = p_perf, width = 8, height = 10)
# convert PDF to PNG (R png device does not support fonts)
# command works on Linux (MacOS not tested)
system(paste('convert -density 300',
             file.path(plot_dir, 'de-novo.performance.cvg.mean.pdf'),
             '-quality 90', file.path(plot_dir, 'de-novo.performance.cvg.mean.png')))

ggsave( file.path( plot_dir, 'de-novo.performance.cvg.pdf'), plot = p_perf, width = 8, height = 10)
# convert PDF to PNG (R png device does not support fonts)
# command works on Linux (MacOS not tested)
system(paste('convert -density 300',
             file.path(plot_dir, 'de-novo.performance.cvg.pdf'),
             '-quality 90', file.path(plot_dir, 'de-novo.performance.cvg.png')))
#ggsave( file.path( plot_dir, 'de-novo.performance.cvg.png'), plot = p_perf, width = 8, height = 10)

# performance by admixture regime
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
# to look up median performance scores manually
df_perf_mix_agg <- df_perf %>% dplyr::filter(!(name_caller %in% noshow)) %>% 
  group_by( name_caller, ttype ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

df <- df_perf %>% dplyr::filter(!(name_caller %in% noshow)) %>% 
  mutate( ttype = fct_recode(ttype, 'high'='us', 'med'='ms', 'low'='hs') )
p_perf_tt <- plot_perf_admix_sig( df )

ggsave( file.path( plot_dir, 'de-novo.performance.admix.mean.pdf'), plot = p_perf_tt, width = 8, height = 10)
# convert PDF to PNG (R png device does not support fonts)
# command works on Linux (MacOS not tested)
system(paste('convert -density 300',
             file.path(plot_dir, 'de-novo.performance.admix.mean.pdf'),
             '-quality 90', file.path(plot_dir, 'de-novo.performance.admix.mean.png')))
#ggsave( file.path( plot_dir, 'de-novo.performance.admix.png'), plot = p_perf_tt, width = 8, height = 10)

ggsave( file.path( plot_dir, 'de-novo.performance.admix.pdf'), plot = p_perf_tt, width = 8, height = 10)
# convert PDF to PNG (R png device does not support fonts)
# command works on Linux (MacOS not tested)
system(paste('convert -density 300',
             file.path(plot_dir, 'de-novo.performance.admix.pdf'),
             '-quality 90', file.path(plot_dir, 'de-novo.performance.admix.png')))
#ggsave( file.path( plot_dir, 'de-novo.performance.admix.png'), plot = p_perf_tt, width = 8, height = 10)

# correlation between F1 score and recall, precision
# ------------------------------------------------------------------------------
cor.test( df$F1, df$recall )
#       Pearson's product-moment correlation
# 
# data:  df$F1 and df$recall
# t = 132.95, df = 9598, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.7978892 0.8119720
# sample estimates:
#   cor 
# 0.805044 
cor.test( df$F1, df$precision )
#       Pearson's product-moment correlation
# 
# data:  df$F1 and df$precision
# t = 92.23, df = 9598, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.6747043 0.6959187
# sample estimates:
#   cor 
# 0.685457 
# ------------------------------------------------------------------------------


# check for significant performance differences
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') ) %>%
  dplyr::filter( !(caller %in% noshow) )

p_perf_f1 <- plot_pairwise_wilcoxon( df_perf %>% dplyr::filter(!(name_caller %in% noshow)) )
ggsave( file.path(plot_dir, 'de-novo.f1.pwt.pdf'), plot = p_perf_f1, width = 8, height = 10)
ggsave( file.path(plot_dir, 'de-novo.f1.pwt.png'), plot = p_perf_f1, width = 8, height = 10)

# Global F1 scores (across conditions)
#--------------------------------------
# one-sided ANOVA (parametric)
aov(F1 ~ caller, data = df_perf) %>% summary()
# Kruskal-Wallis rank sum test (non-parametric)
kruskal.test(F1 ~ caller, data = df_perf)

# F1 scores interaction with factors
#------------------------------------
# two-sided ANOVA (parametric)
aov(F1 ~ caller * cvg, data = df_perf) %>% summary()
aov(F1 ~ caller * ttype, data = df_perf) %>% summary()
# Kruskal-Wallis rank sum test (non-parametric)
df_perf %>%
  group_by( name_caller ) %>%
  kruskal_test( F1 ~ cvg ) %>%
  adjust_pvalue( method = 'BH' )
df_perf %>%
  group_by( name_caller ) %>%
  kruskal_test( F1 ~ ttype ) %>%
  adjust_pvalue( method = 'BH' )

p_cvg_sig <- plot_perf_cvg_sig( df_perf, 'anova' )
ggsave( plot = p_cvg_sig, filename = 'de-novo.perf.cvg.anova.pdf', device = 'pdf', width = 12, height = 10 )
ggsave( plot = p_cvg_sig, filename = 'de-novo.perf.cvg.anova.png', device = 'png', width = 12, height = 10 )

p_cvg_sig <- plot_perf_cvg_sig( df_perf, 'kruskal.test' )
ggsave( plot = p_cvg_sig, filename = 'de-novo.perf.cvg.kruskal.pdf', device = 'pdf', width = 12, height = 10 )
ggsave( plot = p_cvg_sig, filename = 'de-novo.perf.cvg.kruskal.png', device = 'png', width = 12, height = 10 )

################################################################################
# Variant allele freqs for TP, FP, FN variants
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

p_vaf <- plot_vaf_bar_srsv( df_vars %>% dplyr::filter(!(name_caller %in% noshow)), df_rc, df_rep )
ggsave( file.path( plot_dir, 'FigS8.de-novo.vaf.bar.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'FigS8.de-novo.vaf.bar.png'), plot = p_vaf, device = png(), width = 8, height = 8 )

# p_vaf <- plot_vaf_bar_ybreak( df_vars, df_rc, df_rep, 'SNooPer', 3000 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.SNooPer.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.SNooPer.png'), plot = p_vaf, device = png(), width = 8, height = 8 )
# p_vaf <- plot_vaf_bar_ybreak( df_vars, df_rc, df_rep, 'VarDict', 3000 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.VarDict.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.VarDict.png'), plot = p_vaf, device = png(), width = 8, height = 8 )


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
  'Mutect2_multi'
)
## true positives
df_pres_tp <- get_var_pres( df_var, df_caller, 'TP' )
df_jacc_tp <- Jaccard.df( df_pres_tp %>% select(-id_mut)  %>% select(callerorder))
p_jacc_tp <- plot_jacc_idx( df_jacc_tp %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "TP", size = 6)

## false negatives
df_pres_fn <- get_var_pres( df_var, df_caller, 'FN' )
df_jacc_fn <- Jaccard.df( df_pres_fn %>% select(-id_mut)  %>% select(callerorder))
p_jacc_fn <- plot_jacc_idx( df_jacc_fn %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "FN", size = 6)

## false positives
df_pres_fp <- get_var_pres( df_var, df_caller, 'FP' ) %>%
  add_column( CaVEMan = 0, .after = 2 )
df_jacc_fp <- Jaccard.df( df_pres_fp %>% select(-id_mut)  %>% select(callerorder))
p_jacc_fp <- plot_jacc_idx( df_jacc_fp %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "FP", size = 6)

## multi-plot
p_jacc_multi <- plot_jacc_idx_multi( p_jacc_tp, p_jacc_fn, p_jacc_fp )

ggsave( file.path( plot_dir, 'de-novo.jaccard.pdf'), plot = p_jacc_multi, device = pdf(), width = 10.5, height = 4.2 )
ggsave( file.path( plot_dir, 'de-novo.jaccard.png'), plot = p_jacc_multi, device = png(), width = 10.5, height = 4.2 )


# hierarchical clustering based on varcalls
#-------------------------------------------------------------------------------
require(ade4) # dist.binary()
require(ggdendro) # ggdendrogram()
df_pres <- df_pres_tp %>% bind_rows( df_pres_fn ) %>% bind_rows( df_pres_fp )
df_jacc <- Jaccard.df( df_pres %>% select(-id_mut, -MuClone_perf) )

df_jacc_idx <- df_jacc %>% spread( caller1, jaccard_idx ) %>% as.data.frame() 
df_jacc_idx <- df_jacc_idx %>% set_rownames( df_jacc_idx$caller2 ) %>% select( -caller2 )
d <- as.dist( 1-df_jacc_idx )
hc <- hclust(d)

fn_pfx <- 'FigS1.de-novo.hclust.dendro'
pdf( file.path(plot_dir, paste0(fn_pfx, '.pdf')), width = 8, height = 4 )
par( mar = c(2, 1, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()
png( file.path(plot_dir, paste0(fn_pfx, '.png')), width = 8, height = 4, units = 'in', res = 300 )
par( mar = c(2, 0, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()


################################################################################
# UpSet plots
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

df_pres <- get_upset_pres( df_vars, df_rc )
df_pres %>% write_csv( file.path(data_dir, 'de-novo.muts_callsets.csv') )

df_pres <- read.csv( file.path(data_dir, 'de-novo.muts_callsets.csv') )
# annoying, but necessary... only if loaded from file
names(df_pres)[names(df_pres)=='SNV.PPILP'] <- 'SNV-PPILP'

# plot variant calls in relation to TRUE somatic variants
df <- df_pres %>% as.data.frame()
lbl_callers <- setdiff(callers$name_caller, noshow)
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS12.de-novo.upset.som')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()
#system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'FP' ) %>% as.data.frame()
lbl_callers <- setdiff(callers$name_caller, noshow)
n <- c(lbl_callers, 'TRUE_germline')
fn_pfx <- file.path( plot_dir, 'FigS9.de-novo.upset.FP.GL')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()
#system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'TP' ) %>% as.data.frame()
lbl_callers <- setdiff(callers$name_caller, noshow)
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS13.de-novo.upset.TP.som')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()
#system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )


################################################################################
# Admixture (metric: Shannon Entropy of prevalence values)
################################################################################

df_shannon_sample <- get_shannon_sample( df_prev ) 
df_shannon_rep <- df_shannon_sample %>%
  group_by( id_rep ) %>% 
  summarise( H_bar = mean(H) ) %>% 
  ungroup()

# Admixture of tumors by simulation regimes
#-------------------------------------------------------------------------------

df <- df_shannon_rep %>%
  inner_join(df_rep, by = 'id_rep')
p_ttype_shannon <- ggplot(df, aes(x = ttype, y = H_bar)) + 
  geom_boxplot() +
  scale_x_discrete(labels=c('low admixture', 'medium admixture', 'high admixture')) +
  coord_flip() +
  #ggtitle('Mean Shannon entropy by scenario type') +
  theme(axis.title.y = element_blank()) +
  labs(y = 'Mean Shannon entropy')
ggsave( file.path( plot_dir, 'SRSV.ttype_shannon.pdf'), plot = p_ttype_shannon, device = pdf(), height = 2, width = 6 )
ggsave( file.path( plot_dir, 'SRSV.ttype_shannon.png'), plot = p_ttype_shannon, device = png(), height = 2, width = 6 )



# Effect of Admixture on FN
#-------------------------------------------------------------------------------

df_perf_sample <- calculate_performance_sample( df_vars, df_caller, df_rep ) %>%
  dplyr::filter( name_caller != 'Strelka1' )
df <- df_perf_sample %>%
  inner_join( df_shannon_sample, by = c('id_rep', 'id_sample') )
df$name_caller <- factor( df$name_caller, levels = callers$name_caller )
df$class <- factor( df$class, levels = c('marginal', 'joint', 'two-step') )

df_cor_fn <- df %>% nest(-name_caller, -cvg) %>%
  mutate( test = map(data, ~ cor.test(.x$FN, .x$H)),
          tidied = map(test, tidy) ) %>%
  select( name_caller, cvg, tidied ) %>%
  unnest() %>%
  mutate( cor = sprintf('%.2f (%.2g)', estimate, p.value) ) %>%
  select( name_caller, cvg, cor ) %>%
  spread(cvg, cor)
write_csv( df_cor_fn, file.path(data_dir, 'df_cor_fn.csv') )

p_f1_shannon_sample <- ggplot(df, aes(y = F1, x = H)) +
  geom_point(aes(color = class), size = 0.5) +
  geom_label(aes(x = 0.5, y = 0.5, label = round(c, 2)), alpha = 0.5,
             data = df %>% group_by(name_caller, cvg) %>% summarise(c = cor(F1, H)) ) +
  facet_grid(name_caller ~ cvg) +
  labs( x = 'Shannon Entropy' ) + 
  theme_grey() + 
  theme( legend.position = 'top' )
ggsave( file.path( plot_dir, 'SRSV.correlation_F1_ShannonEntropy_sample.pdf'), 
        plot = p_f1_shannon_sample, device = pdf(), width = 4, height = 15 )
ggsave( file.path( plot_dir, 'SRSV.correlation_F1_ShannonEntropy_sample.png'), 
        plot = p_f1_shannon_sample, device = png(), width = 4, height = 15 )

p_rec_shannon_sample <- ggplot(df, aes(y = recall, x = H)) +
  geom_point(aes(color = class), size = 0.5) +
  geom_label(aes(x = 0.5, y = 0.5, label = round(c, 2)), alpha = 0.5,
             data = df %>% group_by(name_caller, cvg) %>% summarise(c = cor(recall, H)) ) +
  facet_grid(name_caller ~ cvg) +
  labs( x = 'Shannon Entropy' ) + 
  theme_grey() + 
  theme( legend.position = 'top' )
ggsave( file.path( plot_dir, 'SRSV.correlation_recall_ShannonEntropy_sample.pdf'), 
        plot = p_rec_shannon_sample, device = pdf(), width = 4, height = 15 )
ggsave( file.path( plot_dir, 'SRSV.correlation_recall_ShannonEntropy_sample.png'), 
        plot = p_rec_shannon_sample, device = png(), width = 4, height = 15 )

p_pre_shannon_sample <- ggplot(df, aes(y = precision, x = H)) +
  geom_point(aes(color = class), size = 0.5) +
  geom_label(aes(x = 0.5, y = 0.5, label = round(c, 2)), alpha = 0.5,
             data = df %>% group_by(name_caller, cvg) %>% summarise(c = cor(precision, H)) ) +
  facet_grid(name_caller ~ cvg) +
  labs( x = 'Shannon Entropy' ) + 
  theme_grey() + 
  theme( legend.position = 'top' )
ggsave( file.path( plot_dir, 'SRSV.correlation_precision_ShannonEntropy_sample.pdf'), 
        plot = p_pre_shannon_sample, device = pdf(), width = 4, height = 15 )
ggsave( file.path( plot_dir, 'SRSV.correlation_precision_ShannonEntropy_sample.png'), 
        plot = p_pre_shannon_sample, device = png(), width = 4, height = 15 )

p_fn_shannon_sample <- ggplot( df, aes(y = FN, x = H) ) +
  geom_point( aes(color = class), size = 0.5 ) +
  geom_label( aes(x = 0.5, y = 75, label = round(c, 2)), alpha = 0.5,
              data = df %>% group_by(name_caller, cvg) %>% summarise(c = cor(FN, H)) ) +
  facet_grid( name_caller ~ cvg ) +
  labs( x = 'Shannon Entropy' ) + 
  theme_grey() + 
  theme( strip.text.y = element_text(angle = 0) ) + 
  theme( legend.position = 'top' )
ggsave( file.path( plot_dir, 'Fig2.SRSV.correlation_FN_ShannonEntropy_sample.pdf'), 
        plot = p_fn_shannon_sample, device = pdf(), width = 6, height = 15 )
ggsave( file.path( plot_dir, 'Fig2.SRSV.correlation_FN_ShannonEntropy_sample.png'), 
        plot = p_fn_shannon_sample, device = png(), width = 6, height = 15 )
