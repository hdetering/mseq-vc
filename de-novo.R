#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Load benchmarking data, perform analyses and generate plots for publication.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-05-23
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
data_dir <- file.path( 'data', 'de-novo' )
plot_dir <- file.path( 'plot', 'de-novo' )

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
    'Mutect2_multi'
  ),
  class = c(rep('marginal', 11), rep('two-step', 3), rep('joint', 3))
)
df_caller <- df_caller %>%
  inner_join( callers, by = 'name_caller' )
# do not show these callers in main plots
noshow <- c( 'MuClone_perf' )

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
df_perf <- calculate_performance_sample( df_vars, df_caller, df_rep )
# write summary stats to file
saveRDS( df_perf, file.path(data_dir, 'df_perf.rds') )

# plot performance metrics
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
# to look up median performance scores manually
df_perf_cvg_agg <- df_perf %>% dplyr::filter(!(name_caller %in% noshow)) %>% 
  group_by( name_caller, cvg ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

p_perf <- plot_perf_cvg( df_perf %>% dplyr::filter(!(name_caller %in% noshow)) )
ggsave( file.path( plot_dir, 'Fig2.de-novo.performance.cvg.pdf'), plot = p_perf, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig2.de-novo.performance.cvg.png'), plot = p_perf, width = 8, height = 10)

p_perf <- plot_perf_cvg_aux( df_perf %>% dplyr::filter( name_caller %in% c('MuClone', 'MuClone_perf') ) ) 
ggsave( file.path( plot_dir, 'FigS17.de-novo.performance.cvg.MuClone.pdf'), plot = p_perf, width = 8, height = 10)
ggsave( file.path( plot_dir, 'FigS17.de-novo.performance.cvg.Muclone.png'), plot = p_perf, width = 8, height = 10)


# performance by admixture regime
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
# to look up median performance scores manually
df_perf_mix_agg <- df_perf %>% dplyr::filter(!(name_caller %in% noshow)) %>% 
  group_by( name_caller, ttype ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

df <- df_perf %>% dplyr::filter(!(name_caller %in% noshow)) %>% 
  mutate( ttype = fct_recode(ttype, 'high'='us', 'med'='ms', 'low'='hs') )
p_perf_tt <- plot_perf_admix( df )
ggsave( file.path( plot_dir, 'Fig3.de-novo.performance.admix.pdf'), plot = p_perf_tt, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig3.de-novo.performance.admix.png'), plot = p_perf_tt, width = 8, height = 10)

# check for significant performance differences
# ------------------------------------------------------------------------------

# Kruskal-Wallis rank sum test
kruskal.test(F1 ~ caller, data = df_perf)

p_perf_f1 <- plot_pairwise_wilcoxon( df_perf %>% dplyr::filter(!(name_caller %in% noshow)) )
ggsave( file.path(plot_dir, 'de-novo.f1.pwt.pdf'), plot = p_perf_f1, width = 14, height = 4.5)
ggsave( file.path(plot_dir, 'de-novo.f1.pwt.png'), plot = p_perf_f1, width = 14, height = 4.5)


################################################################################
# Variant allele freqs for TP, FP, FN variants
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

p_vaf <- plot_vaf_bar_srsv( df_vars %>% dplyr::filter(!(name_caller %in% noshow)), df_rc, df_rep )
ggsave( file.path( plot_dir, 'Fig4.de-novo.vaf.bar.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'Fig4.de-novo.vaf.bar.png'), plot = p_vaf, device = png(), width = 8, height = 8 )

# p_vaf <- plot_vaf_bar_ybreak( df_vars, df_rc, df_rep, 'SNooPer', 3000 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.SNooPer.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.SNooPer.png'), plot = p_vaf, device = png(), width = 8, height = 8 )
# p_vaf <- plot_vaf_bar_ybreak( df_vars, df_rc, df_rep, 'VarDict', 3000 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.VarDict.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
# ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.VarDict.png'), plot = p_vaf, device = png(), width = 8, height = 8 )

# ROC curve for TPR (recall, sensitivity) vs. FDR
#-------------------------------------------------------------------------------
df <- df_vars %>% dplyr::filter( !(name_caller %in% noshow) ) %>% 
  inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
  select(caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt) %>%
  mutate( vaf = round((rc_alt)/(rc_ref+rc_alt), 2) ) %>%
  inner_join( df_rep, by = 'id_rep' )
#df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )

roc <- df %>% 
  group_by( caller, vaf, type ) %>% 
  summarize( n = n() ) %>%
  spread( type, n, fill = 0 ) %>%
  arrange( vaf ) %>%
  mutate( TP_FN = sum(TP+FN), TP_FP = sum(TP+FP) ) %>%
  mutate( recall = sum(TP)/sum(TP+FN), 
          precision = sum(TP)/sum(TP+FP) ) %>%
  replace_na( list(precision = 1, recall = 1) )

df_vaf <- df %>% 
  group_by( caller, bin = cut_width(vaf, width = 0.02, boundary = 0 ), type ) %>% 
  summarize( n = n() ) %>%
  spread( type, n, fill = 0 ) %>%
#  mutate( TP_FN = sum(TP+FN), TP_FP = sum(TP+FP) ) %>%
  mutate( recall = sum(TP)/sum(TP+FN), 
          precision = sum(TP)/sum(TP+FP) ) %>%
  replace_na( list(precision = 1, recall = 1) ) %>%
  mutate( F1 = 2*recall*precision/(recall+precision) )

#roc %>% mutate(AUC = sum(diff(FDR) * na.omit(lead(TPR) + TPR)) / 2)
df_vaf %>% select( caller, bin, recall, precision, F1 ) %>%
  gather( measure, score, -caller, -bin) %>%
  mutate( vaf_lims = str_extract_all(bin, '([\\d\\.]+)'), 
          vaf_min = (vaf_lims[[1]]),
          vaf_max = (vaf_lims[[2]])) %>%
  ggplot( aes(bin, score, color = measure) ) +
  geom_point() + labs( x = 'Variant allele frequency' ) +
  scale_x_discrete( breaks = as.character(seq(0, 1, 0.02)) ) +
  facet_wrap( ~ caller )


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

ggsave( file.path( plot_dir, 'FigS3.de-novo.jaccard.pdf'), plot = p_jacc_multi, device = pdf(), width = 10.5, height = 4.2 )
ggsave( file.path( plot_dir, 'FigS3.de-novo.jaccard.png'), plot = p_jacc_multi, device = png(), width = 10.5, height = 4.2 )


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
df <- df_pres
lbl_callers <- setdiff(callers$name_caller, noshow)
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS2.de-novo.upset.som')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()
#system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'FP' )
lbl_callers <- setdiff(callers$name_caller, noshow)
n <- c(lbl_callers, 'TRUE_germline')
fn_pfx <- file.path( plot_dir, 'FigS3.de-novo.upset.FP.GL')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset( df, n )
dev.off()
#system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'TP' )
lbl_callers <- setdiff(callers$name_caller, noshow)
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS4.de-novo.upset.TP.som')
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
