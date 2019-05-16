#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Load benchmarking data, perform analyses and generate plots for publication.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-05-05
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
df_caller <- df_caller %>% mutate( name_caller = fct_recode( name_caller, 'Mutect2_ms' = 'Mutect2_mseq' ) )
df_caller <- df_caller %>%
  inner_join( callers, by = 'name_caller' )

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
# to look up median performance scores manually
df_perf_agg <- df_perf %>% 
  group_by( name_caller, cvg ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

# plot performance metrics
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
# to look up median performance scores manually
df_perf_cvg_agg <- df_perf %>% 
  group_by( name_caller, cvg ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

p_perf <- plot_perf_cvg( df_perf )
ggsave( file.path( plot_dir, 'Fig2.de-novo.performance.cvg.pdf'), plot = p_perf, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig2.de-novo.performance.cvg.png'), plot = p_perf, width = 8, height = 10)

# performance by admixture regime
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
# to look up median performance scores manually
df_perf_mix_agg <- df_perf %>% 
  group_by( name_caller, ttype ) %>% 
  summarise( med_rec = median(recall), med_pre = median(precision), med_F1 = median(F1) )

df <- df_perf %>% mutate( ttype = fct_recode(ttype, 'high'='us', 'med'='ms', 'low'='hs') )
p_perf_tt <- plot_perf_admix( df )
ggsave( file.path( plot_dir, 'Fig3.de-novo.performance.admix.pdf'), plot = p_perf_tt, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig3.de-novo.performance.admix.png'), plot = p_perf_tt, width = 8, height = 10)

# check for significant performance differences
# ------------------------------------------------------------------------------
# plot Box- and Violinplots of F1 scores for callers
p_perf_f1_box <- ggviolin(df_perf, x = "caller", y = "F1", color = "class",
                         main = "A) Distribution of F1 scores",
                         add = "boxplot") +
  rotate_x_text(45)

# Kruskal-Wallis rank sum test
kruskal.test(F1 ~ caller, data = df_perf)

# perform pairwise Wilcoxon rank sum test (Mann-Whitney test?)
pwt <- pairwise.wilcox.test( df_perf$F1, df_perf$caller, p.adjust.method = "BH" )
df_pwt <- pwt$p.value %>% 
  as_tibble(rownames = "id1") %>% 
  gather(id2, p.adj, -id1) %>% dplyr::filter(!is.na(p.adj)) %>%
  mutate(significance = cut(p.adj, 
                            breaks = c(0.0, 0.0001, 0.001, 0.01, 0.05, 1), 
                            labels = c("****", "***", "**", "*", "ns")))
lst_callers <- callers
#df_pwt$id1 <- factor(df_pwt$id1, levels = lst_callers)
#df_pwt$id2 <- factor(df_pwt$id2, levels = lst_callers)
p_perf_f1_pwt <- ggplot(df_pwt, aes(x = id1, y = id2)) + 
  geom_tile(aes(fill = significance), color = "white") +
  geom_text(aes(label = round(-log(p.adj)))) +
  ggtitle("B) Pairwise Wilcoxon rank sum test (values: -log(p.adj))") + 
  scale_fill_manual(values = c(rev(brewer.pal(4, "Greens")), "grey")) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


p_perf_f1 <- grid.arrange(grobs = list(p_perf_f1_box, p_perf_f1_pwt), 
                         layout_matrix = rbind(c(1,2), c(1,2)))
ggsave( file.path(plot_dir, 'de-novo.f1.pwt.pdf'), plot = p_perf_f1, width = 14, height = 4.5)
ggsave( file.path(plot_dir, 'de-novo.f1.pwt.png'), plot = p_perf_f1, width = 14, height = 4.5)


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
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.png'), plot = p_vaf, device = png(), width = 8, height = 8 )

p_vaf <- plot_vaf_bar_ybreak( df_vars, df_rc, df_rep, 'Mutect2_mseq', 3000 )
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.Mutect2_ms.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.Mutect2_ms.png'), plot = p_vaf, device = png(), width = 8, height = 8 )
p_vaf <- plot_vaf_bar_ybreak( df_vars, df_rc, df_rep, 'SNooPer', 3000 )
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.SNooPer.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.SNooPer.png'), plot = p_vaf, device = png(), width = 8, height = 8 )
p_vaf <- plot_vaf_bar_ybreak( df_vars, df_rc, df_rep, 'VarDict', 3000 )
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.VarDict.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'Fig3.de-novo.vaf.bar.VarDict.png'), plot = p_vaf, device = png(), width = 8, height = 8 )


################################################################################
# Similarity between call sets
################################################################################

df_var <- readRDS( file.path(data_dir, 'df_vars.rds') )
callerorder = c(
  'Bcftools', 
  'CaVEMan', 
  'MuTect1', 
  'Mutect2', 
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
  'Mutect2_ms'
)
## true positives
df_pres_tp <- get_var_pres( df_var, df_caller, 'TP' )
df_jacc_tp <- Jaccard.df( df_pres_tp %>% select(-id_mut)  %>% select(callerorder))
p_jacc_tp <- plot_jacc_idx( df_jacc_tp %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "TP", size = 10)

## false negatives
df_pres_fn <- get_var_pres( df_var, df_caller, 'FN' )
df_jacc_fn <- Jaccard.df( df_pres_fn %>% select(-id_mut)  %>% select(callerorder))
p_jacc_fn <- plot_jacc_idx( df_jacc_fn %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "FN", size = 10)

## false positives
df_pres_fp <- get_var_pres( df_var, df_caller, 'FP' ) %>%
  add_column( CaVEMan = 0, .after = 2 )
df_jacc_fp <- Jaccard.df( df_pres_fp %>% select(-id_mut)  %>% select(callerorder))
p_jacc_fp <- plot_jacc_idx( df_jacc_fp %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "MuTect1", y = "MultiSNV", label = "FP", size = 10)

## multi-plot
p_jacc_multi <- plot_jacc_idx_multi( p_jacc_tp, p_jacc_fn, p_jacc_fp )
ggsave( file.path( plot_dir, 'Fig4.de-novo.jaccard.pdf'), plot = p_jacc_multi, device = pdf(), width = 10, height = 4.5 )
ggsave( file.path( plot_dir, 'Fig4.de-novo.jaccard.png'), plot = p_jacc_multi, device = png(), width = 10, height = 4.5 )

# hierarchical clustering based on varcalls
#-------------------------------------------------------------------------------
require(ade4) # dist.binary()
require(ggdendro) # ggdendrogram()
df_pres <- df_pres_tp %>% bind_rows( df_pres_fn ) %>% bind_rows( df_pres_fp )
df_jacc <- Jaccard.df( df_pres %>% select(-id_mut, -Strelka1) )
df_jacc_idx <- df_jacc %>% 
  spread( caller1, jaccard_idx ) %>% 
  as.data.frame() %>% set_rownames( df_jacc_dist$caller2 ) %>% select( -caller2 )
d <- as.dist( 1-df_jacc_idx )
hc <- hclust(d)

fn_pfx <- 'FigS1.de-novo.hclust.dendro'
pdf( file.path(plot_dir, paste0(fn_pfx, '.pdf')), width = 8, height = 4 )
par( mar = c(2, 0, 0, 6) )
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
lbl_callers <- setdiff(callers$name_caller, c('Strelka1'))
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS2.de-novo.upset.som')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'FP' )
lbl_callers <- setdiff(callers$name_caller, c('Strelka1'))
n <- c(lbl_callers, 'TRUE_germline')
fn_pfx <- file.path( plot_dir, 'FigS3.de-novo.upset.FP.GL')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )

# plot FP variant calls in relation to germline vars
df <- df_pres %>% dplyr::filter( type == 'TP' )
lbl_callers <- setdiff(callers$name_caller, c('Strelka1'))
n <- c(lbl_callers, 'TRUE_somatic')
fn_pfx <- file.path( plot_dir, 'FigS4.de-novo.upset.TP.som')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset( df, n )
dev.off()
system( sprintf('pdftoppm %s.pdf %s -png', fn_pfx, fn_pfx) )


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