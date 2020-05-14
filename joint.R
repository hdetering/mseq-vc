#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Load benchmarking data, perform analyses and generate plots for publication.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-01-10
#------------------------------------------------------------------------------

require(tidyverse)
require(ggpubr)
require(rstatix)
require(RColorBrewer)
require(gridExtra)

data_dn = 'data/de-novo'
data_si = 'data/spike-in'
data_em = 'data/empirical'
plot_dir = 'plot/joint'

source( file.path('analysis', 'similarity.analysis.R') )
source( file.path('plotting', 'similarity.plotting.R') )

source( file.path('plotting', 'performance.plotting.R') )

# load variant call sets
#-------------------------------------------------------------------------------

df_pres_dn <- read.csv( file.path(data_dn, 'de-novo.muts_callsets.csv') )
df_pres_si <- read.csv( file.path(data_si, 'spike-in.muts_callsets.csv') )
df_pres_em <- read.csv( file.path(data_em, 'df_pres.csv') )
# rename SNV-PPILP column (necessary if loaded from csv file)
names(df_pres_dn)[names(df_pres_dn)=='SNV.PPILP'] <- 'SNV-PPILP'
names(df_pres_si)[names(df_pres_si)=='SNV.PPILP'] <- 'SNV-PPILP'
names(df_pres_em)[names(df_pres_em)=='SNV.PPILP'] <- 'SNV-PPILP'
#common_cols <- Reduce( intersect, list(names(df_pres_dn), names(df_pres_si), names(df_pres_em)) )

callers = c(
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

df_jacc_dn <- Jaccard.df( df_pres_dn %>% select(sort(callers)) )
df_jacc_si <- Jaccard.df( df_pres_si %>% select(sort(callers)) )
df_jacc_em <- Jaccard.df( df_pres_em %>% select(sort(callers)) )

p_jacc_dn <- plot_jacc_idx( 
  df_jacc_dn %>% 
    mutate(caller1 = factor(caller1, levels = callers), 
           caller2 = factor(caller2, levels = callers)) )
p_jacc_si <- plot_jacc_idx( 
  df_jacc_si %>% 
    mutate(caller1 = factor(caller1, levels = callers), 
           caller2 = factor(caller2, levels = callers)) )
p_jacc_em <- plot_jacc_idx( 
  df_jacc_em %>% 
    mutate(caller1 = factor(caller1, levels = callers), 
           caller2 = factor(caller2, levels = callers)) )

p_jacc_dn <- plot_jacc_idx( df_jacc_dn )
p_jacc_si <- plot_jacc_idx( df_jacc_si )
p_jacc_em <- plot_jacc_idx( df_jacc_em )

# Jaccard index by empirical patient
#-------------------------------------------------------------------------------
df_sample <- read.csv( file.path(data_em, 'df_samples.csv') )
df <- df_pres_em %>% inner_join(df_sample)
df_jacc_em_patient <- df %>%
  group_split( patient ) %>%
  map_dfr( ~ {
    Jaccard.df.lbl(.x %>% select(sort(callers)), .x$patient[1] ) 
  })

df_jacc_em_patient %>%
  group_split( label ) %>%
  map_dfr( ~ {
    tibble( label = .x$label[1] ) %>%
      bind_cols(tidy(cor.test(df_jacc_dn$jaccard_idx, .x$jaccard_idx))
    )
  })
df_jacc_em_patient %>%
  group_split( label ) %>%
  map_dfr( ~ {
    tibble( label = .x$label[1] ) %>%
      bind_cols(tidy(cor.test(df_jacc_si$jaccard_idx, .x$jaccard_idx))
      )
  })

require(ggcorrplot)
df_jacc_em_wide <- df_jacc_em_patient %>%
  pivot_wider( names_from = label, values_from = jaccard_idx )
corr <- df_jacc_em_wide %>% select( -caller1, -caller2 ) %>% cor()
pmat <- df_jacc_em_wide %>% select( -caller1, -caller2 ) %>% cor_pmat()
ggcorrplot(corr, hc.order = TRUE, type = "lower",
          outline.col = "white", lab = TRUE)

# calculate correlation between Jaccard indices for simulated vs. empirical data
#-------------------------------------------------------------------------------

cor.test(df_jacc_dn$jaccard_idx, df_jacc_si$jaccard_idx)
cor.test(df_jacc_dn$jaccard_idx, df_jacc_em$jaccard_idx)
cor.test(df_jacc_si$jaccard_idx, df_jacc_em$jaccard_idx)


# hierarchical clustering based on varcalls
#-------------------------------------------------------------------------------
require(ggdendro) # ggdendrogram()

p_dendro_dn <- plot_dendrogram( df_jacc_dn )
p_dendro_si <- plot_dendrogram( df_jacc_si )
p_dendro_em <- plot_dendrogram( df_jacc_em )

p_jacc_multi <- ggarrange(
  p_jacc_dn, p_dendro_dn,
  p_jacc_si, p_dendro_si,
  p_jacc_em, p_dendro_em,
  labels = "auto",
  ncol = 2,
  nrow = 3,
  common.legend = TRUE,
  legend = "bottom")
ggsave( file.path( plot_dir, 'jacc_idx_dendro.pdf' ), width = 8, height = 12 )
ggsave( file.path( plot_dir, 'jacc_idx_dendro.png' ), width = 8, height = 12 )

fn_pfx <- 'FigS1.de-novo.hclust.dendro'
pdf( file.path(plot_dir, paste0(fn_pfx, '.pdf')), width = 8, height = 4 )
par( mar = c(2, 1, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()
png( file.path(plot_dir, paste0(fn_pfx, '.png')), width = 8, height = 4, units = 'in', res = 300 )
par( mar = c(2, 0, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()



# Jaccard indices for TP, FN and FP calls in simulated datasets
#-------------------------------------------------------------------------------

df_caller_dn  <- readRDS( file.path(data_dn, 'df_caller.rds') )
df_var_dn <- readRDS( file.path(data_dn, 'df_vars.rds') )
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
df_pres_tp_dn <- get_var_pres( df_var_dn, df_caller_dn, 'TP' )
df_jacc_tp_dn <- Jaccard.df( df_pres_tp_dn %>% select(-id_mut)  %>% select(callerorder))

## false negatives
df_pres_fn_dn <- get_var_pres( df_var_dn, df_caller_dn, 'FN' )
df_jacc_fn_dn <- Jaccard.df( df_pres_fn_dn %>% select(-id_mut)  %>% select(callerorder))

## false positives
df_pres_fp_dn <- get_var_pres( df_var_dn, df_caller_dn, 'FP' ) %>%
  add_column( CaVEMan = 0, .after = 2 )
df_jacc_fp_dn <- Jaccard.df( df_pres_fp_dn %>% select(-id_mut)  %>% select(callerorder))

# spike-in
#-------------------------------------------------------------------------------

df_caller_si <- readRDS( file.path(data_si, 'RRSV.callers.rds') )
df_var_si <- readRDS( file.path(data_si, 'df_vars.rds') )
## true positives
df_pres_tp_si <- get_var_pres( df_var_si, df_caller_si, 'TP' )
df_jacc_tp_si <- Jaccard.df( df_pres_tp_si %>% select(-id_mut)  %>% select(callerorder))

## false positives
df_pres_fp_si <- get_var_pres( df_var_si, df_caller_si, 'FP' )
df_jacc_fp_si <- Jaccard.df( df_pres_fp_si %>% select(-id_mut)  %>% select(callerorder))

## false negatives
df_pres_fn_si <- get_var_pres( df_var_si, df_caller_si, 'FN' )
df_jacc_fn_si <- Jaccard.df( df_pres_fn_si %>% select(-id_mut) %>% select(callerorder))

# generate plots
#-------------------------------------------------------------------------------

p_jacc_tp_dn <- plot_jacc_idx( df_jacc_tp_dn %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                        caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "CaVEMan", y = "MultiSNV", label = "TP", size = 6)
p_jacc_fp_dn <- plot_jacc_idx( df_jacc_fp_dn %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                        caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "CaVEMan", y = "MultiSNV", label = "FP", size = 6)

p_jacc_fn_dn <- plot_jacc_idx( df_jacc_fn_dn %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                        caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "CaVEMan", y = "MultiSNV", label = "FN", size = 6)
p_jacc_tp_si <- plot_jacc_idx( df_jacc_tp_si %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                        caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "CaVEMan", y = "MultiSNV", label = "TP", size = 6)
p_jacc_fn_si <- plot_jacc_idx( df_jacc_fn_si %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                        caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "CaVEMan", y = "MultiSNV", label = "FN", size = 6)
p_jacc_fp_si <- plot_jacc_idx( df_jacc_fp_si %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                        caller2 = factor(caller2, levels = callerorder))) +  
  annotate("text", x = "CaVEMan", y = "MultiSNV", label = "FP", size = 6)

## multi-plot
p_jacc <- plot_jacc_idx_joint( p_jacc_tp_dn, p_jacc_fn_dn, p_jacc_fp_dn,
                               p_jacc_tp_si, p_jacc_fn_si, p_jacc_fp_si )
ggsave( file.path( plot_dir, 'FigS10.joint.jaccard.pdf'), plot = p_jacc, device = pdf(), width = 10.5, height = 8.2 )
ggsave( file.path( plot_dir, 'FigS10.joint.jaccard.png'), plot = p_jacc, device = png(), width = 10.5, height = 8.2 )



# PERFORMANCE SIGNIFICANCE TESTING
#===============================================================================

# de-novo sims
#-------------------------------------------------------------------------------
df_perf_dn <- readRDS( file.path(data_dn, 'df_perf.rds') ) %>%
  dplyr::filter( caller %in% callers )  %>% 
  mutate( ttype = fct_recode(ttype, 'high'='us', 'med'='ms', 'low'='hs') )

# Global F1 scores (across conditions)
#--------------------------------------
# one-sided ANOVA (parametric)
aov(F1 ~ caller, data = df_perf_dn) %>% summary()
# Kruskal-Wallis rank sum test (non-parametric)
kruskal.test(F1 ~ caller, data = df_perf_dn)

# F1 scores in relation to factors
#----------------------------------
# two-sided ANOVA (parametric)
dn_anova_f1 <- df_perf_dn %>% group_by( caller ) %>%
  anova_test( F1 ~ cvg * ttype ) %>%
  adjust_pvalue( method = 'bonferroni' ) %>%
  add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
# Kruskal-Wallis rank sum test (non-parametric)
dn_kruskal_f1_cvg <- df_perf_dn %>% group_by( caller ) %>%
  kruskal_test( F1 ~ cvg ) %>%
  adjust_pvalue( method = 'bonferroni' ) %>%
  add_significance( p.col = 'p.adj', output.col = 'sig_sym' ) %>%
  mutate( factor = 'cvg' )
dn_kruskal_f1_ttype <- df_perf_dn %>% group_by( caller ) %>%
  kruskal_test( F1 ~ ttype ) %>%
  adjust_pvalue( method = 'bonferroni' ) %>%
  add_significance( p.col = 'p.adj', output.col = 'sig_sym' ) %>%
  mutate( factor = 'ttype' )

dn_signif <- dn_anova_f1 %>% 
  dplyr::filter( Effect %in% c('cvg', 'ttype') ) %>%
  select( caller, factor = Effect, p_anova = p.adj, sig_anova = sig_sym ) %>%
  left_join( 
    dn_kruskal_f1_cvg %>% 
      bind_rows(dn_kruskal_f1_ttype) %>%
      select( caller, factor, p_kruskal = p.adj, sig_kruskal = sig_sym), 
    by = c('caller', 'factor') )

p_cvg_dn <- plot_perf_cvg_sig( df_perf_dn, 'anova' )
ggsave(plot = p_cvg_dn, filename = 'de-novo.perf-by-cvg.pdf', width = 12, height = 9)
ggsave(plot = p_cvg_dn, filename = 'de-novo.perf-by-cvg.png', width = 12, height = 9)

p_admix_dn <- plot_perf_admix_sig( df_perf_dn, 'anova' )
ggsave(plot = p_admix_dn, filename = 'de-novo.perf-by-admix.pdf', width = 12, height = 9)
ggsave(plot = p_admix_dn, filename = 'de-novo.perf-by-admix.png', width = 12, height = 9)


# spike-in sims
#-------------------------------------------------------------------------------
df_perf_si <- readRDS( file.path(data_si, 'df_perf.rds') ) 

# Global F1 scores (across conditions)
#--------------------------------------
# one-sided ANOVA (parametric)
aov(F1 ~ caller, data = df_perf_si) %>% summary()
# Kruskal-Wallis rank sum test (non-parametric)
kruskal.test(F1 ~ caller, data = df_perf_si)

# F1 scores in relation to factors
#----------------------------------
# two-sided ANOVA (parametric)
aov(F1 ~ caller * ttype, data = df_perf_si) %>% summary()

si_anova_f1 <- df_perf_si %>% group_by( caller ) %>%
  anova_test( F1 ~ ttype ) %>%
  adjust_pvalue( method = 'bonferroni' ) %>%
  add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
# Kruskal-Wallis rank sum test (non-parametric)
si_kruskal_f1 <- df_perf_si %>% group_by( caller ) %>%
  kruskal_test( F1 ~ ttype ) %>%
  adjust_pvalue( method = 'bonferroni' ) %>%
  add_significance( p.col = 'p.adj', output.col = 'sig_sym' ) %>%
  mutate( factor = 'ttype' )

si_signif <- si_anova_f1 %>%
  select( caller, factor = Effect, p_anova = p.adj, sig_anova = sig_sym ) %>%
  left_join( 
    si_kruskal_f1 %>%
      select( caller, factor, p_kruskal = p.adj, sig_kruskal = sig_sym), 
    by = c('caller', 'factor') )

p_perf_si <- plot_perf_admix_sig( df_perf_si, 'anova' )
ggsave(plot = p_perf_si, filename = 'spike-in.perf-by-admix.pdf', width = 12, height = 9)
ggsave(plot = p_perf_si, filename = 'spike-in.perf-by-admix.png', width = 12, height = 9)

# empirical data
#-------------------------------------------------------------------------------

# generate plots
#-------------------------------------------------------------------------------
