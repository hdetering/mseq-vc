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

# define input/output paths
data_dir <- file.path( 'data', 'SRSV' )
plot_dir <- file.path( 'plot', 'SRSV' )

# source required scripts
source( file.path('analysis', 'performance.analysis.R') )
source( file.path('plotting', 'performance.plotting.R') )
source( file.path('analysis', 'admixture.analysis.R') )
source( file.path('plotting', 'rc_alt.plotting.R') )
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
    'HaplotypeCaller', 
    'MultiSNV', 
    'Mutect2_mseq',
    'MuClone', 
    'SNV-PPILP'
  ),
  class = c(rep('marginal', 12), rep('joint', 3), rep('two-step', 2))
)
df_caller <- df_caller %>%
  inner_join( callers, by = 'name_caller' )

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
p_perf <- plot_perf( df_perf )
ggsave( file.path( plot_dir, 'Fig1.SRSV.performance.pdf'), plot = p_perf, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig1.SRSV.performance.png'), plot = p_perf, width = 8, height = 10)

# performance by admixture regime
# ------------------------------------------------------------------------------
df_perf <- readRDS( file.path(data_dir, 'df_perf.rds') )
df <- df_perf %>% mutate( ttype = fct_recode(ttype, 'high'='us', 'med'='ms', 'low'='hs') )
p_perf_tt <- plot_perf_admix( df )
ggsave( file.path( plot_dir, 'Fig2.SRSV.performance.admix.pdf'), plot = p_perf_tt, width = 8, height = 10)
ggsave( file.path( plot_dir, 'Fig2.SRSV.performance.admix.png'), plot = p_perf_tt, width = 8, height = 10)


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

################################################################################
# ALT read counts for TP, FP, FN variants
################################################################################

df_vars <- readRDS( file.path(data_dir, 'df_vars.rds') )
df_vars <- df_caller %>% inner_join( df_vars, by = 'id_caller' )

p_rc_alt <- plot_rc_alt_ridges_srsv( df_vars, df_rc, df_rep )
ggsave( file.path( plot_dir, 'Fig2.SRSV.rc_alt.ridges.pdf'), plot = p_rc_alt, device = pdf(), width = 8, height = 10 )
ggsave( file.path( plot_dir, 'Fig2.SRSV.rc_alt.ridges.png'), plot = p_rc_alt, device = png(), width = 8, height = 10 )
# p_rc_alt <- plot_rc_alt_bar_srsv( df_vars, df_rc, df_rep )
# ggsave( file.path( plot_dir, 'SRSV.vaf.bar.pdf'), plot = p_rc_alt, device = pdf(), width = 8, height = 10 )
# ggsave( file.path( plot_dir, 'SRSV.vaf.bar.png'), plot = p_rc_alt, device = png(), width = 8, height = 10 )


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


