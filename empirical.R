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
require(RColorBrewer)
require(UpSetR)
require(gridExtra)
#require(cowplot)
require(ggpubr)

data_dir = 'data/empirical'
plot_dir = 'plot/empirical'
source( file.path('analysis', 'similarity.analysis.R') )
source( file.path('plotting', 'similarity.plotting.R') )
source( file.path('plotting', 'performance.plotting.R') )
source( file.path('plotting', 'rc_alt.plotting.R') )
source( file.path('plotting', 'upset.plotting.R') )

df_caller <- tibble(
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
    'SNV-PPILP',
    'HaplotypeCaller', 
    'MultiSNV', 
    'Mutect2_multi'
  ),
  class = c(rep('marginal', 11), rep('two-step', 2), rep('joint', 3))
)

# LOAD DATA
# ------------------------------------------------------------------------------

df_sample <- read.csv( file.path(data_dir, 'df_samples.csv') )
df_vars <- read.csv( file.path(data_dir, 'df_vars.csv') )
df_pres <- read.csv( file.path(data_dir, 'df_pres.csv') )
# rename SNV-PPILP column (necessary if loaded from file)
names(df_pres)[names(df_pres)=='SNV.PPILP'] <- 'SNV-PPILP'

################################################################################
# cross-reference with published SNVs
################################################################################

df_snv_pr1 <- read_csv( file.path(data_dir, 'Liu2017.PR1.snvs.sid.csv'),
                        col_types = cols(.default = col_character(),
                                         pos = col_integer()) )
df_snv_pr2 <- read_csv( file.path(data_dir, 'Liu2017.PR2.snvs.sid.csv'),
                        col_types = cols(.default = col_character(),
                                         pos = col_integer()) )
df_snv_mss1 <- read_csv( file.path(data_dir, 'Lim2016.MSS1.snvs.sid.csv'),
                         col_types = cols(.default = col_character(),
                                          pos = col_integer()) )
df_snv_msih1 <- read_csv( file.path(data_dir, 'Lim2016.MSI-H1.snvs.sid.csv'),
                          col_types = cols(.default = col_character(),
                                           pos = col_integer()) )
df_snv_pub <- df_snv_pr1 %>% bind_rows( df_snv_pr2 ) %>%
  bind_rows( df_snv_mss1 ) %>% bind_rows( df_snv_msih1 ) %>%
  mutate( id_sample = sample ) %>%
  unite( id_mut, sample, chrom, pos )

df_snv_calls <- df_vars %>% as_tibble() %>%
  mutate( sample = id_sample ) %>% 
  unite( id_mut, sample, chrom, pos )

df_tp <- df_snv_calls %>% select(name_caller, id_sample, id_mut) %>% 
  inner_join( df_snv_pub %>% select(id_mut) ) %>%
  mutate( type = 'TP' )
df_fp <- df_snv_calls %>% select(name_caller, id_sample, id_mut) %>% 
  anti_join( df_snv_pub %>% select(id_mut) ) %>%
  mutate( type = 'FP' )
df_fn <- df_snv_calls %>% select(name_caller) %>% 
  crossing( df_snv_pub %>% select( id_sample, id_mut ) ) %>%
  anti_join( df_snv_calls ) %>%
  mutate( type = 'FN' )
df_vars_typed <- df_tp %>%
  bind_rows( df_fp ) %>%
  bind_rows( df_fn )

# calculate performance metrics
# ------------------------------------------------------------------------------
df_eval <- df_vars_typed %>%
  group_by( name_caller, id_sample, type ) %>%
  summarise( n = n() ) %>%
  ungroup() %>%
  complete( name_caller, id_sample, type, fill = list(n = 0) ) %>%
  spread( type, n ) %>%
  mutate( recall = TP/(TP+FN), precision = TP/(TP+FP) ) %>%
  mutate( F1 = 2*recall*precision/(recall+precision) )
# add caller information
df_eval <- df_eval %>%
  inner_join( df_caller, by = c('name_caller') ) %>%
  inner_join( df_sample, by = c('id_sample') ) %>%
  replace_na( list(precision = 0, F1 = 0) ) %>%
  mutate( caller = name_caller )

p_perf_box <- plot_perf_empirical( df_eval )
ggsave( file.path( plot_dir, 'empirical.performance.pdf'), plot = p_perf_box, width = 8, height = 10)
ggsave( file.path( plot_dir, 'empirical.performance.png'), plot = p_perf_box, width = 8, height = 10)

p_perf_pts <- plot_perf_empirical_by_patient( df_eval )
ggsave( file.path( plot_dir, 'empirical.performance.by_patient.pdf'), plot = p_perf_pts, width = 8, height = 10)
ggsave( file.path( plot_dir, 'empirical.performance.by_patient.png'), plot = p_perf_pts, width = 8, height = 10)

df <- df_pres %>% inner_join( df_snv_pub %>% select(id_mut) )

# Barplots of recalled variants by caller and sample
p1 <- df %>% dplyr::filter( grepl('^MSS1', id_sample) ) %>% 
  gather( "caller", "present", -id_sample, -id_mut ) %>%
  group_by( id_sample, caller ) %>% summarise( n = sum(present) ) %>%
  ggplot( aes(x = caller, y = n) ) + 
  geom_bar( stat = 'identity' ) +
  geom_hline( data = df %>% dplyr::filter(grepl('^MSS1', id_sample)) %>% group_by(id_sample) %>% tally(), aes(yintercept=n), lty=2 ) +
  facet_wrap( ~id_sample, ncol = 1 ) +
  ggtitle( 'MSS1' ) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1) )
p2 <- df %>% dplyr::filter( grepl('^MSI-H1', id_sample) ) %>% 
  gather( "caller", "present", -id_sample, -id_mut ) %>%
  group_by( id_sample, caller ) %>% summarise( n = sum(present) ) %>%
  ggplot( aes(x = caller, y = n) ) + 
  geom_bar( stat = 'identity' ) +
  geom_hline( data = df %>% dplyr::filter(grepl('^MSI-H1', id_sample)) %>% group_by(id_sample) %>% tally(), aes(yintercept=n), lty=2 ) +
  facet_wrap( ~id_sample, ncol = 1 ) +
  ggtitle( 'MSI-H1' ) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1) )
p3 <- df %>% dplyr::filter( grepl('^PR1', id_sample) ) %>% 
  gather( "caller", "present", -id_sample, -id_mut ) %>%
  group_by( id_sample, caller ) %>% summarise( n = sum(present) ) %>%
  ggplot( aes(x = caller, y = n) ) + 
  geom_bar( stat = 'identity' ) +
  geom_hline( data = df %>% dplyr::filter(grepl('^PR1', id_sample)) %>% group_by(id_sample) %>% tally(), aes(yintercept=n), lty=2 ) +
  facet_wrap( ~id_sample, ncol = 1 ) +
  ggtitle( 'PR1' ) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1) )
p4 <- df %>% dplyr::filter( grepl('^PR2', id_sample) ) %>% 
  gather( "caller", "present", -id_sample, -id_mut ) %>%
  group_by( id_sample, caller ) %>% summarise( n = sum(present) ) %>%
  ggplot( aes(x = caller, y = n) ) + 
  geom_bar( stat = 'identity' ) +
  geom_hline( data = df %>% dplyr::filter(grepl('^PR2', id_sample)) %>% group_by(id_sample) %>% tally(), aes(yintercept=n), lty=2 ) +
  facet_wrap( ~id_sample, ncol = 1 ) +
  ggtitle( 'PR2' ) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1) )
pdf( file.path(plot_dir, 'empirical.recall.published.pdf'), 12, 12, onefile = TRUE )
#grid.arrange( p1, p2, p3, p4, nrow = 1 )
ggarrange( p1, p2, p3, p4, nrow = 1, ncol = 4,
           labels = c("A)", "B)", "C)", "D)") )
dev.off()

# correlation between F1 score and recall, precision
# ------------------------------------------------------------------------------
cor.test( df_eval$F1, df_eval$recall )
#       Pearson's product-moment correlation
# 
# data:  df_eval$F1 and df_eval$recall
# t = 6.1244, df = 238, p-value = 3.733e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2542223 0.4734815
# sample estimates:
#   cor 
# 0.3689741
cor.test( df_eval$F1, df_eval$precision )
#       Pearson's product-moment correlation
# 
# data:  df_eval$F1 and df_eval$precision
# t = 46.333, df = 238, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9344251 0.9600707
# sample estimates:
#   cor 
# 0.9487885 
# ------------------------------------------------------------------------------

# UpSet plot
n <- c('Bcftools', 'CaVEMan', 'HaplotypeCaller', 'MuClone', 'MultiSNV', 'MuTect1', 
       'Mutect2_multi', 'Mutect2_single', 'NeuSomatic', 'Shimmer', 'SNV-PPILP', 
       'SNooPer', 'SomaticSniper', 'Strelka2', 'VarDict', 'VarScan')
pdf( file.path(plot_dir, 'empirical.upset.pub.pdf'), width = 16, height = 8, onefile = FALSE )
#png( 'plots/empirical.upset.pub.png', width = 10, height = 4.5 ) # does produce empty output...
upset(df, sets = n, 
      keep.order = FALSE, 
      order.by = 'freq',
      mb.ratio = c(0.5, 0.5),
      number.angles = 30,
      text.scale = 1.5)
dev.off()


################################################################################
# Similarity between call sets
################################################################################

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
## Jaccard distance plot
df_jacc <- Jaccard.df( df_pres %>% select(-id_sample, -id_mut) %>% select(callerorder))
p_jacc <- plot_jacc_idx( df_jacc %>% mutate(caller1 = factor(caller1, levels = callerorder), 
                                                  caller2 = factor(caller2, levels = callerorder)))

ggsave( file.path( plot_dir, 'empirical.jaccard.pdf'), plot = p_jacc, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'empirical.jaccard.png'), plot = p_jacc, device = png(), width = 8, height = 8 )


# hierarchical clustering based on varcalls
#-------------------------------------------------------------------------------
require(ade4) # dist.binary()
require(ggdendro) # ggdendrogram()
#df_pres <- df_pres_tp %>% bind_rows( df_pres_fn ) %>% bind_rows( df_pres_fp )
df_jacc <- Jaccard.df( df_pres %>% select(-id_sample, -id_mut) )

df_jacc_idx <- df_jacc %>% spread( caller1, jaccard_idx ) %>% as.data.frame() 
df_jacc_idx <- df_jacc_idx %>% set_rownames( df_jacc_idx$caller2 ) %>% select( -caller2 )
d <- as.dist( 1-df_jacc_idx )
hc <- hclust(d)

fn_pfx <- 'empirical.hclust.dendro'
pdf( file.path(plot_dir, paste0(fn_pfx, '.pdf')), width = 8, height = 4 )
par( mar = c(2, 1, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()
png( file.path(plot_dir, paste0(fn_pfx, '.png')), width = 8, height = 4, units = 'in', res = 300 )
par( mar = c(2, 0, 0, 6) )
plot(as.dendrogram(hc), horiz = TRUE)
dev.off()

# UPSET PLOT
# ------------------------------------------------------------------------------
df_pres <- read.csv( file.path(data_dir, 'df_pres.csv') )
# annoying, but necessary... only if loaded from file
names(df_pres)[names(df_pres)=='SNV.PPILP'] <- 'SNV-PPILP'
n <- c('Bcftools', 'CaVEMan', 'HaplotypeCaller', 'MuClone', 'MultiSNV', 'MuTect1', 
       'Mutect2_multi', 'Mutect2_single', 'NeuSomatic', 'Shimmer', 'SNV-PPILP', 
       'SNooPer', 'SomaticSniper', 'Strelka2', 'VarDict', 'VarScan')
#png( 'plots/empirical.upset.all.png', width = 10, height = 4.5 ) # does produce empty output...
fn_pfx <- file.path( plot_dir, 'empirical.upset')
pdf( paste0(fn_pfx, '.pdf'), width = 8, height = 6, onefile = FALSE )
plot_upset_empirical( df_pres, n )
dev.off()
png( paste0(fn_pfx, '.png'), width = 8, height = 6, units = 'in', res = 300 )
plot_upset_empirical( df_pres, n )
dev.off()

# ratio of private variants for each caller
# ------------------------------------------------------------------------------
df_uniq <- df_pres %>% mutate( n = rowSums(select(., -c(1:2))) ) %>% dplyr::filter( n == 1 ) %>%
  gather( "caller", "present", -id_sample, -id_mut, -n ) %>% group_by( caller ) %>% summarise( n_unique = sum(present) )
df_total <- df_pres %>% 
  gather( "caller", "present", -id_sample, -id_mut ) %>% group_by( caller ) %>% summarise( n_total = sum(present) )
df_uniq <- df_uniq %>% inner_join( df_total ) %>% mutate( r_unique = n_unique/n_total )
p_uniq_calls <- ggplot( df_uniq ) + 
  geom_bar( aes(x = reorder(caller, -n_total), y = r_unique), stat = 'identity' ) +
  geom_text( aes(x = reorder(caller, -n_total), y = 0.75, label = n_total, color = (r_unique > 0.75)) ) +
  labs( x = 'caller', y = 'unique/total calls' ) +
  guides( color = 'none' ) +
  coord_flip() +
  theme_minimal() +
  scale_color_manual( values = c('black', 'white') )
ggsave( plot = p_uniq_calls, filename = file.path(plot_dir, 'empirical.private_calls.pdf' ), device = pdf() )
ggsave( plot = p_uniq_calls, filename = file.path(plot_dir, 'empirical.private_calls.png' ), device = png() )

# by sample
df_uniq <- df_pres %>% mutate( n = rowSums(select(., -c(1:2))) ) %>% dplyr::filter( n == 1 ) %>%
  gather( "caller", "present", -id_sample, -id_mut, -n ) %>% group_by( id_sample, caller ) %>% summarise( n_unique = sum(present) )
df_total <- df_pres %>% 
  gather( "caller", "present", -id_sample, -id_mut ) %>% group_by( id_sample, caller ) %>% summarise( n_total = sum(present) )
df_uniq <-df_uniq %>% inner_join( df_total ) %>% mutate( r_unique = n_unique/n_total )
ggplot( df_uniq ) + geom_bar( aes(x = reorder(caller, -r_unique), y = r_unique), stat = 'identity' ) +
  labs( x = 'caller', y = 'unique/total calls' ) +
  coord_flip() +
  facet_wrap( ~id_sample, ncol = 3 )

# repeat uniqueness tally pooling Mutect2 and Mutect1 calls
df_pres_mod <- df_pres %>% mutate( Mutect_pooled = ifelse(MuTect1 + Mutect2_single + Mutect2_multi == 0, 0, 1) ) %>%
  select( -MuTect1, -Mutect2_single, -Mutect2_multi )
df_uniq <- df_pres_mod %>% mutate( n = rowSums(select(., -c(1:2))) ) %>% dplyr::filter( n == 1 ) %>%
  gather( "caller", "present", -id_sample, -id_mut, -n ) %>% group_by( caller ) %>% summarise( n_unique = sum(present) )
df_total <- df_pres_mod %>% 
  gather( "caller", "present", -id_sample, -id_mut ) %>% group_by( caller ) %>% summarise( n_total = sum(present) )
df_uniq <- df_uniq %>% inner_join( df_total ) %>% mutate( r_unique = n_unique/n_total )
p_uniq_calls <- ggplot( df_uniq ) + geom_bar( aes(x = reorder(caller, -n_total), y = r_unique), stat = 'identity' ) +
  labs( x = 'caller', y = 'unique/total calls' ) +
  coord_flip() +
  theme_minimal()
ggsave( plot = p_uniq_calls, filename = file.path(plot_dir, 'empirical.private_calls.mutect_pooled.pdf' ), device = pdf() )
ggsave( plot = p_uniq_calls, filename = file.path(plot_dir, 'empirical.private_calls.mutect_pooled.png' ), device = png() )



## inspect vars manually
df <- df_pres %>% mutate(n=rowSums(.[3:16])) %>% dplyr::filter(n==13 & MultiSNV==0)

# load data from files
df_sample <- read_csv( file.path(data_dir, 'df_samples.csv') )
df_vars <- read_csv( file.path(data_dir, 'df_vars.csv'),
                     col_types = cols(
                       .default = col_character(),
                       pos = col_integer()) )

# number of variants per sample
ggplot( df_vars, aes(id_sample) ) + geom_histogram(stat="count") +
  coord_flip() +
  facet_wrap( ~name_caller, scales = "free_x" )


################################################################################
# Variant allele frequencies
################################################################################

# HELPER FUNCTIONS
# ------------------------------------------------------------------------------
parse_readcounts <- function( fn ) {
  read_delim( fn, delim = '\t', 
              comment = '@',
              col_types = cols(CONTIG = 'c' ) ) %>% 
    mutate( filename = basename(fn) )
}
# ------------------------------------------------------------------------------

# VAF spectrum for each caller
fn_rc <- file.path(
  data_dir, 
  sprintf('%s.vars.alleliccounts.tsv', 
          df_sample %>% dplyr::filter( source == 'tumour' ) %>% 
            pull(name_sample)) )
df_rc <- fn_rc %>% map_df( ~parse_readcounts(.) ) %>%
  mutate( name_sample = str_extract(filename, '^([^\\.])+') ) %>%
  inner_join( df_sample ) %>%
  #unite( id_mut, c('id_sample', 'CONTIG', 'POSITION'), remove = FALSE ) %>%
  mutate( vaf = ALT_COUNT/(REF_COUNT+ALT_COUNT),
          depth = REF_COUNT+ALT_COUNT) %>%
  select( id_sample, chrom = CONTIG, pos = POSITION, rc_ref = REF_COUNT, rc_alt = ALT_COUNT, depth, vaf  )

p_vaf <- plot_vaf_bar_empirical( df_vars, df_rc )
ggsave( file.path( plot_dir, 'empirical.vaf.bar.pdf'), plot = p_vaf, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'empirical.vaf.bar.png'), plot = p_vaf, device = png(), width = 8, height = 8 )

p_vaf_sample <- plot_vaf_bar_empirical_sample( df_vars, df_rc, df_sample )
ggsave( file.path( plot_dir, 'empirical.vaf.bar.sample.pdf'), plot = p_vaf_sample, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'empirical.vaf.bar.sample.png'), plot = p_vaf_sample, device = png(), width = 8, height = 8 )

# VAF spectrum for private vs. shared calls
p_vaf_priv <- plot_vaf_bar_empirical_private( df_vars, df_rc, df_pres )
ggsave( file.path( plot_dir, 'empirical.vaf.bar.private.pdf'), plot = p_vaf_priv, device = pdf(), width = 8, height = 8 )
ggsave( file.path( plot_dir, 'empirical.vaf.bar.private.png'), plot = p_vaf_priv, device = png(), width = 8, height = 8 )

