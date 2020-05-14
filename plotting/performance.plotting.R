plot_perf_min <- function ( df )
{
  # exclude these callers from plots
  blacklist <- c( 'Strelka1' )
  
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
    'Bcftools', 
    'CaVEMan', 
    'MuTect1', 
    'Mutect2_single', 
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
    'Mutect2_multi_F'))
  df$class <- factor( df$class, levels = c('marginal', 'two-step', 'joint') )
  
  # format caller names for better plotting
  df <- df %>% mutate( lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T) )
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi_F',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # remove data for some callers (older versions of certain methods)
  df <- df %>% dplyr::filter( !(caller %in% blacklist) )
  
  # subplots for grid layout
  p_rec <- ggplot( df, aes(x = caller, y = recall) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(caller) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = caller, y = mrec), fill = "gold", shape = 23) + 
    labs( x = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )
  
  p_pre <- ggplot( df, aes(x = caller, y = precision) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(caller) %>% summarise(mpre = median(precision)) %>% arrange(desc(mpre)) %>% dplyr::filter(mpre==max(mpre)), aes(x = caller, y = mpre), fill = "gold", shape = 23) + 
    labs( x = 'caller' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )
  
  p_f1 <- ggplot( df, aes(x = caller, y = F1) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(caller) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x = caller, y = mF1), fill = "gold", shape = 23) + 
    labs( x = '', fill = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    theme(legend.position = 'right')
    
    # theme( legend.position = 'bottom' ) +
    # guides( colour = 'none' )+
    # scale_fill_discrete(name="")
  
  p_perf <- ggarrange(
    p_rec, p_pre, p_f1,
    ncol = 3,
    nrow = 1,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "top"
  )
  
  return( p_perf )
}

plot_perf_cvg <- function ( df )
{
  # exclude these callers from plots
  blacklist <- c( 'Strelka1', 'MuClone_perf' )
  
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
    'Bcftools', 
    'CaVEMan', 
    'MuTect1', 
    'Mutect2_single', 
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
    'Mutect2_multi_F'))
  df$class <- factor(df$class, levels = c('marginal', 'two-step', 'joint'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi_F',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # remove data for some callers (older versions of certain methods)
  df <- df %>% dplyr::filter( !(caller %in% blacklist) )
  
  # subplots for grid layout
  p_r_cvg <- ggplot(df, aes(x = as.factor(cvg), y = recall)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = factor(cvg), y=mrec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(.~lbl) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = 'none', alpha = "none")
  
  p_p_cvg <- ggplot(df, aes(x = as.factor(cvg), y = precision)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(cvg), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth')  +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_cvg <- ggplot(df, aes(x = as.factor(cvg), y = F1)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth', y = 'F1 score', fill = '') +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    theme(strip.text.x = element_text(size = 6)) +
    theme(legend.position = 'bottom') +
    guides(colour = "none", alpha = "none")
  
  p_perf <- ggarrange(
    p_r_cvg, p_p_cvg, p_f_cvg,
    ncol = 1,
    nrow = 3,
    labels = "auto", 
    common.legend = TRUE, 
    legend = "bottom"
  )
  
  return( p_perf )
}

# plot performance at different depth of coverage
# test significance of difference between depth levels for each caller
# params:
#  df     - caller, recall, precision, F1, cvg
#  method - one of 'anova', 'kruskal'
plot_perf_cvg_sig <- function ( df, method )
{
  # exclude these callers from plots
  blacklist <- c( 'Strelka1', 'MuClone_perf' )
  
  # remove data for some callers (older versions of certain methods)
  df <- df %>% dplyr::filter( !(caller %in% blacklist) )
  
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
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
    'Mutect2_multi_F'))
  df$class <- factor(df$class, levels = c('marginal', 'two-step', 'joint'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi_F',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # significance levels
  sig_cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1)
  sig_symbols = c("****","***", "**", "*", "ns")
  sig_legend = sprintf('Group-wise comparison of means by factor; method: "%s"; p.adjust: "bonferroni"\n(%s)', method, paste(sprintf('%s: p<%s', sig_symbols, sig_cutpoints[-1]), collapse = ', '))
  sig_legend = sprintf('Group-wise comparison of means by factor\nmethod: "A: ANOVA, K: Kruskal-Wallis"; p.adjust: "bonferroni"\n(%s)', paste(sprintf('%s: p<%s', sig_symbols, sig_cutpoints[-1]), collapse = ', '))
  
  # subplots for grid layout
  p_r_cvg <- ggboxplot(df, x = 'cvg', y = 'recall', fill = 'class', yticks.by = 0.25) +
    #geom_text( data = gwc_rec, aes(label = sig_sym, y = 1.05, x = 2) ) +
    facet_wrap( ~lbl, nrow = 1) +
    rotate_x_text()
    #stat_compare_means(label = '..p.signif..', method = method, label.y = 1.05)
  
  p_p_cvg <- ggboxplot(df, x = 'cvg', y = 'precision', fill = 'class', yticks.by = 0.25) %>%
    facet(facet.by = 'lbl', nrow = 1) +
    rotate_x_text()
    #geom_text( data = gwc_rec, aes(label = sig_sym, y = 1.05, x = 2) )
    #stat_compare_means(label = '..p.signif..', method = method, label.y = 1.05)

    # group-wise comparison of differences in means
  # gwc_f1 <- df %>% group_by( lbl ) %>%
  # { if (method == 'anova' ) anova_test(., F1 ~ cvg) else . } %>%
  # { if (method == 'kruskal' ) kruskal_test(., F1 ~ cvg) else . } %>%
  #   adjust_pvalue( method = 'bonferroni' ) %>%
  #   add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  gwc_f1_anova <- df %>% group_by( lbl ) %>%
    anova_test( F1 ~ cvg) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  gwc_f1_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., F1 ~ cvg) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  p_f_cvg <- ggboxplot(df, x = 'cvg', y = 'F1', fill = 'class', yticks.by = 0.25) %>%
    facet(facet.by = 'lbl', nrow = 1) +
    labs(x = 'sequencing depth')  +
    rotate_x_text() +
    geom_text( data = gwc_f1_kruskal, aes(label = sprintf('K: %s', sig_sym), y = 1.05, x = 2.5) ) +
    geom_text( data = gwc_f1_anova, aes(label = sprintf('A: %s', sig_sym), y = 1.15, x = 2.5) )
    #geom_text( data = gwc_f1, aes(label = sig_sym, y = 1.05, x = 2.5) )
    #stat_compare_means(label = '..p.signif..', method = method, label.y = 1.05) # Add pairwise comparisons p-value
  #font("xy.text", size = 12)
  #ggsave(plot = p_f_cvg, filename = 'de-novo.f1.kruskal-wallis.pdf', width = 12, height = 9)
  #ggsave(plot = p_f_cvg, filename = 'de-novo.f1.kruskal-wallis.png', width = 12, height = 9)
  
  
  # combine subplots
  p_perf <- ggarrange(
    p_r_cvg, p_p_cvg, p_f_cvg,
    ncol = 1,
    nrow = 3,
    labels = "auto", 
    common.legend = TRUE, 
    legend = "bottom"
  ) %>%
  annotate_figure(
    bottom = text_grob(sig_legend, hjust = 1, x = 1, face = "italic", size = 10))
  
  return( p_perf )
}

# plot performance at different depth of coverage
# test significance of difference between depth levels for each caller
# params:
#  df     - caller, recall, precision, F1, cvg
#  method - one of 'anova', 'kruskal'
plot_perf_admix_sig <- function ( df, method )
{
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
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
    'Mutect2_multi_F'))
  df$class <- factor(df$class, levels = c('marginal', 'two-step', 'joint'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi_F',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # significance levels
  sig_cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1)
  sig_symbols = c("****","***", "**", "*", "ns")
  sig_legend = sprintf('Group-wise comparison of means by factor; method: "%s"; p.adjust: "bonferroni"\n(%s)', method, paste(sprintf('%s: p<%s', sig_symbols, sig_cutpoints[-1]), collapse = ', '))
  sig_legend = sprintf('Group-wise comparison of means by factor\nmethod: "A: ANOVA, K: Kruskal-Wallis"; p.adjust: "bonferroni"\n(%s)', paste(sprintf('%s: p<%s', sig_symbols, sig_cutpoints[-1]), collapse = ', '))
  
  # subplots for grid layout
  p_rec <- ggboxplot(df, x = 'ttype', y = 'recall', fill = 'class', yticks.by = 0.25) +
    #geom_text( data = gwc_rec, aes(label = sig_sym, y = 1.05, x = 2) ) +
    facet_wrap( ~lbl, nrow = 1) +
    rotate_x_text()
  #stat_compare_means(label = '..p.signif..', method = method, label.y = 1.05)
  
  p_pre <- ggboxplot(df, x = 'ttype', y = 'precision', fill = 'class', yticks.by = 0.25) %>%
    facet(facet.by = 'lbl', nrow = 1) +
    rotate_x_text()
  #geom_text( data = gwc_rec, aes(label = sig_sym, y = 1.05, x = 2) )
  #stat_compare_means(label = '..p.signif..', method = method, label.y = 1.05)
  
  # group-wise comparison of differences in means
  # gwc_f1 <- df %>% group_by( lbl ) %>%
  # { if (method == 'anova' ) anova_test(., F1 ~ cvg) else . } %>%
  # { if (method == 'kruskal' ) kruskal_test(., F1 ~ cvg) else . } %>%
  #   adjust_pvalue( method = 'bonferroni' ) %>%
  #   add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  gwc_f1_anova <- df %>% group_by( lbl ) %>%
    anova_test( F1 ~ ttype) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  gwc_f1_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., F1 ~ ttype) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  p_f1 <- ggboxplot(df, x = 'ttype', y = 'F1', fill = 'class', yticks.by = 0.25) %>%
    facet(facet.by = 'lbl', nrow = 1) +
    labs(x = 'admixture')  +
    rotate_x_text() +
    geom_text( data = gwc_f1_kruskal, aes(label = sprintf('K: %s', sig_sym), y = 1.05, x = 2.5) ) +
    geom_text( data = gwc_f1_anova, aes(label = sprintf('A: %s', sig_sym), y = 1.15, x = 2.5) )
  #geom_text( data = gwc_f1, aes(label = sig_sym, y = 1.05, x = 2.5) )
  #stat_compare_means(label = '..p.signif..', method = method, label.y = 1.05) # Add pairwise comparisons p-value
  #font("xy.text", size = 12)
  #ggsave(plot = p_f_cvg, filename = 'de-novo.f1.kruskal-wallis.pdf', width = 12, height = 9)
  #ggsave(plot = p_f_cvg, filename = 'de-novo.f1.kruskal-wallis.png', width = 12, height = 9)
  
  
  # combine subplots
  p_perf <- ggarrange(
    p_rec, p_pre, p_f1,
    ncol = 1,
    nrow = 3,
    labels = "auto", 
    common.legend = TRUE, 
    legend = "bottom"
  ) %>%
    annotate_figure(
      bottom = text_grob(sig_legend, hjust = 1, x = 1, face = "italic", size = 10))
  
  return( p_perf )
}

plot_perf_rrsv <- function ( df )
{
  # exclude these callers from plots
  blacklist <- c( 'Strelka1' )
  
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
    'Bcftools', 
    'CaVEMan', 
    'MuTect1', 
    'Mutect2_single', 
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
    'Mutect2_multi'))
  df$class <- factor(df$class, levels = c('marginal', 'two-step', 'joint'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # remove data for some callers (older versions of certain methods)
  df <- df %>% dplyr::filter( !(caller %in% blacklist) )
  
  # subplots for grid layout
  p_r_cvg <- ggplot( df, aes(x = as.factor(cvg), y = recall) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = as.factor(cvg), y=mrec), fill = "gold", shape = 23) + 
    labs( x = 'caller', fill = '' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank()) + 
    guides( colour = "none", fill = 'none' )
  
  p_p_cvg <- ggplot( df, aes(x = as.factor(cvg), y = precision) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mpre = median(precision)) %>% arrange(desc(mpre)) %>% dplyr::filter(mpre==max(mpre)), aes(x = as.factor(cvg), y=mpre), fill = "gold", shape = 23) + 
    labs( x = 'caller' ) + 
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) + 
    guides( colour = "none", fill = 'none' )
  
  p_f_cvg <- ggplot( df, aes(x = as.factor(cvg), y = F1) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=as.factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    labs( x = 'caller', y = 'F1 score' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) +
    #theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
    theme( legend.position = 'bottom' ) +
    guides( colour = 'none' )+
    scale_fill_discrete(name="")
  
  p_perf <- ggarrange(
    p_r_cvg, p_p_cvg, p_f_cvg,
    ncol = 1,
    nrow = 3,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "bottom"
 

  )
  
  return( p_perf )
}

plot_perf_pairs <- function ( df )
{
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-|_)([A-Z])", "\n\\1", df$caller, perl = T))
  # df <- df %>% 
  #   mutate(lbl = fct_recode(caller, 
  #                           'Haplotype\nCaller' = 'HaplotypeCaller',
  #                           'Mutect2\nmulti_F' = 'Mutect2_multi_F',
  #                           'Mutect2\nsingle' = 'Mutect2_single',
  #                           'Neu\nSomatic' = 'NeuSomatic',
  #                           'Somatic\nSniper' = 'SomaticSniper',
  #                           'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # reorder callers by median F1 score
  df$lbl <- factor( df$lbl, levels = df %>% group_by(lbl) %>% summarise(med_F1 = median(F1)) %>% arrange(desc(med_F1)) %>% pull(lbl) )
  
  # subplots for grid layout
  p_r_cvg <- ggplot( df, aes(x = as.factor(cvg), y = recall) ) + 
    theme_minimal()+
    geom_boxplot() + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = as.factor(cvg), y=mrec), fill = "gold", shape = 23) + 
    labs( x = 'caller', fill = '' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank()) + 
    guides( colour = "none", fill = 'none' )
  
  p_p_cvg <- ggplot( df, aes(x = as.factor(cvg), y = precision) ) + 
    theme_minimal()+
    geom_boxplot() + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mpre = median(precision)) %>% arrange(desc(mpre)) %>% dplyr::filter(mpre==max(mpre)), aes(x = as.factor(cvg), y=mpre), fill = "gold", shape = 23) + 
    labs( x = 'caller' ) + 
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) + 
    guides( colour = "none", fill = 'none' )
  
  p_f_cvg <- ggplot( df, aes(x = as.factor(cvg), y = F1) ) + 
    theme_minimal()+
    geom_boxplot() + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=as.factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    labs( x = 'caller', y = 'F1 score' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) +
    #theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
    theme( legend.position = 'bottom' ) +
    guides( colour = 'none' )+
    scale_fill_discrete(name="")
  
  p_perf <- ggarrange(
    p_r_cvg, p_p_cvg, p_f_cvg,
    ncol = 1,
    nrow = 3,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "bottom"
    
    
  )
  
  return( p_perf )
}

plot_perf_pairs_hm <- function ( df )
{
  # format caller names for better plotting
  df <- df %>% mutate(
    lbl1 = gsub("(?<=\\w{5}|-|_)([A-Z])(?=\\w{3})", "\n\\1", df$caller1, perl = T),
    lbl2 = gsub("(?<=\\w{5}|-|_)([A-Z])(?=\\w{3})", "\n\\1", df$caller2, perl = T)
  )
  
  # define plotting  order of callers (by id)
  df$lbl1 <- factor( df$lbl1, levels = df %>% select(id_caller1, lbl1) %>% unique() %>% arrange((id_caller1)) %>% pull(lbl1) )
  df$lbl2 <- factor( df$lbl2, levels = df %>% select(id_caller2, lbl2) %>% unique() %>% arrange(desc(id_caller2)) %>% pull(lbl2) )
  
  p_mF1_hm <- ggplot( df, aes(x = lbl1, y = lbl2, fill = med_F1) ) + 
    geom_tile() +
    geom_text( aes(label = round(med_F1, 2)) ) +
    scale_fill_distiller( palette = 'Spectral' ) +
    labs( x = 'caller1', y = 'caller2', fill = 'median\nF1' ) +
    theme_minimal() +
    theme( axis.text.x = element_text(angle = 20) )
  
  return( p_mF1_hm )
}

plot_perf_trip_hm <- function ( df )
{
  library(ggpubr)  # as_ggplot()
  
  # format caller names for better plotting
  regex_str <- '(?<=\\w{5}|-|_)([A-Z])(?=\\w{3})'
  df <- df %>% mutate(
    #lbl1 = gsub( regex_str, '\n\\1', df$caller1, perl = T),
    lbl1 = df$caller1,
    lbl2 = gsub( regex_str, '\n\\1', df$caller2, perl = T),
    lbl3 = gsub( regex_str, '\n\\1', df$caller3, perl = T)
  )
  
  # define plotting  order of callers (by id)
  df$lbl1 <- factor( df$lbl1, levels = df %>% select(id_caller1, lbl1) %>% unique() %>% arrange((id_caller1)) %>% pull(lbl1) )
  df$lbl2 <- factor( df$lbl2, levels = df %>% select(id_caller2, lbl2) %>% unique() %>% arrange(desc(id_caller2)) %>% pull(lbl2) )
  df$lbl3 <- factor( df$lbl3, levels = df %>% select(id_caller3, lbl3) %>% unique() %>% arrange(desc(id_caller3)) %>% pull(lbl3) )
  
  p <- ggplot( df, aes(x = lbl1, y = med_F1, fill = med_F1) ) + 
    geom_col() +
    ylim( c(0, 1) ) +
    #geom_text( aes(label = round(med_F1, 2)) ) +
    #scale_fill_distiller( palette = 'Spectral' ) +
    labs( x = 'caller1', fill = 'median\nF1' ) +
    facet_grid( lbl3 ~ lbl2 ) +
    theme_minimal() +
    theme( axis.text.x = element_text(angle = 90) )
  
  # g <- ggplotGrob( p )
  #gtable::gtable_show_layout(g)
  
  # get gtable columns corresponding to the facets (5 & 9, in this case)
  # facet.columns <- g$layout$l[grepl("panel", g$layout$name)] %>% unique()
  # get the number of unique x-axis values per facet (1 & 3, in this case)
  # x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                  # function(l) length(l$range$range))
  # change the relative widths of the facet columns based on
  # how many unique x-axis values are in each facet
  # g$widths[facet.columns] <- g$widths[facet.columns] * x.var
  
  # remove empty panels
  # pos <- g$layout %>% 
  #   rowid_to_column( 'row_name' ) %>%
  #   dplyr::filter( str_detect(name, 'panel') ) %>% 
  #   separate( name, c('type', 'col', 'row') ) %>%
  #   dplyr::filter( row > col ) %>%
  #   pull( row_name )
  # hide <- 1:nrow(g$layout) %in% pos
  # g$grobs <- g$grobs[!hide]
  # g$layout <- g$layout[!hide, ]
  
  # Draw the plot
  # grid.newpage()
  # grid.draw(g)
  
  # convert gtable back into ggplot
  # p_mF1_trio <- as_ggplot( g )
  
  return( p )
}

plot_perf_empirical <- function ( df )
{
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
    'Bcftools', 
    'CaVEMan', 
    'MuTect1', 
    'Mutect2_single', 
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
    'Mutect2_multi'))
  df$class <- factor(df$class, levels = c('marginal', 'two-step', 'joint'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # insert dummy column for coverage
  df$cvg = 100
  
  # subplots for grid layout
  p_r_cvg <- ggplot( df, aes(x = as.factor(cvg), y = recall) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = as.factor(cvg), y=mrec), fill = "gold", shape = 23) + 
    labs( x = 'caller', fill = '' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank()) + 
    guides( colour = "none", fill = 'none' )
  
  p_p_cvg <- ggplot( df, aes(x = as.factor(cvg), y = precision) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mpre = median(precision)) %>% arrange(desc(mpre)) %>% dplyr::filter(mpre==max(mpre)), aes(x = as.factor(cvg), y=mpre), fill = "gold", shape = 23) + 
    labs( x = 'caller' ) + 
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) + 
    guides( colour = "none", fill = 'none' )
  
  p_f_cvg <- ggplot( df, aes(x = as.factor(cvg), y = F1) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=as.factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    labs( x = 'caller', y = 'F1 score' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) +
    #theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
    theme( legend.position = 'bottom' ) +
    guides( colour = 'none' )+
    scale_fill_discrete(name="")
  
  p_perf <- ggarrange(
    p_r_cvg, p_p_cvg, p_f_cvg,
    ncol = 1,
    nrow = 3,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "bottom"
  )
  
  return( p_perf )
}

plot_perf_empirical_by_patient <- function ( df )
{
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
    'Bcftools', 
    'CaVEMan', 
    'MuTect1', 
    'Mutect2_single', 
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
    'Mutect2_multi'))
  df$class <- factor(df$class, levels = c('marginal', 'two-step', 'joint'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # subplots for grid layout
  p_r_cvg <- ggplot( df, aes(x = patient, y = recall) ) + 
    theme_minimal() +
    #geom_violin( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( aes(fill = class), shape = 21, stroke = 0.5, alpha = 0.6 ) + ylim( 0, 1 ) +
    #geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = as.factor(cvg), y=mrec), fill = "gold", shape = 23) + 
    labs( x = 'patient', fill = '' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    guides( colour = "none", fill = 'none' )
  
  p_p_cvg <- ggplot( df, aes(x = patient, y = precision) ) + 
    theme_minimal() +
    #geom_violin( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( aes(fill = class), shape = 21, stroke = 0.5, alpha = 0.6 ) + ylim( 0, 1 ) +
    #geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mpre = median(precision)) %>% arrange(desc(mpre)) %>% dplyr::filter(mpre==max(mpre)), aes(x = as.factor(cvg), y=mpre), fill = "gold", shape = 23) + 
    labs( x = 'patient' ) + 
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    guides( colour = "none", fill = 'none' )
  
  p_f_cvg <- ggplot( df, aes(x = patient, y = F1) ) + 
    theme_minimal() +
    #geom_violin( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( aes(fill = class), shape = 21, stroke = 0.5, alpha = 0.6 ) + ylim( 0, 1 ) +
    #geom_point( data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=as.factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    labs( x = 'patient', y = 'F1 score' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    #theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
    theme( legend.position = 'bottom' ) +
    guides( colour = 'none' )+
    scale_fill_discrete(name="")
  
  p_perf <- ggarrange(
    p_r_cvg, p_p_cvg, p_f_cvg,
    ncol = 1,
    nrow = 3,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "bottom"
  )
  
  return( p_perf )
}

plot_perf_admix <- function ( df )
{
  # exclude these callers from plots
  blacklist <- c( 'Strelka1' )
  
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
    'Bcftools', 
    'CaVEMan', 
    'MuTect1', 
    'Mutect2_single', 
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
    'Mutect2_multi_F'))
  df$class <- factor(df$class, levels = c('marginal', 'two-step', 'joint'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi_F',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # remove data for some callers (older versions of certain methods)
  df <- df %>% dplyr::filter( !(caller %in% blacklist) )
  
  # subplots for grid layout
  p_r_ttype <- ggplot(df, aes(x = as.factor(ttype), y = recall)) + 
    geom_boxplot(aes(alpha = ttype, fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = factor(ttype), y=mrec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = 'none', alpha = "none")
  
  p_p_ttype <- ggplot(df, aes(x = as.factor(ttype), y = precision)) + 
    geom_boxplot(aes(alpha = factor(ttype), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(ttype), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture') + 
    facet_grid(.~lbl) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_ttype <- ggplot(df, aes(x = as.factor(ttype), y = F1)) + 
    geom_boxplot(aes(alpha = factor(ttype), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(ttype), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture', y = 'F1 score', fill = '') + 
    facet_grid(.~lbl) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    theme(strip.text.x = element_text(size = 6)) +
    theme(legend.position = 'bottom') +
    guides(colour = "none", alpha = "none")
  
  p_perf <- ggarrange(
    p_r_ttype, p_p_ttype, p_f_ttype,
    ncol = 1,
    nrow = 3,
    labels = "auto", 
    common.legend = TRUE, 
    legend = "bottom"
  )
  return( p_perf )
}

plot_perf_pairs_admix <- function ( df )
{
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-|_)([A-Z])", "\n\\1", df$caller, perl = T))
  df$lbl <- factor( df$lbl, levels = df %>% select(lbl, rank_mF1) %>% unique() %>% arrange(rank_mF1) %>% pull(lbl) )
  
  # subplots for grid layout
  p_r_ttype <- ggplot(df, aes(x = as.factor(ttype), y = recall)) + 
    geom_boxplot(aes(alpha = ttype)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = factor(ttype), y=mrec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = 'none', alpha = "none")
  
  p_p_ttype <- ggplot(df, aes(x = as.factor(ttype), y = precision)) + 
    geom_boxplot(aes(alpha = factor(ttype))) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(ttype), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture') + 
    facet_grid(.~lbl) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_ttype <- ggplot(df, aes(x = as.factor(ttype), y = F1)) + 
    geom_boxplot(aes(alpha = factor(ttype))) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(ttype), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture', y = 'F1 score', fill = '') + 
    facet_grid(.~lbl) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    theme(strip.text.x = element_text(size = 6)) +
    theme(legend.position = 'bottom') +
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_perf <- ggarrange(
    p_r_ttype, p_p_ttype, p_f_ttype,
    ncol = 1,
    nrow = 3,
    labels = "auto", 
    common.legend = TRUE, 
    legend = "bottom"
  )
  return( p_perf )
}

plot_pairwise_wilcoxon <- function( df ) {
  # plot Box- and Violinplots of F1 scores for callers
  p_perf_f1_box <- ggviolin( df, x = 'caller', y = 'F1', color = 'class',
                             main = 'a) Distribution of F1 scores',
                             add = 'boxplot' ) +
    rotate_x_text( 45 )
  
  # perform Wilcoxon rank sum test (Mann-Whitney test)
  pwt <- paired.wilcox.test( df$F1, df$caller, p.adjust.method = "BH" )
  df_pwt <- pwt$p.value %>% 
    as_tibble( rownames = 'id1' ) %>% 
    gather( id2, p.adj, -id1) %>% dplyr::filter(!is.na(p.adj) ) %>%
    mutate( significance = cut(p.adj, 
                               breaks = c(0.0, 0.0001, 0.001, 0.01, 0.05, 1), 
                               labels = c('****', '***', '**', '*', 'ns')) )

  p_perf_f1_pwt <- ggplot( df_pwt, aes(x = id1, y = id2) ) + 
    geom_tile( aes(fill = significance), color = "white" ) +
    geom_text( aes(label = round(-log(p.adj))) ) +
    ggtitle( 'b) Pairwise Wilcoxon rank sum test',
             subtitle = 'values: -log(p.adj)') + 
    scale_fill_manual( values = c(rev(brewer.pal(4, 'Greens')), 'grey') ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  p <- grid.arrange(
    grobs = list(p_perf_f1_box, p_perf_f1_pwt), 
    layout_matrix = rbind(c(1,2), c(1,2))
  )
  
  return( p )
}

plot_perf_cvg_aux <- function ( df )
{
  # subplots for grid layout
  p_r_cvg <- ggplot(df, aes(x = as.factor(cvg), y = recall)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, caller) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = factor(cvg), y=mrec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth') + ggtitle( 'a' ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(.~caller) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = 'none', alpha = "none")
  
  p_p_cvg <- ggplot(df, aes(x = as.factor(cvg), y = precision)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, caller) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(cvg), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth') + ggtitle( 'b' ) +
    facet_grid(.~caller) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_cvg <- ggplot(df, aes(x = as.factor(cvg), y = F1)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, caller) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth', y = 'F1 score', fill = '') + ggtitle( 'c' ) +
    facet_grid(.~caller) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    theme(strip.text.x = element_text(size = 6)) +
    theme(legend.position = 'bottom') +
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_perf <- grid.arrange(
    grobs = list(p_r_cvg, p_p_cvg, p_f_cvg), 
    ncol = 1,
    nrow = 3,
    labels = "auto", 
    common.legend = TRUE, 
    legend = "bottom"
  )
  
  return( p_perf )
}