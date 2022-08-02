require(rstatix)

plot_perf_min_OLD <- function ( df )
{
  # define order of variant callers (will affect plots)
  df = df %>% filter(caller!= "MuClone_perf")
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
    'Mutect2_multi_F',
    'Mutect2_MOSS'))
  df$class <- factor( df$class, levels = c('marginal', 'two-step', 'joint') )

  # subplots for grid layout
  p_rec <- ggplot( df, aes(x = caller, y = recall) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(mrec = median(recall, na.rm = T)) %>% dplyr::filter(mrec==max(mrec, na.rm = T)), aes(x = caller, y = mrec), fill = 'gold', shape = 23) + 
    stat_summary(fun='mean', shape=15) + # TODO: make changes permanent or remove
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    labs( x = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )

  p_pre <- ggplot( df, aes(x = caller, y = precision) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(mpre = median(precision, na.rm = T)) %>% dplyr::filter(mpre==max(mpre, na.rm = T)), aes(x = caller, y = mpre), fill = 'gold', shape = 23) + 
    stat_summary(fun='mean', shape=15) + # TODO: make changes permanent or remove
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    labs( x = 'caller' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )

  p_f1 <- ggplot( df, aes(x = caller, y = F1) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(mF1 = median(F1, na.rm = T)) %>% dplyr::filter(mF1==max(mF1, na.rm = T)), aes(x = caller, y = mF1), fill = 'gold', shape = 23) + 
    stat_summary(fun='mean', shape=15) + # TODO: make changes permanent or remove
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
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

plot_perf_min_mean <- function ( df )
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
    'Mutect2_multi_F',
    'Mutect2_Moss', 'Mutect2_MOSS', 'Strelka2_MOSS'))
  df$class <- factor( df$class, levels = c('marginal', 'two-step', 'joint') )
  
  # subplots for grid layout
  p_rec <- ggplot( df, aes(x = reorder(caller, -F1, na.rm = T), y = recall) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    labs( x = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )
  
  p_pre <- ggplot( df, aes(x = reorder(caller, -F1, na.rm = T), y = precision) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    labs( x = 'caller' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )
  
  p_f1 <- ggplot( df, aes(x = reorder(caller, -F1, na.rm = T), y = F1) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    labs( x = '', fill = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    theme( legend.position = 'right' )
  
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

plot_perf_min_sig <- function ( df )
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
    'Mutect2_multi_F',
    'Mutect2_MOSS'))
  df$class <- factor( df$class, levels = c('marginal', 'two-step', 'joint') )

  # reorder callers by overall F1 score
  df <- df %>% mutate( caller = fct_reorder(caller, -F1, .fun = mean, na.rm = T) )

  # perform group-wise comparisons for differences between distributions
  pairs <- cbind(head(levels(df$caller), -1), levels(df$caller)[-1]) %>%
    t() %>% as_tibble() %>% as.list()
  # gwc_rec_wilcox <- df %>%
  #   wilcox_test( ., recall ~ caller, comparisons = pairs ) %>%
  #   adjust_pvalue( method = 'bonferroni' ) %>%
  #   add_significance( p.col = 'p.adj', output.col = 'sig_sym' ) %>%
  #   add_xy_position( x = 'caller', step.increase = 0 ) %>%
  #   mutate( y.position = y.position + 0.02 )
  # gwc_pre_wilcox <- df %>%
  #   wilcox_test( ., precision ~ caller, comparisons = pairs ) %>%
  #   adjust_pvalue( method = 'bonferroni' ) %>%
  #   add_significance( p.col = 'p.adj', output.col = 'sig_sym' ) %>%
  #   add_xy_position( x = 'caller', step.increase = 0 ) %>%
  #   mutate( y.position = y.position + 0.02 )
  # gwc_f1_wilcox <- df %>%
  #   wilcox_test( ., F1 ~ caller, comparisons = pairs ) %>%
  #   adjust_pvalue( method = 'bonferroni' ) %>%
  #   add_significance( p.col = 'p.adj', output.col = 'sig_sym' ) %>%
  #   add_xy_position( x = 'caller', step.increase = 0 ) %>%
  #   mutate( y.position = y.position + 0.02 )

  # subplots for grid layout
  p_rec <- ggplot( df, aes(x = caller, y = recall) ) +
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) +
    #stat_pvalue_manual( data = gwc_f1_wilcox, label = 'p.adj.signif', tip.length = .01, size = 2 ) +
    scale_y_continuous( expand = expansion(mult = c(0.00, 0.05)) ) +
    labs( x = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )

  p_pre <- ggplot( df, aes(x = caller, y = precision) ) +
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) +
    #stat_pvalue_manual( data = gwc_f1_wilcox, label = 'p.adj.signif', tip.length = .01, size = 2 ) +
    scale_y_continuous( expand = expansion(mult = c(0.00, 0.05)) ) +
    labs( x = 'caller' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )

  p_f1 <- ggplot( df, aes(x = caller, y = F1) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    # stat_pvalue_manual( data = gwc_f1_wilcox, label = 'p.adj.signif', tip.length = .01, size = 2 ) +
    scale_y_continuous( expand = expansion(mult = c(0.00, 0.05)) ) +
    labs( x = '', fill = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    theme( legend.position = 'right' )

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
 
plot_perf_freq_OLD <- function ( df )
{
  # define order of variant callers (will affect plots)
  df = df %>% filter(caller!= "MuClone_perf")
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
  df$class <- factor( df$class, levels = c('marginal', 'two-step', 'joint') )
  
  # format caller names for better plotting
  #df <- df %>% mutate( lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T) )
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Mutect2\nmulti_F' = 'Mutect2_multi_F',
                            'Mutect2\nsingle' = 'Mutect2_single',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # reorder callers by overall F1 score
  df <- df %>% mutate( lbl = fct_reorder(lbl, -F1, .fun = mean, na.rm = T) )
  
  # subplots for grid layout
  p_rec <- ggplot( df, aes(x = freq_bin, y = recall, alpha = freq_bin) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class, middle = mean(recall)) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22, alpha = 1 ) +
    geom_point( data = df %>% group_by(freq_bin, lbl) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = freq_bin, y = m), color = 'gold', shape = 15, alpha = 1) + 
    labs( x = 'AF range')  +
    facet_wrap(~lbl, nrow = 1)+
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95, size = 6))+
    guides( fill = 'none' , alpha = 'none')
  
  p_pre <- ggplot( df, aes(x = freq_bin, y = precision, alpha = freq_bin) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class, middle = mean(precision)) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22, alpha = 1 ) +
    geom_point( data = df %>% group_by(freq_bin, lbl) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = freq_bin, y = m), color = 'gold', shape = 15, alpha = 1) + 
    labs( x = 'AF range')  +
    facet_wrap(~lbl, nrow = 1)+
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95, size = 6) ) +
    guides( fill = 'none', alpha = 'none' )
  
  p_f1 <- ggplot( df, aes(x = freq_bin, y = F1, alpha = freq_bin) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class, middle = mean(F1)) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22, alpha = 1 ) +
    geom_point( data = df %>% group_by(freq_bin, lbl) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = freq_bin, y = m), color = 'gold', shape = 15, alpha = 1) + 
    labs( x = 'AF range')  +
    facet_wrap(~lbl, nrow = 1)+
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95, size = 6) ) +
    guides(alpha = 'none')+
    theme(legend.position = 'right')
  
  # theme( legend.position = 'bottom' ) +
  # guides( colour = 'none' )+
  # scale_fill_discrete(name="")
  
  p_perf <- ggarrange(
    p_rec, p_pre, p_f1,
    ncol = 1,
    nrow = 3,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "top"
  )
  
  return( p_perf )
}

plot_perf_freq <- function ( df )
{
  # define order of variant callers (will affect plots)
  df = df %>% filter(caller!= "MuClone_perf")
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
  
  # reorder callers by overall F1 score
  df <- df %>% mutate( lbl = fct_reorder(lbl, -recall, .fun = mean, na.rm = T) )
  
  # subplots for grid layout
  p_perf <- ggplot( df, aes(x = freq_bin, y = recall, alpha = freq_bin) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class, middle = mean(recall)) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22, alpha = 1 ) +
    geom_point( data = df %>% group_by(freq_bin, lbl) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = freq_bin, y = m), color = 'gold', shape = 15, alpha = 1) + 
    labs( x = 'AF range')  +
    facet_wrap(~lbl, nrow = 1)+
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95, size = 6))+
    guides( fill = 'none' , alpha = 'none')+
    xlab("Simulated AF")
  
  
  
  # theme( legend.position = 'bottom' ) +
  # guides( colour = 'none' )+
  # scale_fill_discrete(name="")
  
 
  
  
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
    'Mutect2_multi_F',
    'SNooPerGermres',
    'SNooPerGermres.70'))

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
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mrec = median(recall, na.rm = T)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec, na.rm = T)), aes(x = factor(cvg), y=mrec), fill = "gold", shape = 23) + 
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
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mprec = median(precision, na.rm = T)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec, na.rm = T)), aes(x = factor(cvg), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth')  +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_cvg <- ggplot(df, aes(x = as.factor(cvg), y = F1)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1, na.rm = T)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1, na.rm = T)), aes(x=factor(cvg), y=mF1), fill = "gold", shape = 23) + 
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
plot_perf_cvg_sig <- function ( df )
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
  sig_legend = sprintf('Group-wise comparison of means by factor; method: Kruskal-Wallis; p.adjust: "bonferroni"')
  sig_legend = paste(sig_legend, sprintf('(%s)', paste(sprintf('%s: p<%s', sig_symbols, sig_cutpoints[-1]), collapse = ', ')), sep = ' ')
  
  # perform group-wise comparisons for differences between distributions
  gwc_rec_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., recall ~ cvg) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  gwc_pre_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., precision ~ cvg) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  gwc_f1_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., F1 ~ cvg) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  # reorder callers by overall F1 score
  df <- df %>% mutate( lbl = fct_reorder(lbl, -F1, .fun = mean, na.rm = T) )
  
  # subplots for grid layout
  p_r_cvg <- ggplot( df, aes(x = factor(cvg), y = recall), yticks.by = 0.25) +
    geom_boxplot( aes(fill = class) ) +
    facet_wrap( ~lbl, nrow = 1 ) +
    labs( x = 'sequencing depth' ) +
    theme_minimal() +
    rotate_x_text() +
    font( 'xy.text', size = 8 ) +
    theme( strip.text.x = element_text(family='Helvetica-Narrow', size = 8) ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = factor(cvg), y = m), color = 'gold', shape = 15) + 
    geom_text( data = gwc_rec_kruskal, aes(label = sig_sym, y = 1.05, x = 2.5) )
  
  p_p_cvg <- ggplot( df, aes(x = factor(cvg), y = precision), yticks.by = 0.25) +
    geom_boxplot( aes(fill = class) ) +
    facet_wrap( ~lbl, nrow = 1 ) +
    labs( x = 'sequencing depth' ) +
    theme_minimal() +
    rotate_x_text() +
    font( 'xy.text', size = 8 ) +
    theme( strip.text.x = element_text(family='Helvetica-Narrow', size = 8) ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = factor(cvg), y = m), color = 'gold', shape = 15) + 
    geom_text( data = gwc_pre_kruskal, aes(label = sig_sym, y = 1.05, x = 2.5) )
  
  p_f_cvg <- ggplot( df, aes(x = factor(cvg), y = F1), yticks.by = 0.25) +
    geom_boxplot( aes(fill = class) ) +
    facet_wrap( ~lbl, nrow = 1 ) +
    labs( x = 'sequencing depth' ) +
    theme_minimal() +
    rotate_x_text() +
    font( 'xy.text', size = 8 ) +
    theme( strip.text.x = element_text(family='Helvetica-Narrow', size = 8) ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(cvg, lbl) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = factor(cvg), y = m), color = 'gold', shape = 15) + 
    geom_text( data = gwc_f1_kruskal, aes(label = sig_sym, y = 1.05, x = 2.5) )
  
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
    bottom = text_grob(sig_legend, hjust = 1, x = 1, face = "italic", size = 8))
  
  return( p_perf )
}

# plot performance at different depth of coverage
# test significance of difference between depth levels for each caller
# params:
#  df     - caller, recall, precision, F1, ttype
plot_perf_admix_sig <- function ( df )
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
  sig_legend = sprintf('Group-wise comparison of means by factor; method: Kruskal-Wallis; p.adjust: "bonferroni"')
  sig_legend = paste(sig_legend, sprintf('(%s)', paste(sprintf('%s: p<%s', sig_symbols, sig_cutpoints[-1]), collapse = ', ')), sep = ' ')
  
  # perform group-wise comparisons for different means
  gwc_rec_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., recall ~ ttype) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  gwc_pre_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., precision ~ ttype) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  gwc_f1_kruskal <- df %>% group_by( lbl ) %>%
    kruskal_test(., F1 ~ ttype) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' )
  
  # reorder callers by overall F1 score
  df <- df %>% mutate( lbl = fct_reorder(lbl, -F1, .fun = mean, na.rm = T) )
  
  # subplots for grid layout
  p_r_mix <- ggboxplot(df, x = 'ttype', y = 'recall', fill = 'class', yticks.by = 0.25) +
    #geom_text( data = gwc_rec, aes(label = sig_sym, y = 1.05, x = 2) ) +
    facet_wrap( ~lbl, nrow = 1) +
    labs( x = 'admixture' ) +
    theme_minimal() +
    rotate_x_text() +
    font( 'xy.text', size = 8 ) +
    theme( strip.text.x = element_text(family='Helvetica-Narrow', size = 8) ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(ttype, lbl) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = ttype, y = m), color = 'gold', shape = 15) + 
    geom_text( data = gwc_rec_kruskal, aes(label = sig_sym, y = 1.05, x = 2) )
  
  p_p_mix <- ggboxplot(df, x = 'ttype', y = 'precision', fill = 'class', yticks.by = 0.25) %>%
    facet(facet.by = 'lbl', nrow = 1) +
    labs( x = 'admixture' ) +
    theme_minimal() +
    rotate_x_text() +
    font( 'xy.text', size = 8 ) +
    theme( strip.text.x = element_text(family='Helvetica-Narrow', size = 8) ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(ttype, lbl) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = ttype, y = m), color = 'gold', shape = 15) + 
    geom_text( data = gwc_pre_kruskal, aes(label = sig_sym, y = 1.05, x = 2) )
  
  p_f_mix <- ggboxplot(df, x = 'ttype', y = 'F1', fill = 'class', yticks.by = 0.25) %>%
    facet(facet.by = 'lbl', nrow = 1) +
    labs( x = 'admixture' ) +
    theme_minimal() +
    rotate_x_text() +
    font( 'xy.text', size = 8 ) +
    theme( strip.text.x = element_text(family='Helvetica-Narrow', size = 8) ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(ttype, lbl) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = ttype, y = m), color = 'gold', shape = 15) + 
    geom_text( data = gwc_f1_kruskal, aes(label = sig_sym, y = 1.05, x = 2) )
  
  # combine subplots
  p_perf <- ggarrange(
    p_r_mix, p_p_mix, p_f_mix,
    ncol = 1,
    nrow = 3,
    labels = "auto", 
    common.legend = TRUE, 
    legend = "bottom"
  ) %>%
    annotate_figure(
      bottom = text_grob(sig_legend, hjust = 1, x = 1, face = "italic", size = 8))
  
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
    'Mutect2_multi_F',
    'SNooPerGermres',
    'SNooPerGermres.70'))
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
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mrec = median(recall, na.rm = T)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec, na.rm = T)), aes(x = factor(ttype), y=mrec), fill = "gold", shape = 23) + 
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
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mprec = median(precision, na.rm = T)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec, na.rm = T)), aes(x = factor(ttype), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture') + 
    facet_grid(.~lbl) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_ttype <- ggplot(df, aes(x = as.factor(ttype), y = F1)) + 
    geom_boxplot(aes(alpha = factor(ttype), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mF1 = median(F1, na.rm = T)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1, na.rm = T)), aes(x=factor(ttype), y=mF1), fill = "gold", shape = 23) + 
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

plot_pairwise_wilcoxon <- function( df ) {
  # plot Box- and Violinplots of F1 scores for callers
  p_perf_f1_box <- ggviolin( df, x = 'caller', y = 'F1', color = 'class',
                             main = 'a) Distribution of F1 scores',
                             add = 'boxplot' ) +
    rotate_x_text( 45 )
  
  # perform Wilcoxon rank sum test (Mann-Whitney test)
  pwt <- pairwise.wilcox.test( df$F1, df$caller, p.adjust.method = "BH" )
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
    layout_matrix = rbind(c(1,1), c(1,1), c(2,2), c(2,2), c(2,2))
  )
  
  return( p )
}

# Create subplots for recall, precision, F1 score.
# Returns: 
#   vector of subplots c(p_rec, p_pre, p_f1)
#   to be combined in subsequent step
plot_perf_min_internal <- function ( df )
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
    'Mutect2_multi_F',
    'Mutect2_MOSS'))
  df$class <- factor( df$class, levels = c('marginal', 'two-step', 'joint') )
  
  # reorder callers by overall F1 score
  df <- df %>% mutate( caller = fct_reorder(caller, -F1, .fun = mean, na.rm = T) )
  
  # perform group-wise comparisons for differences between distributions
  pairs <- cbind(head(levels(df$caller), -1), levels(df$caller)[-1]) %>% 
    t() %>% as_tibble() %>% as.list()
  gwc_f1_wilcox <- df %>%
    wilcox_test( ., F1 ~ caller, comparisons = pairs ) %>%
    adjust_pvalue( method = 'bonferroni' ) %>%
    add_significance( p.col = 'p.adj', output.col = 'sig_sym' ) %>%
    add_xy_position( x = 'caller', step.increase = 0 ) %>%
    mutate( y.position = y.position + 0.02 )

  # subplots for grid layout
  p_rec <- ggplot( df, aes(x = caller, y = recall) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    scale_y_continuous( expand = expansion(mult = c(0.00, 0.05)) ) +
    labs( x = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )
  
  p_pre <- ggplot( df, aes(x = caller, y = precision) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    scale_y_continuous( expand = expansion(mult = c(0.00, 0.05)) ) +
    labs( x = 'caller' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    guides( fill = 'none' )
  
  p_f1 <- ggplot( df, aes(x = caller, y = F1) ) + 
    theme_minimal() +
    geom_boxplot( aes(fill = class) ) + ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_point( data = df %>% group_by(caller) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = caller, y = m), color = 'gold', shape = 15) + 
    stat_pvalue_manual( data = gwc_f1_wilcox, label = 'p.adj.signif', tip.length = .01, size = 2 ) +
    scale_y_continuous( expand = expansion(mult = c(0.00, 0.05)) ) +
    labs( x = '', fill = '' )  +
    theme( axis.text.x = element_text(angle = 45, hjust = 0.95) ) +
    theme( legend.position = 'right' )

  p_perf <- ggarrange(
    p_rec, p_pre, p_f1,
    ncol = 3,
    nrow = 1,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "top"
  )
  
  return( list(p_rec, p_pre, p_f1, p_perf) )  
}

plot_perf_min <- function ( df1, df2 )
{
  subplots1 <- plot_perf_min_internal( df1 )
  subplots2 <- plot_perf_min_internal( df2 )
  
  p_perf <- ggarrange(
    subplots1[[1]], subplots1[[2]], subplots1[[3]],
    subplots2[[1]], subplots2[[2]], subplots2[[3]],
    ncol = 3,
    nrow = 2,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "top"
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
  
  # reorder callers by mean F1 score
  df$lbl <- factor( df$lbl, levels = df %>% group_by(lbl) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% arrange(desc(m)) %>% pull(lbl) )
  
  # subplots for grid layout
  p_rec <- ggplot( df, aes(x = as.factor(cvg), y = recall) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = group) ) + 
    ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_hline( color = 'red', yintercept = df %>% dplyr::filter(group == 'single') %>% summarise(m = mean(recall, na.rm = T)) %>% pull(m) ) + 
    geom_point( data = df %>% group_by(cvg, lbl) %>% dplyr::summarise(m = mean(recall, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = factor(cvg), y = m), color = 'gold', shape = 15) + 
    labs( x = 'caller', fill = '' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank()) + 
    scale_fill_manual( name = '', values = c('tomato', 'white') ) +
    guides( colour = "none", fill = 'none' )
  
  p_pre <- ggplot( df, aes(x = as.factor(cvg), y = precision) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = group) ) + 
    ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_hline( color = 'red', yintercept = df %>% dplyr::filter(group == 'single') %>% summarise(m = mean(precision, na.rm = T)) %>% pull(m) ) + 
    geom_point( data = df %>% group_by(cvg, lbl) %>% dplyr::summarise(m = mean(precision, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = factor(cvg), y = m), color = 'gold', shape = 15) + 
    labs( x = 'caller' ) + 
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) + 
    scale_fill_manual( name = '', values = c('tomato', 'white') ) +
    guides( colour = "none", fill = 'none' )
  
  p_f1 <- ggplot( df, aes(x = as.factor(cvg), y = F1) ) + 
    theme_minimal()+
    geom_boxplot( aes(fill = group) ) + 
    ylim( 0, 1 ) +
    stat_summary( fun = 'mean', fill = 'white', shape = 22 ) +
    geom_hline( color = 'red', yintercept = df %>% dplyr::filter(group == 'single') %>% summarise(m = mean(F1, na.rm = T)) %>% pull(m) ) + 
    geom_point( data = df %>% group_by(cvg, lbl) %>% dplyr::summarise(m = mean(F1, na.rm = T)) %>% dplyr::filter(m==max(m, na.rm = T)), aes(x = factor(cvg), y = m), color = 'gold', shape = 15) + 
    labs( x = 'caller', y = 'F1 score' )  +
    facet_wrap( ~ lbl, nrow = 1 ) +
    theme( strip.text.x = element_text(size = 6) ) +
    theme( axis.text.x = element_blank(),
           axis.ticks.x = element_blank() ) +
    #theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
    theme( legend.position = 'bottom' ) +
    scale_fill_manual( name = '', values = c('tomato', 'white') ) +
    guides( colour = 'none' )
    
  p_perf <- ggarrange(
    p_rec, p_pre, p_f1,
    ncol = 1,
    nrow = 3,
    labels = "auto" , 
    common.legend = TRUE, 
    legend = "bottom"
  )
  
  return( p_perf )
}

plot_perf_pairs_heatmap <- function ( df )
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

plot_perf_trip_heatmap <- function ( df )
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
  
  return( p )
}
