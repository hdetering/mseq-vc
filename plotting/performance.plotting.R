plot_perf_cvg <- function ( df )
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

plot_pairwise_wilcoxon <- function( df ) {
  # plot Box- and Violinplots of F1 scores for callers
  p_perf_f1_box <- ggviolin( df, x = 'caller', y = 'F1', color = 'class',
                             main = 'a) Distribution of F1 scores',
                             add = 'boxplot' ) +
    rotate_x_text( 45 )
  
  # perform pairwise Wilcoxon rank sum test (Mann-Whitney test?)
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
    ggtitle( 'b) Pairwise Wilcoxon rank sum test (values: -log(p.adj))' ) + 
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
    'MuClone_perf', 
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
                            'Mutect2\nmulti' = 'Mutect2_multi',
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
    labs(x = 'sequencing depth') + ggtitle( 'a' ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = 'none', alpha = "none")
  
  p_p_cvg <- ggplot(df, aes(x = as.factor(cvg), y = precision)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(cvg), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth') + ggtitle( 'b' ) +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 6)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_cvg <- ggplot(df, aes(x = as.factor(cvg), y = F1)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth', y = 'F1 score', fill = '') + ggtitle( 'c' ) +
    facet_grid(.~lbl) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    theme(strip.text.x = element_text(size = 6)) +
    theme(legend.position = 'bottom') +
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_perf <- grid.arrange(
    grobs = list(p_r_cvg, p_p_cvg, p_f_cvg), 
    layout_matrix = rbind(c(1,1), c(2,2), c(3,3), c(3,3)),
    #top = textGrob("Performance of variant callers at different sequencing depths", gp = gpar(fontsize=16,font=3))
  )
  
  return( p_perf )
}