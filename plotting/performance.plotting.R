plot_perf_cvg <- function ( df )
{
  # exclude these callers from plots
  blacklist <- c( 'Strelka1' )
  
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
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
    'SNV-PPILP'))
  df$class <- factor(df$class, levels = c('marginal', 'joint', 'two-step'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'Mutect2\n(mseq)' = 'Mutect2_mseq',
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
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 7)) + 
    guides(colour = "none", fill = 'none', alpha = "none")
  
  p_p_cvg <- ggplot(df, aes(x = as.factor(cvg), y = precision)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(cvg), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth') + ggtitle( 'b' ) +
    facet_grid(.~lbl) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 7)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_cvg <- ggplot(df, aes(x = as.factor(cvg), y = F1)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'sequencing depth', y = 'F1 score') + ggtitle( 'c' ) +
    facet_grid(.~lbl) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    theme(strip.text.x = element_text(size = 7)) +
    theme(legend.position = 'bottom') +
    guides(colour = "none", alpha = "none")
  
  p_perf <- grid.arrange(
    grobs = list(p_r_cvg, p_p_cvg, p_f_cvg), 
    layout_matrix = rbind(c(1,1), c(2,2), c(3,3), c(3,3)),
    #top = textGrob("Performance of variant callers at different sequencing depths", gp = gpar(fontsize=16,font=3))
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
    'SNV-PPILP'))
  df$class <- factor(df$class, levels = c('marginal', 'joint', 'two-step'))
  
  # format caller names for better plotting
  df <- df %>% mutate(lbl = gsub("(?<=[a-z]{5}|-)([A-Z])", "\n\\1", df$caller, perl = T))
  df <- df %>% 
    mutate(lbl = fct_recode(caller, 
                            'Haplotype\nCaller' = 'HaplotypeCaller',
                            'Neu\nSomatic' = 'NeuSomatic',
                            'Somatic\nSniper' = 'SomaticSniper',
                            'Mutect2\n(mseq)' = 'Mutect2_mseq',
                            'SNV-\nPPILP' = 'SNV-PPILP'))
  
  # remove data for some callers (older versions of certain methods)
  df <- df %>% dplyr::filter( !(caller %in% blacklist) )
  
  # subplots for grid layout
  p_r_ttype <- ggplot(df, aes(x = as.factor(ttype), y = recall)) + 
    geom_boxplot(aes(alpha = ttype, fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = factor(ttype), y=mrec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture class') + ggtitle( 'a' ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(.~lbl) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 7)) + 
    guides(colour = "none", fill = 'none', alpha = "none")
  
  p_p_ttype <- ggplot(df, aes(x = as.factor(ttype), y = precision)) + 
    geom_boxplot(aes(alpha = factor(ttype), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(ttype), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture class') + ggtitle( 'b' ) +
    facet_grid(.~lbl) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          strip.text.x = element_text(size = 7)) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_ttype <- ggplot(df, aes(x = as.factor(ttype), y = F1)) + 
    geom_boxplot(aes(alpha = factor(ttype), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(ttype, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(ttype), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = 'admixture class', y = 'F1 score') + ggtitle( 'c' ) +
    facet_grid(.~lbl) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    theme(strip.text.x = element_text(size = 7)) +
    theme(legend.position = 'bottom') +
    guides(colour = "none", alpha = "none")
  
  p_perf <- grid.arrange(
    grobs = list(p_r_ttype, p_p_ttype, p_f_ttype), 
    layout_matrix = rbind(c(1,1), c(2,2), c(3,3), c(3,3)),
    #top = textGrob("Performance of variant callers at different sequencing depths", gp = gpar(fontsize=16,font=3))
  )
  
  return( p_perf )
}