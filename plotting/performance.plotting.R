plot_perf <- function ( df )
{
  # exclude these callers from plots
  blacklist <- c( 'Strelka1' )
  
  # define order of variant callers (will affect plots)
  df$caller = factor(df$caller, levels = c(
    "Bcftools", 
    "HaplotypeCaller", 
    "CaVEMan", 
    "MuTect1", 
    "Mutect2", 
    "NeuSomatic", 
    "Shimmer", 
    "SNooPer", 
    "SomaticSniper", 
    "Strelka1", 
    "Strelka2", 
    "VarDict", 
    "VarScan", 
    "MuClone", 
    "MultiSNV", 
    "Mutect2_mseq", 
    "SNV-PPILP"))
  df$class <- factor(df$class, levels = c("general", "tumor-normal", "multi-sample"))
  
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
  
  # smaller plots for grid layout
  p_r_cvg_m <- ggplot(df, aes(x = as.factor(cvg), y = recall)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, caller) %>% summarise(mrec = median(recall)) %>% arrange(desc(mrec)) %>% dplyr::filter(mrec==max(mrec)), aes(x = factor(cvg), y=mrec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = "sequencing depth") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(.~caller) +
    theme_gray() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_p_cvg_m <- ggplot(df, aes(x = as.factor(cvg), y = precision)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, caller) %>% summarise(mprec = median(precision)) %>% arrange(desc(mprec)) %>% dplyr::filter(mprec==max(mprec)), aes(x = factor(cvg), y = mprec), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = "sequencing depth") +
    facet_grid(.~caller) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank()) + 
    guides(colour = "none", fill = "none", alpha = "none")
  
  p_f_cvg_m <- ggplot(df, aes(x = as.factor(cvg), y = F1)) + 
    geom_boxplot(aes(alpha = factor(cvg), fill = class)) + ylim(0, 1) +
    geom_point(data = df %>% group_by(cvg, lbl) %>% summarise(mF1 = median(F1)) %>% arrange(desc(mF1)) %>% dplyr::filter(mF1==max(mF1)), aes(x=factor(cvg), y=mF1), fill = "gold", shape = 23) + 
    scale_alpha_manual(values = c(0.2, 0.5, 0.8, 1)) +
    labs(x = "sequencing depth") +
    facet_grid(.~lbl) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(strip.text.x = element_text(size = 7)) +
    theme(legend.position = "top") +
    guides(colour = "none", alpha = "none")
  
  p_perf <- grid.arrange(
    grobs = list(p_r_cvg_m, p_p_cvg_m, p_f_cvg_m), 
    layout_matrix = rbind(c(1,2), c(3,3), c(3,3)),
    #top = textGrob("Performance of variant callers at different sequencing depths", gp = gpar(fontsize=16,font=3))
  )
  
  return( p_perf )
}