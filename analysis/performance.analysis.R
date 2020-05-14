classify_variants <- function ( 
  df_varcall, 
  df_mut, 
  df_mut_sample,
  df_caller
)
{
  # join detected variants and true mutations
  df_tp <- df_varcall %>%
    inner_join( df_mut, by = c('id_rep', 'chrom', 'pos') ) %>%
    inner_join( df_mut_sample %>% dplyr::filter( as.logical( is_present ) ), 
                by = c('id_rep', 'id_sample', 'id_mut') ) %>%
    select( id_caller, id_rep, id_sample, chrom, pos ) %>%
    mutate( type = 'TP' )
  
  # join detected variants and true mutations
  df_fp <- df_varcall %>% 
    anti_join( df_mut_sample %>% dplyr::filter( as.logical( is_present ) ) %>%
                 inner_join( df_mut, by = c('id_rep', 'id_mut') ), 
               by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
    select( id_caller, id_rep, id_sample, chrom, pos ) %>%
    mutate( type = 'FP' )
  
  # join detected variants and true mutations
  df_fn <- df_mut_sample %>% dplyr::filter( as.logical( is_present ) ) %>%
    inner_join( df_mut, by = c('id_rep', 'id_mut') ) %>%
    crossing( tibble( id_caller = df_caller$id_caller ) ) %>%
    anti_join( df_varcall, by = c('id_rep', 'id_caller', 'id_sample', 'chrom', 'pos') ) %>%
    select( id_caller, id_rep, id_sample, chrom, pos ) %>%
    mutate( type = 'FN' )
  
  df_vars <- df_tp %>% rbind( df_fp ) %>% rbind( df_fn )
  
  return( df_vars )
}

calculate_performance <- function (
  df_vars, 
  df_caller, 
  df_rep
)
{
  df_eval <- df_vars %>% select( id_caller, id_rep, type ) %>%
    group_by( id_caller, id_rep, type ) %>%
    summarise( n = n() ) %>%
    ungroup() %>%
    complete( id_caller, id_rep, type, fill = list(n = 0) ) %>%
    spread( type, n ) %>%
    mutate( recall = TP/(TP+FN), precision = TP/(TP+FP) ) %>%
    mutate( F1 = 2*recall*precision/(recall+precision) )
  # add caller information
  df_eval <- df_eval %>%
    inner_join( df_caller, by = c('id_caller') ) %>%
    inner_join( df_rep, by = c('id_rep') ) %>%
    replace_na( list(precision = 0, F1 = 0) ) %>%
    mutate( caller = name_caller )
  
  return( df_eval )
}

calculate_performance_sample <- function (
  df_vars, 
  df_caller, 
  df_rep
)
{
  df_eval <- df_vars %>% select( id_caller, id_rep, id_sample, type ) %>%
    group_by( id_caller, id_rep, id_sample, type ) %>%
    summarise( n = n() ) %>%
    ungroup() %>%
    complete( id_caller, id_rep, id_sample, type, fill = list(n = 0) ) %>%
    spread( type, n ) %>%
    mutate( recall = TP/(TP+FN), precision = TP/(TP+FP) ) %>%
    mutate( F1 = 2*recall*precision/(recall+precision) )
  # add caller information
  df_eval <- df_eval %>%
    inner_join( df_caller, by = c('id_caller') ) %>%
    inner_join( df_rep, by = c('id_rep') ) %>%
    replace_na( list(precision = 0, F1 = 0) ) %>%
    mutate( caller = name_caller )
  
  return( df_eval )
}

performance_kruskal_wallis <- function ( df_perf )
{
  require(rstatix) # wilcox_test()
  
  kruskwal <- df_perf %>%
    group_by( name_caller ) %>%
    kruskal_test( F1 ~ cvg ) %>%
    adjust_pvalue( method = 'BH' )
  # Pairwise comparisons between condition levels
  # (this is overkill, though...)
  # wilcox <- df %>%
  #   group_by(lbl) %>%
  #   wilcox_test(F1 ~ cvg) %>%
  #   adjust_pvalue(method = 'BH') %>%
  #   mutate(y.position = case_when(
  #     (group1 == '30' & group2 == '50') ~ 1.02,
  #     (group1 == '30' & group2 == '100') ~ 1.10,
  #     (group1 == '30' & group2 == '300') ~ 1.18,
  #     (group1 == '50' & group2 == '100') ~ 1.06,
  #     (group1 == '50' & group2 == '300') ~ 1.14,
  #     (group1 == '100' & group2 == '300') ~ 1.02,
  #     TRUE ~ 1.35
  #   ))
  # 
  # p_f_cvg <- ggboxplot(df, x = 'cvg', y = 'F1', fill = 'class', yticks.by = 0.25) %>%
  #   facet(facet.by = 'lbl', nrow = 1) +
  #   rotate_x_text() +
  #   #ylim(0.0, 1.16) +
  #   stat_pvalue_manual(wilcox, label = 'p.adj.signif', tip.length = 0.01) # Add pairwise comparisons p-value
  #     #font("xy.text", size = 12)
  # ggsave(plot = p_f_cvg, filename = 'de-novo.f1.wilcox.pdf', width = 12, height = 9)
  # ggsave(plot = p_f_cvg, filename = 'de-novo.f1.wilcox.png', width = 12, height = 9)
  
}

performance_anova <- function ( df_perf )
{
  require(rstatix) # wilcox_test()
  
  kruskwal <- df_perf %>%
    group_by( name_caller ) %>%
    kruskal_test( F1 ~ cvg ) %>%
    adjust_pvalue( method = 'BH' )
  # Pairwise comparisons between condition levels
  # (this is overkill, though...)
  # wilcox <- df %>%
  #   group_by(lbl) %>%
  #   wilcox_test(F1 ~ cvg) %>%
  #   adjust_pvalue(method = 'BH') %>%
  #   mutate(y.position = case_when(
  #     (group1 == '30' & group2 == '50') ~ 1.02,
  #     (group1 == '30' & group2 == '100') ~ 1.10,
  #     (group1 == '30' & group2 == '300') ~ 1.18,
  #     (group1 == '50' & group2 == '100') ~ 1.06,
  #     (group1 == '50' & group2 == '300') ~ 1.14,
  #     (group1 == '100' & group2 == '300') ~ 1.02,
  #     TRUE ~ 1.35
  #   ))
  # 
  # p_f_cvg <- ggboxplot(df, x = 'cvg', y = 'F1', fill = 'class', yticks.by = 0.25) %>%
  #   facet(facet.by = 'lbl', nrow = 1) +
  #   rotate_x_text() +
  #   #ylim(0.0, 1.16) +
  #   stat_pvalue_manual(wilcox, label = 'p.adj.signif', tip.length = 0.01) # Add pairwise comparisons p-value
  #     #font("xy.text", size = 12)
  # ggsave(plot = p_f_cvg, filename = 'de-novo.f1.wilcox.pdf', width = 12, height = 9)
  # ggsave(plot = p_f_cvg, filename = 'de-novo.f1.wilcox.png', width = 12, height = 9)
  
}