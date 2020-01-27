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
