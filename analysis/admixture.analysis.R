get_shannon_sample <- function( df_prev )
{
  df <- df_prev %>% 
    dplyr::filter( id_sample!='RN', prev>0 ) %>% 
    group_by( id_rep, id_sample ) %>% 
    summarise( H = -sum(prev*log(prev)) ) %>%
    ungroup()
  
  return( df )
}