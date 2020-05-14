#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Analysis functions to quantify overlaps between variant callers.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-05-09
#------------------------------------------------------------------------------

# calculate the Jaccard Index between two binary vectors
Jaccard <- function (x, y) {
  M.11 <- sum( x == 1 & y == 1 )
  M.10 <- sum( x == 1 & y == 0 )
  M.01 <- sum( x == 0 & y == 1 )
  if ( (M.11 + M.10 + M.01) == 0 )
    return ( 1 )
  else
    return ( M.11 / (M.11 + M.10 + M.01) )
}

# calculate the Jaccard Index for a data frame of multiple callers
Jaccard.df <- function (data) {
  lbl <- names( data )
  df_jacc <- tibble(
    caller1 = character(),
    caller2 = character(),
    jaccard_idx = numeric()
  )
  for ( r in 1:length(lbl) ) {
    for ( c in 1:length(lbl) ) {
      if ( c == r ) {
        df_jacc <- df_jacc %>% add_row(
          caller1 = lbl[r], 
          caller2 = lbl[c], 
          jaccard_idx = 1.0
        )
      }
      else if ( c > r ) {
        jacc_idx <- Jaccard( data[, r], data[, c] )
        df_jacc <- df_jacc %>% add_row(
          caller1 = lbl[r], 
          caller2 = lbl[c],
          jaccard_idx = jacc_idx
        )
      }
    }
  }
  return( df_jacc )
}

Jaccard.df.lbl <- function (data, id) {
  lbl <- names( data )
  df_jacc <- tibble(
    label = character(),
    caller1 = character(),
    caller2 = character(),
    jaccard_idx = numeric()
  )
  for ( r in 1:length(lbl) ) {
    for ( c in 1:length(lbl) ) {
      if ( c == r ) {
        df_jacc <- df_jacc %>% add_row(
          label = id,
          caller1 = lbl[r], 
          caller2 = lbl[c], 
          jaccard_idx = 1.0
        )
      }
      else if ( c > r ) {
        jacc_idx <- Jaccard( data[, r], data[, c] )
        df_jacc <- df_jacc %>% add_row(
          label = id,
          caller1 = lbl[r], 
          caller2 = lbl[c],
          jaccard_idx = jacc_idx
        )
      }
    }
  }
  return( df_jacc )
}

# get presence / absence info for given variant type
get_var_pres <- function( df_var, df_caller, var_type ) 
{
  df_var_pres <- df_var %>% dplyr::filter( type == var_type ) %>% 
    inner_join( df_caller, by = 'id_caller' ) %>%
    mutate(id_mut = sprintf("%s:%s:%s:%d", id_rep, id_sample, chrom, pos), present = 1) %>%
    select(name_caller, id_mut, present) %>%
    spread(name_caller, present, fill = 0)
  
  return( df_var_pres )
}
