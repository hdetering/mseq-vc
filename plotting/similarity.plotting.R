#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Plotting functions to visualize overlaps between variant callers.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-05-09
#------------------------------------------------------------------------------

# plot single similarity matrix
plot_jacc_idx <- function( df ) {
  p_jacc <- df %>% dplyr::filter(caller1 != 'Strelka1' & caller2 != 'Strelka1') %>%
    ggplot( aes(x = caller2, y = caller1)) + geom_tile(aes(fill = jaccard_idx) ) +
    scale_fill_distiller( palette = "Spectral", limits = c(0, 1), name = 'Jaccard\nindex' ) +
    theme_bw() +
    theme( axis.title = element_blank(),
           axis.text.x = element_text(angle = 30, hjust = 1, size = 6) )
  
  return( p_jacc )
}

# combine three plots into one
plot_jacc_idx_multi <- function( p1, p2, p3 ) {
  
  require( cowplot ) # get_legend
  
  p_jacc_multi_graphs <- plot_grid( 
    p1 + theme(legend.position = 'none'),
    p2 + theme(legend.position = 'none'),
    p3 + theme(legend.position = 'none'),
    labels = c('a', 'b', 'c'), nrow = 1
  )
  p_jacc_multi_legend <- get_legend( p1 + theme(legend.position = 'bottom') )
  p_jacc_multi <- plot_grid(
    p_jacc_multi_graphs, 
    p_jacc_multi_legend, 
    ncol = 1, rel_heights = c(1, .2) )
  
  return( p_jacc_multi )
}