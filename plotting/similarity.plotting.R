#!/usr/bin/env Rscript
# vim: syntax=R tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Plotting functions to visualize overlaps between variant callers.
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-01-10
#------------------------------------------------------------------------------

# plot single similarity matrix
plot_jacc_idx <- function( df ) {
  p_jacc <- df %>% dplyr::filter(caller1 != 'Strelka1' & caller2 != 'Strelka1') %>%
    ggplot( aes(x = caller2, y = caller1)) + geom_tile(aes(fill = jaccard_idx) ) +
    scale_fill_distiller( palette = "Spectral", limits = c(0, 1), name = 'Jaccard\nindex' ) +
    theme_bw() +
    theme( axis.title = element_blank(),
           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) ,
           plot.margin = margin(t = 20, l = 10))
  
  return( p_jacc )
}

# combine three plots into one
plot_jacc_idx_multi <- function( p1, p2, p3 ) {
  
  # require( cowplot ) # get_legend
  # 
  # p_jacc_multi_graphs <- plot_grid( 
  #   p1 + theme(legend.position = 'none'),
  #   p2 + theme(legend.position = 'none'),
  #   p3 + theme(legend.position = 'none'),
  #   labels = c('a', 'b', 'c'), nrow = 1
  # )
  # p_jacc_multi_legend <- get_legend( p1 + theme(legend.position = 'bottom') )
  # p_jacc_multi <- plot_grid(
  #   p_jacc_multi_graphs, 
  #   p_jacc_multi_legend, 
  #   ncol = 1, rel_heights = c(1, .2) )
  
  p_jacc_multi = ggarrange( p1, p2, p3,
                            labels = "auto",
                            ncol = 3,
                            nrow = 1,
                            common.legend = TRUE,
                            legend = "bottom" )
  return( p_jacc_multi )
}

plot_dendrogram <- function( df_jacc_idx ) {
  require(ggdendro) # dendro_data()
  
  # reformat data frame into
  df <- df_jacc_idx %>% 
    pivot_wider( names_from = caller1, values_from = jaccard_idx ) %>% 
    column_to_rownames( 'caller2' )
  # convert Jaccard similarity indices into distances
  d <- as.dist( 1-df )
  # c;uster callers by distance
  hc <- hclust(d)
  # plot dendrogram
  dd <- dendro_data( hc, type = 'rectangle' )
  p <- ggplot() +
    geom_segment( data = segment(dd), aes(x, y, xend = xend, yend = yend) ) +
    geom_text( data = label(dd), aes(x, y, label = label, hjust = 0), size = 3 ) +
    coord_flip() + scale_y_reverse( expand = expand_scale(add = c(0, 0.35)) ) +
    theme_dendro()
  
  p
}