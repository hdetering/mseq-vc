plot_rc_alt_bar_srsv <- function( df_vars, df_rc, df_rep ) {
  
  df <- df_vars %>% dplyr::filter( name_caller != 'Strelka1' ) %>% 
    inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
    select(caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt) %>%
    mutate( vaf = (rc_alt)/(rc_ref+rc_alt) ) %>%
    inner_join( df_rep, by = 'id_rep' )
  df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )
  
  p_rc_alt <- ggplot(df, aes(x = vaf, fill = type)) + 
    geom_histogram( position = 'dodge', bins = 50 ) +
    scale_fill_brewer( palette = 'Set1' ) +
    facet_grid( caller~cvg, scales = 'free_y' ) +
    ggtitle( 'VAF distribution for classes of variant calls' ) + xlab( 'variant allele frequency' ) +
    theme_grey() +
    theme(
      strip.text.y = element_text(angle=0),
      legend.position = 'bottom'
    )
}

plot_rc_alt_ridges_srsv <- function( df_vars, df_rc, df_rep ) {
  
  require(ggridges)
  
  df <- df_vars %>% dplyr::filter( name_caller != 'Strelka1' ) %>% 
    inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
    select(caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt) %>%
    inner_join( df_rep, by = 'id_rep' ) %>%
    group_by(caller, cvg, type, rc_alt) %>%
    summarize( n = n() )
  df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )
  
  p_rc_alt <- ggplot(df, aes(x = rc_alt, y = type)) + 
    geom_density_ridges(aes(height = n, fill = type), stat="identity", alpha = 0.4) +
    scale_fill_brewer(palette = 'Set1' ) +
    facet_grid(caller~cvg, scales = 'free_x') +
    #ggtitle('ALT allele read count distribution for classes of variant calls') + 
    xlab( 'alternative allele read count' ) +
    theme_grey() +
    theme(
      strip.text.y = element_text(angle=0),
      legend.position = 'bottom'
    )
}

plot_rc_alt_ridges_rrsv <- function( df_vars, df_rc, df_rep ) {
  
  require(ggridges)
  
  # df <- df_vars %>% dplyr::filter( name_caller != 'Strelka1' ) %>%
  #   inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
  #   select(caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt) %>%
  #   inner_join( df_rep, by = 'id_rep' ) %>%
  #   group_by(caller, cvg, type, rc_alt) %>%
  #   summarize( n = n() )
  df <- df_vars %>% dplyr::filter( name_caller != 'Strelka1' ) %>%
    inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
    select( caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt ) %>%
    mutate( vaf = (rc_alt+1)/(rc_ref+rc_alt+1) ) %>%
    inner_join( df_rep, by = 'id_rep' )
  df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )
  
  p_rc_alt <- ggplot(df, aes(x = vaf, y = type)) + 
    geom_density_ridges( aes(fill = type), alpha = 0.4 ) +
    scale_fill_brewer( palette = 'Set1' ) +
    facet_wrap( ~caller, scales = 'free_y', ncol = 4 ) +
    ggtitle('VAF distribution for classes of variant calls') + xlab( 'variant allele frequency' ) +
    theme_grey() +
    theme(
      strip.text.y = element_text(angle=0),
      legend.position = 'bottom'
    )
}