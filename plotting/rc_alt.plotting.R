plot_vaf_bar_srsv <- function( df_vars, df_rc, df_rep ) {
  
  df <- df_vars %>% dplyr::filter( name_caller != 'Strelka1' ) %>% 
    inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
    select(caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt) %>%
    mutate( vaf = if_else(rc_ref+rc_alt==0, 0, (rc_alt)/(rc_ref+rc_alt) )) %>%
    inner_join( df_rep, by = 'id_rep' )
  df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )
  
  p_rc_alt <- ggplot(df, aes(x = vaf, fill = type)) + 
    geom_histogram( position = 'dodge', bins = 50 ) +
    scale_fill_brewer( palette = 'Set1' ) +
    facet_wrap( ~caller, ncol = 4, scales = 'free_y' ) +
    #ggtitle( 'VAF distribution for classes of variant calls' ) + 
    labs( x = 'variant allele frequency', fill = '' ) +
    theme_pubclean() +
    theme( legend.position = 'bottom' )
}

plot_vaf_bar_empirical <- function( df_vars, df_rc ) {
  
  df <- df_vars %>%
    inner_join( df_rc, by = c('id_sample', 'chrom', 'pos') ) %>%
    select( caller = name_caller, id_sample, chrom, pos, rc_ref, rc_alt) %>%
    mutate( vaf = if_else(rc_ref+rc_alt==0, 0, (rc_alt)/(rc_ref+rc_alt) ))
  
  p_rc_alt <- ggplot( df, aes(x = vaf) ) + 
    geom_histogram( bins = 50 ) +
    facet_wrap( ~caller, ncol = 4, scales = 'free_y' ) +
    #ggtitle( 'VAF distribution for classes of variant calls' ) + 
    labs( x = 'variant allele frequency' ) +
    theme_pubclean()
}

plot_vaf_bar_empirical_sample <- function( df_vars, df_rc, df_sample ) {
  
  df <- df_vars %>%
    inner_join( df_rc, by = c('id_sample', 'chrom', 'pos') ) %>%
    select( caller = name_caller, id_sample, chrom, pos, rc_ref, rc_alt) %>%
    mutate( vaf = if_else(rc_ref+rc_alt==0, 0, (rc_alt)/(rc_ref+rc_alt) )) %>%
    inner_join( df_sample, by = 'id_sample' )
  
  p_rc_alt <- ggplot( df, aes(x = vaf) ) + 
    geom_histogram( aes(y = ..density..), binwidth = 0.02, fill = 'grey' ) +
    geom_density( aes(color = id_sample) ) +
    scale_color_manual( values = c('steelblue', 'red') ) +
    facet_grid( caller~patient, scale = 'free_y' ) +
    theme_dark() +
    theme( strip.text.y = element_text(angle = 0) )
  
}

# VAF spectrum for private vs. shared calls
plot_vaf_bar_empirical_private <- function( df_vars, df_rc, df_pres ) {
  
  df <- df_vars %>%
    inner_join( df_rc, by = c('id_sample', 'chrom', 'pos') ) %>%
    unite( id_mut, id_sample, chrom, pos ) %>%
    inner_join( df_pres, by = c('id_mut') ) %>%
    mutate( private = rowSums(select(., -c(1:7))) == 1 ) %>%
    select( caller = name_caller, id_mut, rc_ref, rc_alt, private ) %>%
    mutate( vaf = if_else(rc_ref+rc_alt==0, 0, (rc_alt)/(rc_ref+rc_alt) ))
  
  p <- ggplot( df, aes(x = vaf) ) + 
    geom_histogram( aes(fill = private), bins = 50, position = 'dodge' ) +
    facet_wrap( ~caller, ncol = 4, scales = 'free_y' ) +
    #ggtitle( 'VAF distribution for classes of variant calls' ) + 
    labs( x = 'variant allele frequency' ) +
    theme_minimal() +
    scale_fill_manual(
      labels = c('shared', 'private'),
      values = c('steelblue', 'red')
    ) +
    theme( legend.position = 'bottom' ) +
    labs( fill = '' )
  
  p
}

plot_vaf_bar_ybreak <- function( df_vars, df_rc, df_rep, lbl_caller, break_at = 5000 ) {
  
  #Function to transform data to y positions
  trans <- function(x){ pmin(x, break_at) + 0.05*pmax(x-break_at,0) }
  yticks <- c(0, 1000, 2000, 3000, 5000, 10000)
  
  df <- df_vars %>% dplyr::filter( name_caller == lbl_caller ) %>% 
    inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
    select(caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt) %>%
    mutate( vaf = (rc_alt)/(rc_ref+rc_alt) ) %>%
    mutate( bin = cut_width(vaf, 0.02) ) %>%
    group_by( caller, type, bin ) %>%
    summarise( n = n() ) %>%
    complete( caller, nesting(type, bin), fill = list(n=0) ) %>%
    mutate( n_t = trans(n) )
  df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )
  
  p_rc_alt <- ggplot(df, aes(x = bin, y = n_t, fill = type)) + 
    geom_histogram( stat = 'identity', position = 'dodge' ) +
    geom_rect( aes(xmin = 0, xmax = 50, ymin = break_at, ymax = break_at+100), fill = 'grey' ) +
    scale_y_continuous(limits = c(0,NA), breaks = trans(yticks), labels = yticks ) +
    #scale_x_discrete(limits = c(1,51), breaks = 1:51, labels = seq(0.0, 1.0, 0.02) ) +
    scale_fill_brewer( palette = 'Set1' ) +
    facet_wrap( ~caller, ncol = 4, scales = 'free_y' ) +
    #ggtitle( 'VAF distribution for classes of variant calls' ) + 
    xlab( 'variant allele frequency' ) +
    theme_pubclean() +
    theme( legend.position = 'bottom' )
}


plot_vaf_dens_srsv <- function( df_vars, df_rc, df_rep ) {

  df <- df_vars %>% dplyr::filter( name_caller != 'Strelka1' ) %>% 
    inner_join( df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
    select(caller = name_caller, id_rep, id_sample, chrom, pos, type, rc_ref, rc_alt) %>%
    mutate( vaf = (rc_alt+1)/(rc_ref+rc_alt+1) ) %>%
    inner_join( df_rep, by = 'id_rep' )
  df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )
  
  p_rc_alt <- ggplot(df, aes(x = vaf)) + 
    geom_density(aes(fill = type), alpha = 0.3) + xlim(0.0, 1.0) +
    scale_fill_brewer(palette = 'Set1' ) +
    facet_wrap( ~ caller, ncol = 4, scales = 'free_y' ) +
    #ggtitle('ALT allele read count distribution for classes of variant calls') + 
    xlab( 'variant allele frequency' ) +
    theme_pubclean() +
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
    mutate( vaf = (rc_alt+1)/(rc_ref+rc_alt+1) ) %>%
    inner_join( df_rep, by = 'id_rep' )
  df$type <- factor( df$type, levels = c('FP', 'FN', 'TP') )
  
  p_rc_alt <- ggplot(df, aes(x = vaf, y = type)) + 
    geom_density_ridges(aes(fill = type), alpha = 0.3) + 
    coord_cartesian( xlim = c(0.0, 1.0) ) +
    scale_fill_brewer(palette = 'Set1' ) +
    facet_wrap( ~ caller, ncol = 4 ) +
    #ggtitle('ALT allele read count distribution for classes of variant calls') + 
    xlab( 'variant allele frequency' ) +
    theme_pubclean() +
    theme(
      strip.text.y = element_text(angle=0),
      legend.position = 'bottom'
    )
}

plot_rc_alt_ridges_srsv2 <- function( df_vars, df_rc, df_rep ) {
  
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
    facet_wrap( ~caller, scales = 'free', ncol = 4 ) +
    #ggtitle('VAF distribution for classes of variant calls') + 
    xlab( 'variant allele frequency' ) +
    theme_pubclean() +
    theme(
      strip.text.y = element_text(angle=0),
      legend.position = 'bottom'
    )
}