require(tidyverse)
require(vcfR)

getVariantCalls <- function( filename, id_caller, mseq = FALSE ) {
  vcf <- read.vcfR( fn_vcf ) %>% vcfR2tidy()
  vars <- vcf %>% filterVcf( id_caller, mseq )
  vars
}

filterVcf <- function( vcf, caller, mseq ) {
  # --- Germline Callers ------------------------------------------------------
  if ( caller %in% c('HaplotypeCaller', 'Bcftools') ) {
    vcf$fix %>% 
      inner_join( vcf$gt, by = c( 'ChromKey', 'POS' ) ) %>%
      dplyr::filter( Indiv != 'RN' & gt_GT != '0/0' ) %>%
      mutate(id_sample = Indiv)
  } 
  # --- Pre-filtered variants -------------------------------------------------
  else if ( caller %in% c('CaVEMan', 'Shimmer', 'SomaticSniper') ) {
    vcf$fix %>%
      filterSetField( mseq )
  }
  # --- Mutect, Strelka -------------------------------------------------------
  else if ( caller %in% c('MuTect1', 'Mutect2', 'Strelka1', 'Strelka2', 'NeuSomatic') ) {
    vcf$fix %>%
      dplyr::filter( FILTER == 'PASS' ) %>%
      filterSetField( mseq )
  }
  # --- SNooPer ---------------------------------------------------------------
  else if ( caller == 'SNooPer' ) {
    vcf$fix %>%
      dplyr::filter( FILTER == 'PASS' ) %>%
      inner_join(vcf$gt, by = c('ChromKey', 'POS')) %>%
      mutate(id_sample = str_sub(Indiv, 1, 2)) %>%
      dplyr::filter( id_sample != 'RN' )
  }
  # --- VarDict ---------------------------------------------------------------
  else if ( caller == 'VarDict' ) {
    vcf$fix %>%
      dplyr::filter( FILTER == 'PASS' & grepl('Somatic$', STATUS) ) %>%
      filterSetField( mseq )
  }
  # --- VarScan ---------------------------------------------------------------
  else if ( caller == 'VarScan' ) {
    vcf$fix %>% 
      dplyr::filter( FILTER == 'PASS' ) %>%
      inner_join( vcf$gt, by=c('ChromKey', 'POS') ) %>%
      dplyr::filter( gt_SS == 2 ) %>%
      mutate( id_sample = str_sub(Indiv, 1, 2) )
  }
  # --- MultiSNV --------------------------------------------------------------
  else if ( caller == 'MultiSNV' ) {
    vcf$fix %>% 
      dplyr::filter( FILTER == 'PASS' ) %>%
      inner_join( vcf$gt, by=c('ChromKey', 'POS') ) %>%
      dplyr::filter( Indiv != 'N' & gt_SS == 2 ) %>%
      mutate( id_sample = str_replace(Indiv, 'T', 'R') )
  }
  # --- SNV-PPILP -------------------------------------------------------------
  else if ( caller == 'SNV-PPILP' ) {
    vcf$fix %>% 
      dplyr::filter( FILTER == 'PASS' ) %>%
      inner_join(vcf$gt, by=c('ChromKey', 'POS')) %>%
      dplyr::filter(Indiv != 'RN' & gt_GT == '0/1') %>%
      mutate( id_sample = Indiv )
  }
  # --- MuClone ---------------------------------------------------------------
  else if ( caller == 'MuClone' ) {
    vcf$fix %>% 
      inner_join(vcf$gt, by=c('ChromKey', 'POS')) %>%
      dplyr::filter( gt_GT == '0/1' ) %>%
      mutate( id_sample = Indiv )
  }
  # --- Mutect2 multi-sample --------------------------------------------------
  else if ( caller == 'Mutect2_multi_F' ) {
    vcf$fix %>% 
      dplyr::filter( FILTER == 'PASS' ) %>%
      inner_join(vcf$gt, by=c('ChromKey', 'POS')) %>%
      dplyr::filter(Indiv != 'RN' & gt_GT == '0/1') %>%
      mutate( id_sample = Indiv ) %>%
      rowwise() %>% 
      mutate( rc_alt = as.numeric(str_split(gt_AD, ',', simplify=T)[2]) ) %>%
      dplyr::filter( rc_alt > 0 ) %>%
      select(CHROM, POS, REF, ALT, id_sample, gt_AD, rc_alt)
  }
  else {
    stop( sprintf( '[ERROR] Unknown caller: "%s"', caller) )
  }
}

# analyze GATK CombineVariants' "set" field
# mseq: if TRUE, accept filtered variants if unfiltered in at least one sample
filterSetField <- function( vcf_data, mseq = FALSE ) {
  if ( mseq ) {
    vcf_data %>%
      mutate(set_smp = ifelse(set=='Intersection', 'variant-variant2-variant3-variant4-variant5', set)) %>%
      separate_rows(set_smp) %>%
      dplyr::filter(set_smp != 'filteredInAll') %>%
      mutate(id_sample = ifelse(set_smp=="variant", "R1", paste0("R", str_sub(set_smp, nchar(set_smp)))))
  }
  else {
    vcf_data %>%
      mutate(set_smp = ifelse(set=='Intersection', 'variant-variant2-variant3-variant4-variant5', set)) %>%
      separate_rows(set_smp) %>%
      dplyr::filter(!startsWith(set_smp, 'filterIn')) %>%
      mutate(id_sample = ifelse(set_smp=="variant", "R1", paste0("R", str_sub(set_smp, nchar(set_smp)))))  
  }
}
