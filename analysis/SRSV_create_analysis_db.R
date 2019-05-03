require(tidyverse)
require(dbplyr)
require(RSQLite)
require(yaml)
require(vcfR)

source('SRSV_create_analysis_db.util.R')

# connection to database
db <- 'data/SRSV/analysis.db'
con <- DBI::dbConnect(RSQLite::SQLite(), db)

# path to simulations output
sims_root <- "/home/uvi/be/hde/lustre/m-seq_varcall/sims"
sims_reps <- dir( sims_root, include.dirs = F )
# path to variant calling results
data_root <- "/home/uvi/be/hde/lustre/m-seq_varcall/data"
data_reps <- dir( data_root, include.dirs = F )

#------------------------------------------------------------------------------
# REPLICATES
#------------------------------------------------------------------------------

if ( !DBI::dbExistsTable( con, 'replicates' ) ) {
  df_rep <- tibble(
    id_rep   = integer(),
    name_rep = character(),
    nclones  = numeric(), 
    nsamples = numeric(),
    ttype    = character(),
    cvg      = numeric()
  )
  n_reps <- 0
  # loop over replicate directories
  for ( rep in sims_reps ) {
    n_reps <- n_reps + 1
    # read meta info about replicate
    fn_meta <- file.path( sims_root, rep, 'meta.yml' )
    meta <- read_yaml( fn_meta )
    # read config file for replicate
    fn_conf <- file.path( sims_root, rep, 'config.yml' )
    conf <- read_yaml( fn_conf )
    
    # marshal relevant info for replicate
    df_rep <- df_rep %>% rbind( tibble(
      id_rep   = as.integer( n_reps ),
      name_rep = rep, 
      nclones  = meta$nclones,
      nsamples = meta$nsamples, 
      ttype    = meta$ttype,
      cvg      = conf$`seq-coverage` ) )
  }
  # write data to table
  copy_to( con, df_rep, 'replicates',
           temporary = FALSE,
           index = list( c( 'id_rep' ) )
  )
} # if table doesn't exist

# get reference to replicates table
tbl_rep <- tbl( con, 'replicates' )
df_rep <- tbl_rep %>% collect()

#------------------------------------------------------------------------------
# SOMATIC MUTATIONS (truth)
#------------------------------------------------------------------------------

if ( !DBI::dbExistsTable( con, 'mutations' ) ) {
  df_mut <- tibble(
    id_rep = integer(),
    id_mut = character(),
    chrom  = character(),
    pos    = integer(),
    ref    = character(),
    alt    = character()
  )
  # loop over replicate directories
  for ( rep in sims_reps ) {
    # read mutated positions
    fn_mut_som <- file.path( sims_root, rep, 'output', 'somatic.vcf' )
    vcf <- read.vcfR( fn_mut_som ) %>% vcfR2tidy()
    mut_som <- vcf$fix %>%
      mutate( name_rep = rep ) %>% 
      inner_join( df_rep, by = c('name_rep') ) %>%
      select( id_rep, id_mut = ID, chrom = CHROM, pos = POS, ref = REF, alt = ALT )
    df_mut <- df_mut %>% rbind( mut_som )
  }
  # write data to table
  copy_to( con, df_mut, 'mutations',
           temporary = FALSE,
           index = list( c( 'id_rep', 'id_mut' ) )
  )
} # if table doesn't exist

#------------------------------------------------------------------------------
# GERMLINE MUTATIONS (truth)
#------------------------------------------------------------------------------

if ( !DBI::dbExistsTable( con, 'snps' ) ) {
  df_snp <- tibble(
    id_rep = integer(),
    id_mut = character(),
    chrom  = character(),
    pos    = integer(),
    ref    = character(),
    alt    = character(),
    is_hom = logical()
  )
  # loop over replicate directories
  for ( rep in sims_reps ) {
    # read mutated positions
    fn_mut_som <- file.path( sims_root, rep, 'output', 'germline.vcf' )
    vcf <- read.vcfR( fn_mut_som ) %>% vcfR2tidy()
    mut_gl <- vcf$fix %>%
      inner_join(vcf$gt, by=c('ChromKey', 'POS')) %>%
      mutate( name_rep = rep, is_hom = (gt_GT=='1/1') ) %>% 
      inner_join( df_rep, by = c('name_rep') ) %>%
      select( id_rep, id_mut = ID, chrom = CHROM, pos = POS, ref = REF, alt = ALT, is_hom )
    df_snp <- df_snp %>% rbind( mut_gl )
  }
  # write data to table
  copy_to( con, df_snp, 'snps',
           temporary = FALSE,
           index = list( c( 'id_rep', 'id_mut' ) )
  )
} # if table doesn't exist

# get reference to replicates table
tbl_snp <- tbl( con, 'snps' )
df_snp <- tbl_snp %>% collect()

#------------------------------------------------------------------------------
# CLONE PREVALENCES
#------------------------------------------------------------------------------

if ( !DBI::dbExistsTable( con, 'clones_prev' ) ) {
  df_prev <- tibble(
    id_rep    = integer(),
    id_sample = character(),
    id_clone  = character(),
    prev      = numeric()
  )
  # loop over replicate directories
  for ( rep in sims_reps ) {
    # read prevalence matrix
    fn_prev <- file.path( sims_root, rep, 'clone_prev.csv' )
    prev <- read.table( fn_prev, sep = ',', header = T, row.names = 1 )
    # NOTE: clone order mismatch between prevalence and mut matrices!
    #prev <- prev %>% rename_all( funs(c(mm.id_clone)) )
    prev <- prev %>% rownames_to_column( var = 'id_sample' ) %>% 
      mutate( name_rep = rep )
    prev <- prev %>% gather( id_clone, prev, -name_rep, -id_sample ) %>%
      inner_join( df_rep, by = c('name_rep') ) %>%
      select( id_rep, id_sample, id_clone, prev )
    df_prev <- df_prev %>% rbind( prev )
  }
  # write data to table
  copy_to( con, df_prev, 'clones_prev',
           temporary = FALSE,
           index = list( c( 'id_rep', 'id_sample', 'id_clone' ) )
  )
} # if table doesn't exist

# get reference to replicates table
tbl_prev <- tbl( con, 'clones_prev' )
df_prev <- tbl_prev %>% collect()

#------------------------------------------------------------------------------
# CLONES - MUTATION mapping
#------------------------------------------------------------------------------

if ( !DBI::dbExistsTable( con, 'clones_mutations' ) ) {
  df_mut_clone <- tibble(
    id_rep     = integer(),
    id_mut     = character(),
    id_clone   = character(),
    is_present = logical()
  )
  # loop over replicate directories
  for ( rep in sims_reps ) {
    # read mutation matrix
    fn_mm <- file.path( sims_root, rep, 'output', 'mm.csv' )
    mm <- read.table( fn_mm, sep = ',', header = F, as.is = T,
                      col.names = c( 'id_clone', paste0( 's', 0:99 ) ) )
    mm.id_clone <- mm %>% subset( id_clone != 'H' ) %>% pull( id_clone ) %>% as.vector()
    mm <- mm %>% mutate( name_rep = rep )
    mm <- mm %>% gather( id_mut, is_present, -name_rep, -id_clone ) %>%
      inner_join( df_rep, by = c('name_rep') ) %>%
      select( id_rep, id_mut, id_clone, is_present )
    df_mut_clone <- df_mut_clone %>% bind_rows( mm )
  }
  # write data to table
  copy_to( con, df_mut_clone, 'clones_mutations',
           temporary = FALSE,
           index = list( c( 'id_rep', 'id_mut', 'id_clone' ) )
  )
} # if table doesn't exist

# get reference to replicates table
tbl_mut_clone <- tbl( con, 'clones_mutations' )
df_mut_clone <- tbl_mut_clone %>% collect()

#------------------------------------------------------------------------------
# MUTATIONS - PRESENCE / ABSENCE
#------------------------------------------------------------------------------

if ( !DBI::dbExistsTable( con, 'mutations_samples' ) ) {
  # calculate expected VAFs for somatic variants in samples
  df_mut_sample <- df_mut %>% 
    inner_join( df_mut_clone, by = c('id_rep', 'id_mut') ) %>% 
    inner_join( df_prev, by = c('id_rep', 'id_clone') ) %>% 
    select( id_rep, id_sample, id_mut, id_clone, prev, is_present ) %>%
    group_by( id_rep, id_sample, id_mut ) %>%
    summarise( vaf_exp = 0.5*sum(is_present*prev) ) %>%
    ungroup() %>%
    mutate( is_present = (vaf_exp > 0.0) )
  # get read counts from simulation output
  fn_csv <- list.files(sims_root, pattern = 'R.\\.vars\\.csv', recursive = T)
  pattern <- '([^/]+)/output/bam/([^\\.]+)\\.vars\\.csv'
  df_rc <- tibble(fp = fn_csv) %>% 
    mutate( csv = map( fp, ~ read_tsv( 
      file.path( sims_root, . ), 
      col_names = c('id_mut', 'rc_tot', 'rc_alt'),
      col_types = cols(
        id_mut = col_character(),
        rc_tot = col_integer(),
        rc_alt = col_integer()
      ) ) ) ) %>% 
    extract( fp, c('name_rep', 'id_sample'), pattern) %>% 
    inner_join( df_rep, by = c('name_rep') ) %>%
    select( id_rep, id_sample, csv) %>%
    unnest() %>%
    dplyr::filter( str_detect( id_mut, '^s' ) )
  # add read count info
  df_mut_sample <- df_mut_sample %>%
    inner_join( df_rc, by = c('id_rep', 'id_sample', 'id_mut') )
  
  # write data to table
  copy_to( con, df_mut_sample, 'mutations_samples',
           temporary = FALSE,
           index = list( c( 'id_rep', 'id_sample', 'id_mut' ) )
  )
} # if table doesn't exist

# get reference to table
tbl_mut_sample <- tbl( con, 'mutations_samples' )
df_mut_sample <- tbl_mut_sample %>% collect()

#===============================================================================

#------------------------------------------------------------------------------
# VARIANT CALLS
#------------------------------------------------------------------------------

# specify for which tools to import data (will delete existing records first)
do_import_for <- c(
  # 'HaplotypeCaller',
  # 'Bcftools',
  # 'CaVEMan',
  # 'MuTect1',
  # 'Mutect2',
  # 'NeuSomatic',
  # 'Shimmer',
  # 'SNooPer',
  # 'SomaticSniper',
  # 'Strelka1',
  # 'Strelka2',
  # 'VarDict',
  # 'VarScan',
  # 'MultiSNV',
  # 'SNV-PPILP',
  # 'MuClone',
  'Mutect2_mseq'
)

# loop over callers  
for ( caller in do_import_for ) {
  
  df_vars_caller <- tibble(
    id_caller = character(),
    id_rep    = character(),
    id_sample = character(),
    chrom     = character(),
    pos       = integer()
  )
  df_vars_caller_mseq <- tibble(
    id_caller = character(),
    id_rep    = character(),
    id_sample = character(),
    chrom     = character(),
    pos       = integer()
  )
  
  # loop over data replicates
  for ( rep in data_reps ) {
    
    # set filename for caller
    fn_vcf <- file.path( data_root, rep, 
                         switch( caller,
                                 HaplotypeCaller = 'haplotypecaller.filt.vcf',
                                 Bcftools = 'bcftools.vcf',
                                 CaVEMan = 'caveman.vcf',
                                 MuTect1 = 'mutect1.vcf',
                                 Mutect2 = 'mutect2.vcf',
                                 NeuSomatic = 'neusomatic.vcf',
                                 SNooPer = 'snooper.fix.vcf',
                                 SomaticSniper = 'somaticsniper.vcf',
                                 Strelka1 = 'strelka1.fix.vcf',
                                 Strelka2 = 'strelka2.fix.vcf',
                                 Shimmer = 'shimmer.vcf',
                                 VarDict = 'vardict.vcf',
                                 VarScan = 'varscan.vcf',
                                 MultiSNV = 'multisnv.vcf',
                                 `SNV-PPILP` = 'snv-ppilp.vcf',
                                 MuClone = 'muclone.vcf',
                                 Mutect2_mseq = 'mutect2m.vcf'
                         ) )
    
    
    # get variant calls
    if ( file.exists( fn_vcf ) ) {
      print( fn_vcf )
      
      vars <- getVariantCalls( fn_vcf, caller ) %>%
        mutate( name_caller = caller, name_rep = rep ) %>%
        select( name_caller, id_rep, id_sample, chrom = CHROM, pos = POS )
      df_vars_caller <- df_vars_caller %>% rbind( vars )
      
      vars_mseq <- getVariantCalls( fn_vcf, caller, mseq = TRUE ) %>%
        mutate( name_caller = caller, name_rep = rep ) %>%
        select( name_caller, id_rep, id_sample, chrom = CHROM, pos = POS )
      df_vars_caller_mseq <- df_vars_caller_mseq %>% rbind( vars_mseq )
    }
    
  } # replicates
  
  write_csv( df_vars_caller, sprintf( 'df_vars.%s.csv', tolower(caller) ) )
  write_csv( df_vars_caller_mseq, sprintf( 'df_vars_mseq.%s.csv', tolower(caller) ) )
  
} # caller



# write all new data to table
#dbAppendTable( con, 'varcalls', df_vars )
# NOTE: this is way faster than through DBI
callers <- c(
  'HaplotypeCaller',
  'Bcftools',
  'CaVEMan',
  'MuTect1',
  'Mutect2',
  'NeuSomatic',
  'Shimmer',
  'SNooPer',
  'SomaticSniper',
  'Strelka1',
  'Strelka2',
  'VarDict',
  'VarScan',
  'MultiSNV',
  'SNV-PPILP',
  'MuClone',
  'Mutect2_mseq'
)

fn_callers <- sprintf( 'df_vars.%s.csv', tolower( callers ) )
fn_csv <- fn_callers[ file.exists( fn_callers ) ]
df_vars <- lapply( fn_csv, read_csv ) %>% bind_rows()
write_csv(df_vars, 'df_vars.csv')
df_vars <- df_vars %>% 
  inner_join( df_rep, by = c('name_rep') ) %>%
  select( name_caller, name_rep, id_sample, chrom, pos )

fn_callers_m <- sprintf( 'df_vars_mseq.%s.csv', tolower( callers ) )
fn_csv_m <- fn_callers[ file.exists( fn_callers_m ) ]
df_vars_m <- lapply( fn_csv_m, read_csv ) %>% bind_rows()
write_csv(df_vars_m, 'df_vars_mseq.csv')
df_vars_m <- df_vars_m %>% 
  inner_join( df_rep, by = c('name_rep') ) %>%
  select( name_caller, name_rep, id_sample, chrom, pos )

#------------------------------------------------------------------------------
# VARIANT CALLERS
#------------------------------------------------------------------------------

# populate caller table
df_callers <- df_vars %>% group_by( name_caller ) %>%
  summarise( n = n() ) %>%
  rowid_to_column( 'id_caller' ) %>%
  select( id_caller, name_caller )
if ( DBI::dbExistsTable( con, 'callers' ) ) {
  dbExecute( con, "DROP TABLE callers;" )
}
copy_to( con, df_callers, 'callers',
         temporary = FALSE,
         index = list( c('id_caller') )
)

#------------------------------------------------------------------------------
# VARIANT CALLS
#------------------------------------------------------------------------------

df_varcall <- df_vars %>% 
  inner_join( df_rep, by = c('name_rep') ) %>%
  inner_join( df_callers, by = c('name_caller') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos )

# delete existing records
dbExecute( con, "DROP TABLE varcalls;" )
# create table anew
copy_to( con, df_varcall, 'varcalls',
         temporary = FALSE,
         index = list( c( 'id_caller', 'id_rep', 'id_sample' ) )
)

# get reference to varcalls table
tbl_varcall <- tbl( con, 'varcalls' )
df_varcall <- tbl_varcall %>% collect()

# M-seq filter
#-------------------------------------------------------------------------------

df_varcall_m <- df_vars_m %>% 
  inner_join( df_rep, by = c('name_rep') ) %>%
  inner_join( df_callers, by = c('name_caller') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos )

# delete existing records
dbExecute( con, "DROP TABLE varcalls_mseq;" )
# create table anew
copy_to( con, df_varcall_m, 'varcalls_mseq',
         temporary = FALSE,
         index = list( c( 'id_caller', 'id_rep', 'id_sample' ) )
)

# get reference to varcalls table
tbl_varcall_m <- tbl( con, 'varcalls_mseq' )
df_varcall_m <- tbl_varcall_m %>% collect()

#------------------------------------------------------------------------------
# READ COUNTS
#------------------------------------------------------------------------------

fn_reps <- file.path( data_root, data_reps )
fn_counts <- as.vector(sapply(fn_reps, function(x) file.path(x, sprintf( "R%s.allelicCounts.tsv", c(1:5, 'N') ))))
fn_tsv <- fn_counts[ file.exists( fn_counts ) ]
#df_rc <- lapply( fn_tsv, read_ac_tsv ) %>% bind_rows()
pattern <- '.*/([^/]+)/([^\\.]+)\\.allelicCounts\\.tsv'
df_rc <- tibble(fp = fn_tsv) %>% 
  mutate( tsv = map( fp, ~ read_tsv(
    file = .,
    comment = '@',
    col_types = cols(
      CONTIG = col_character(),
      POSITION = col_integer(),
      REF_COUNT = col_integer(),
      ALT_COUNT = col_integer(),
      REF_NUCLEOTIDE = col_character(),
      ALT_NUCLEOTIDE = col_character()
    ) ) ) ) %>% 
  extract( fp, c('name_rep', 'id_sample'), pattern) %>% 
  inner_join( df_rep, by = c('name_rep') ) %>%
  unnest() %>%
  select( id_rep, id_sample, chrom = CONTIG, pos = POSITION, 
          ref = REF_NUCLEOTIDE, alt = ALT_NUCLEOTIDE, 
          rc_ref = REF_COUNT, rc_alt = ALT_COUNT)

write_csv(df_vars, 'df_rc.csv')

# delete existing records
dbExecute( con, "DROP TABLE readcounts;" )
# create table anew
copy_to( con, df_rc, 'readcounts',
         temporary = FALSE,
         index = list( c( 'id_rep', 'id_sample' ) )
)

# get reference to varcalls table
tbl_rc <- tbl( con, 'readcounts' )
df_rc <- tbl_rc %>% collect()


#------------------------------------------------------------------------------
# PERFORMANCE STATS
#------------------------------------------------------------------------------

# join detected variants and true mutations
df_tp <- df_varcall %>%
  inner_join( df_mut, by = c('id_rep', 'chrom', 'pos') ) %>%
  inner_join( df_mut_sample %>% dplyr::filter( as.logical( is_present ) ), 
              by = c('id_rep', 'id_sample', 'id_mut') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos ) %>%
  mutate( type = 'TP') 
df_fp <- df_varcall %>% 
  anti_join( df_mut_sample %>% dplyr::filter( as.logical( is_present ) ) %>%
               inner_join( df_mut, by = c('id_rep', 'id_mut') ), 
             by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos ) %>%
  mutate(type = 'FP')
df_fn <- df_mut_sample %>% dplyr::filter( as.logical( is_present ) ) %>%
  inner_join( df_mut, by = c('id_rep', 'id_mut') ) %>%
  crossing( tibble( id_caller = seq( length(callers) ) ) ) %>%
  anti_join( df_varcall, by = c('id_rep', 'id_caller', 'id_sample', 'chrom', 'pos') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos ) %>%
  mutate( type = 'FN' )

df_vars <- df_tp %>% rbind( df_fp ) %>% rbind( df_fn )
df_eval <- df_vars %>% select( id_caller, id_rep, type ) %>%
  group_by( id_caller, id_rep, type ) %>%
  summarise( n = n() ) %>%
  ungroup() %>%
  complete( id_caller, id_rep, type, fill = list(n = 0) ) %>%
  spread( type, n ) %>%
  mutate( recall = TP/(TP+FN), precision = TP/(TP+FP) ) %>%
  mutate( F1 = 2*recall*precision/(recall+precision) )

df <- df_eval %>%
  inner_join( df_callers, by = c('id_caller') ) %>% 
  inner_join( df_rep, by = c('id_rep') )
# write all variant calls to file
saveRDS( df_vars, 'df_vars.rds' )
# write summary stats to file
saveRDS( df, 'df_eval.sample.rds' )

# M-seq filter stats
#-------------------------------------------------------------------------------

# join detected variants and true mutations
df_tp_m <- df_varcall_m %>%
  inner_join( df_mut, by = c('id_rep', 'chrom', 'pos') ) %>%
  inner_join( df_mut_sample %>% dplyr::filter( as.logical( is_present ) ), 
              by = c('id_rep', 'id_sample', 'id_mut') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos ) %>%
  mutate( type = 'TP') 
df_fp_m <- df_varcall_m %>% 
  anti_join( df_mut_sample %>% dplyr::filter( as.logical( is_present ) ) %>%
               inner_join( df_mut, by = c('id_rep', 'id_mut') ), 
             by = c('id_rep', 'id_sample', 'chrom', 'pos') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos ) %>%
  mutate(type = 'FP')
df_fn_m <- df_mut_sample %>% dplyr::filter( as.logical( is_present ) ) %>%
  inner_join( df_mut, by = c('id_rep', 'id_mut') ) %>%
  crossing( tibble( id_caller = seq( length(callers) ) ) ) %>%
  anti_join( df_varcall, by = c('id_rep', 'id_caller', 'id_sample', 'chrom', 'pos') ) %>%
  select( id_caller, id_rep, id_sample, chrom, pos ) %>%
  mutate( type = 'FN' )

df_vars_m <- df_tp_m %>% rbind( df_fp_m ) %>% rbind( df_fn_m )
df_eval_m <- df_vars_m %>% select( id_caller, id_rep, type ) %>%
  group_by( id_caller, id_rep, type ) %>%
  summarise( n = n() ) %>%
  ungroup() %>%
  complete( id_caller, id_rep, type, fill = list(n = 0) ) %>%
  spread( type, n ) %>%
  mutate( recall = TP/(TP+FN), precision = TP/(TP+FP) ) %>%
  mutate( F1 = 2*recall*precision/(recall+precision) )

df <- df_eval_m %>%
  inner_join( df_callers, by = c('id_caller') ) %>% 
  inner_join( df_rep, by = c('id_rep') )
# write all variant calls to file
saveRDS( df_vars_m, 'df_vars_mseq.rds' )
# write summary stats to file
saveRDS( df, 'df_eval_mseq.sample.rds' )
