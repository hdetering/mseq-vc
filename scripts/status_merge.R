#!/usr/bin/env Rscript
require(tidyverse, quietly=T)
require(yaml, quietly=T)

# load meta info about replicate sim runs
reps <- list.dirs(path='sims', recursive=F)
df <- sapply(reps, FUN=function(p) {
  rep <- basename(p)
  meta <- read_yaml( file.path(p, 'meta.yml') )
  conf <- read_yaml( file.path(p, 'config.yml') )
  c(rep=rep, nclones=meta$nclones, nsamples=meta$nsamples, ttype=meta$ttype, depth=conf$'seq-coverage')
}) %>% t %>% as.tibble

df.bwa      <- read_tsv('status_bwa.tsv', col_names=c('idx', 'rep', 'bwa_mem'))
df.mutect1  <- read_tsv('status_mutect1.tsv', col_names=c('idx', 'rep', 'mutect1'))
df.mutect2  <- read_tsv('status_mutect2.tsv', col_names=c('idx', 'rep', 'mutect2'))
df.vardict  <- read_tsv('status_vardict.tsv', col_names=c('idx', 'rep', 'vardict'))
df.multisnv <- read_tsv('status_multisnv.tsv', col_names=c('idx', 'rep', 'multisnv'))

df <- df %>% 
      inner_join(df.bwa %>% select(rep, bwa_mem), by='rep') %>%
      inner_join(df.mutect1 %>% select(rep, mutect1), by='rep') %>%
      inner_join(df.mutect2 %>% select(rep, mutect2), by='rep') %>%
      inner_join(df.vardict %>% select(rep, vardict), by='rep') %>%
      inner_join(df.multisnv %>% select(rep, multisnv), by='rep')

writeLines(text=format_delim(df, delim='\t'), con=stdout())
