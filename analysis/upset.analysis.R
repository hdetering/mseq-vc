get_upset_pres <- function(
  df_vars,
  df_rc
)
{
  df_pres <- df_vars %>%
    inner_join(df_rc, by = c('id_rep', 'id_sample', 'chrom', 'pos')) %>%
    unite(id_mut, c('id_rep', 'id_sample', 'chrom', 'pos')) %>% 
    select(id_mut, type, rc_ref, rc_alt, germline, caller = name_caller) %>%
    mutate(present = ifelse(type == 'FN', 0, 1)) %>%
    spread(caller, present, fill = 0) %>%
    mutate(TRUE_somatic = ifelse(type == 'FP', 0, 1), TRUE_germline = as.numeric(germline)) %>%
    select(-germline)
}