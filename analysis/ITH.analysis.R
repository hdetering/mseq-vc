calculate_fst <- function( df_varcall , 
                           df_af ){
  
  t = df_af %>% 
    gather(key = "id_sample", 
           value = "af",
           T1:T5) %>% 
    dplyr::rename(name_caller = software,
                  name_rep = replicate) %>%
    separate(mut_info, into = c("chrom", "pos", "ref", "alt"), sep ="_") %>%
    replace_na(list(af = 0)) %>% # 6659195
    mutate(name_caller = str_replace(name_caller, "Corrected", "")) %>% 
    left_join(df_rep) %>%
    left_join(df_caller) %>% 
    inner_join(df_varcall %>% mutate(pos = as.character(pos)))  # Should be right_join?
  
  
  t2 = t %>% select(chrom_pos, ref, alt, name_caller, name_rep, id_sample, af)
  
  sample_combinations = t2  %>% mutate(id_sample_2 = id_sample) %>% expand(id_sample, id_sample_2) %>% filter(id_sample_2>id_sample)
  
  t3 = t2 %>% mutate(id_sample_2 = id_sample, af_2 = af) %>% full_join(sample_combinations)
  
  t3 = t2 %>% select(-c(ref,alt)) 
  t4 = t2 %>% dplyr::rename(id_sample_2 = id_sample, af_2 = af) %>% select(-c(ref,alt)) 
  
  
  t5 = inner_join(t3,t4, by = c("chrom_pos", "name_caller", "name_rep")) %>% 
    # mutate(id_sample = as.factor(id_sample),
    #        id_sample_2 = as.factor(id_sample_2)) %>%
    tidyr::complete(id_sample, id_sample_2, fill = list(af = 0, af_2 = 0)) %>%    # Missing combination?
    filter(id_sample_2 > id_sample)
  # 
  # inner_join(t3 %>% 
  #              filter(name_caller =="HaplotypeCaller"),
  #            t4 %>% 
  #              filter(name_caller =="HaplotypeCaller"), 
  #            by = c("chrom_pos", "name_caller", "name_rep")) %>%
  #   complete(id_sample, id_sample_2) %>% View()
  # 
  
  # af_pairwise =  data.frame(chrom_pos = character(),
  #                           name_calle = character(),
  #                           name_rep = character(),
  #                           id_sample = character(), 
  #                           af = character(),
  #                           pair = character(),
  #                           stringsAsFactors = FALSE)    
  # 
  # ids_samples = unique(sort(t2$id_sample))
  # for (i in 1:(length(ids_samples)-1)){
  #   for(j in (i+1):length(ids_samples)){
  #     id_sample_1 = ids_samples[i]
  #     id_sample_2 = ids_samples[j]
  #     
  #     pair = t2 %>% unique() %>%
  #       dplyr::filter(id_sample == id_sample_1 | id_sample == id_sample_2) %>% 
  #       mutate(id_sample = case_when(id_sample ==id_sample_1 ~ "A",
  #                                    id_sample ==id_sample_2 ~ "B")) %>%
  #       spread(key = id_sample, value = af, fill = 0) %>%
  #       mutate(pair = paste(id_sample_1, id_sample_2, sep = "_"))
  #       
  #     af_pairwise = rbind(af_pairwise, pair)  
  #   }
  #   
  # }
  
  t5 %>% 
    tidyr::separate(af, into = c("af"), extra = "drop", sep = ",") %>%
    tidyr::separate(af_2, into = c("af_2"), extra = "drop", sep = ",") %>%
    mutate(af = as.numeric(af),
           af_2 = as.numeric(af_2)) %>% 
    # filter(!(af == 0 & af_2 == 0)) %>%
    mutate(het_1 = 2*af*(1-af),
           het_2 = 2*af_2*(1-af_2),
           ht = 2*(((1-af) + (1-af_2)) / 2)*((af + af_2) / 2)) %>%  
    mutate(hs = (het_1+het_2)/2) %>%
    mutate(fst = if_else(ht == 0, NaN, (ht-hs)/ht)) %>% 
    # mutate(fst = if_else(ht == 0, 0, (ht-hs)/ht)) %>%   # Could I do this if I remove the af=0 and af_2=0? or is it better to keep NA and don't take that into account for the mean? Maybe yes, since when comparing two populations you probably don't use the shared variants (how could you even tell there is a variant there without a ref genome?)
    group_by(id_sample, id_sample_2, name_caller, name_rep) %>% # Mean per pair and then per tumor?
    dplyr::summarise(fst = mean(fst, na.rm = TRUE)) %>%            # Remove NA values? Why do they arise? When both are 0 or both are 1
    group_by(name_caller, name_rep) %>%
    dplyr::summarise(fst = mean(fst))
  
  
  
}