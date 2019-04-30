plotEffectPON <- function(varcallspon){


varcallspon$chrom = as.character(varcallspon$chrom)
resultspon = NA  
callerspon = c("3","4","special_Mutect2NoPon","special_Mutect1NoPon" )
for (caller in callerspon){
  print(caller)
  varcallscaller = varcallspon %>% filter(id_caller == caller & id_sample =="T1") %>% mutate(called = TRUE)
  resultscaller  = varcallscaller %>% 
    full_join(mutations_samples %>% 
                inner_join(mutations) %>% 
                filter(is_present & id_sample =="T1")) %>%
    mutate(Type = if_else(!is.na(is_present), if_else(!is.na(called), "TP", "FN"), "FP")) %>%
    mutate(id_caller = caller)
  if(class(resultspon)!="data.frame"){
    resultspon = resultscaller
  }else{
    resultspon = rbind(resultspon, resultscaller)
  }
}



p = resultspon %>%
  group_by(id_rep, id_sample, id_caller) %>%
  summarise(recall = sum(Type=="TP")/(sum(Type=="TP")+sum(Type=="FN")),
            precision = sum(Type=="TP")/(sum(Type=="TP")+sum(Type=="FP"))) %>%
  left_join(callers %>% mutate(id_caller = as.character(id_caller))) %>%
  mutate(name_caller = if_else(id_caller =="special_Mutect2NoPon", "Mutect2NoPON", if_else(id_caller =="special_Mutect1NoPon", "Mutect1NoPON", name_caller))) %>%
  mutate(Fscore = 2*(recall*precision)/(recall+ precision)) %>%
  inner_join(replicates ) %>%
  gather(key = "Measure", value = "Result", recall, precision, Fscore) %>%
  select(c(Measure, Result, id_caller, name_caller)) %>%
  mutate(Caller = if_else(id_caller == "special_Mutect2NoPon"| id_caller =="4", "Mutect2", "Mutect1"),
         Approach = if_else(id_caller == "special_Mutect2NoPon"| id_caller =="special_Mutect1NoPon", "Without PON", "With PON")) %>%
  ggplot() +
  geom_boxplot(aes(x=Approach, y = Result, fill= name_caller)) +
  facet_grid(Caller~Measure) +
  theme(axis.text.x = element_text(angle = 45))+
  xlab("Approach")+
  scale_fill_discrete(name="Variant caller")+
  ylab("")

return(p)
}



