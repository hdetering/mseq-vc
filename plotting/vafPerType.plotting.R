plotVAFperType <- function(results,readcounts, callers){

p_temp= results %>% 
  left_join(readcounts %>% select(id_rep, id_sample, chrom, pos, rc_ref, rc_alt)) %>%     # When joining by ref and alt too, the FPs don't appear
  inner_join(callers) %>%
  mutate(type = factor(type, levels = c("FP", "TP", "FN"))) %>%
  mutate(VAF = if_else(rc_ref+rc_alt ==0, 0,rc_alt/(rc_alt+rc_ref) ))
p= ggplot(p_temp)+
  geom_histogram(aes(x = VAF, fill = type), position = "dodge")+
  facet_wrap(~name_caller, ncol = 4, scales = "free_y")+
   ylab("Number of mutations")+
   theme(strip.text.x = element_text(size = 14))

return(p)

}
