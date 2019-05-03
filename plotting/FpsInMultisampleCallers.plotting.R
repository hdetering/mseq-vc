plot <- function(results){
p=results %>%
  inner_join(callers)%>%
  filter(class=="multi-sample") %>%
  group_by(id_caller, id_rep, pos) %>%
  mutate(TPInOtherRegions = sum(type=="TP"|type=="FN")) %>%
  filter(type=="FP") %>%
  filter(TPInOtherRegions>0) %>%
  group_by(TPInOtherRegions, name_caller) %>%
  summarise(count = n()) %>%
  rbind(data.frame(TPInOtherRegions = 3, name_caller = "MuClone",count = 0 ) %>% group_by(name_caller, TPInOtherRegions)) %>%
  ggplot()+
  geom_bar(aes(x = as.character(TPInOtherRegions), y = count, fill = name_caller), stat = "Identity")+
  facet_wrap(~name_caller, nrow = 1) +
  xlab("Number of regions with the variant as TP or FN")+
  ylab("Number of FPs")+
  theme(legend.position = "None")+
  scale_fill_manual(values = c("#AAA235","#87AC34","#57BE9C","#8696F8"))




}
