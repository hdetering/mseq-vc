plotRecallByAdmixture <- function(results, callers replicates){

temp = results %>%
  group_by(id_rep, id_sample, id_caller) %>%
  summarise(recall = sum(type=="TP")/(sum(type=="TP")+sum(type=="FN")),
            precision = sum(type=="TP")/(sum(type=="TP")+sum(type=="FP")),
            FPs = sum(type=="FP")) %>%
  group_by(id_caller) %>%
  mutate(F1 = 2*(recall*precision)/(recall+ precision)) %>%
  mutate(median_F1 = median(F1)) %>%
  arrange(desc(median_F1)) %>%
  rowid_to_column("f1order") %>%
  inner_join(callers  ) %>%
  inner_join(replicates ) %>%
  mutate(ttype = factor(ttype, levels = c("low", "medium", "high"))) %>%
  gather(key = "Measure", value = "Result", recall) 

p =  ggplot(data = temp) +
  geom_boxplot(aes(x=reorder(name_caller, f1order), y = Result, fill= ttype)) +
  theme(axis.text.x = element_text(angle = 45))+
  xlab("Variant caller")+
  ylab("Recall")+
  scale_fill_discrete(name="Admixture")+
  geom_point(data = temp  %>%
  mutate(Result = if_else(is.na(Result), 0, Result)) %>%
  group_by(Measure, name_caller, ttype) %>% 
  summarise(medianresult=median(Result)) %>% 
  group_by(Measure, ttype) %>% 
  filter((medianresult==max(medianresult) & Measure=="recall") | (medianresult==min(medianresult) & Measure=="FPs")), aes(x=name_caller, y = medianresult),fill = "gold", shape = 23, size = 3)


return(p)

}
