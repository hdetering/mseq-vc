plotRecurrenceFPs <- function(results, callers, replicates){

library(ggpubr)

pairwise_intersections = NA
for(caller in callers$name_caller){
  
  temp = results %>% 
  filter(type =="FP") %>%
  inner_join(callers) %>%
  filter(name_caller ==caller) %>%
  inner_join(replicates)

temp2 = temp %>%
  select(name_rep, id_sample,pos) %>%
  mutate(tempcol = 1) %>%
  mutate(tempcol = as.integer(as.character(tempcol))) %>%
  unite(sample, id_sample, name_rep, sep = ".") %>%
  spread(key = sample, value = tempcol, fill = as.integer(0) ) %>%
  mutate(pos = as.character(pos))

rownames(temp2) = temp2$pos
temp2 = temp2 %>% select(-c(pos))


for (i in 1:ncol(temp2)){
  for (j in i:ncol(temp2)){
    sample1 = colnames(temp2)[i]
    sample2 = colnames(temp2)[j]
  
    t = temp2[,sample1]+ temp2[,sample2]
    interse = length(t[which(t==2)])
    
    if(class(pairwise_intersections)=="data.frame"){
      pairwise_intersections = rbind(pairwise_intersections, data.frame(sample1 = c(sample1), sample2 = c(sample2), intersection_size = c(interse), caller = caller))
    }else{
      pairwise_intersections = data.frame(sample1 = sample1, sample2 = sample2, intersection_size = interse, caller = caller)
    }
    
    
    
  }}}



p=pairwise_intersections %>%
  separate(sample1, into = c("region1", "scenario1", "replicate1")) %>%
  separate(sample2, into = c("region2", "scenario2", "replicate2")) %>%
  mutate(Case = if_else(region1==region2 & scenario1 == scenario2 & replicate1 == replicate2, 
                        "Same sample",
                        if_else(region1 ==region2,
                                "Same region", 
                                "Different regions"))) %>% 
  filter(Case!="Same sample") %>% 
  ggplot()+
  geom_boxplot(aes(x = Case,y = intersection_size,  color = caller), outlier.size = 0.05) +
  facet_wrap(~caller, scales = "free_y")+
  stat_compare_means(aes(x = Case,y = intersection_size*1.1, size = 12) )+
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "None")+
  ylab("FPs intersection size")+
  xlab("Comparison type")


return(p)




}
