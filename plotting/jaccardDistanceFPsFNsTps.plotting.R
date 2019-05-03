plotJaccardDistance <- function(results, callers){



jaccard_fps = NA
jaccard_fns = NA
jaccard_tps = NA


callernames = callers$name_caller
results2 = results %>% inner_join(callers)


for (i in 1:length(callernames)){
  for (j in i:length(callernames)){
    print(paste(callernames[i], callernames[j]))
    caller1=callernames[i]
    caller2=callernames[j]
    
    caller1matrixt = results2 %>% filter(name_caller==caller1 ) %>% select(c(chrom, pos, type, id_rep, id_sample))
    caller2matrixt = results2 %>% filter(name_caller==caller2)  %>% select(c(chrom, pos, type, id_rep, id_sample))    
    
    jaccardvector = c()
    for (tumor in 1:30){
      for (sample in samples){
        caller1matrix= caller1matrixt %>% filter(id_rep ==tumor & id_sample ==sample & type=="FP") %>% select(c(chrom, pos))
        caller2matrix = caller2matrixt %>% filter(id_rep ==tumor & id_sample ==sample & type=="FP") %>% select(c(chrom, pos))
        jaccard = 2*nrow(intersect(caller1matrix, caller2matrix))/nrow(union_all(caller1matrix, caller2matrix))
        jaccardvector = c(jaccardvector, jaccard)
      }
    }
    
    temp_row = data.frame(caller1 = caller1, caller2 = caller2, jaccard = mean(jaccardvector, na.rm = TRUE))
    if(class(jaccard_fps)!="data.frame"){
      jaccard_fps = temp_row
    }else{
      jaccard_fps = rbind(jaccard_fps, temp_row)
    }
    
    jaccardvector = c()
    for (tumor in 1:30){
      for (sample in samples){
        caller1matrix = caller1matrixt %>% filter(id_rep ==tumor & id_sample ==sample) %>% filter(type=="FN") %>% select(c(chrom, pos))
        caller2matrix = caller2matrixt %>% filter(id_rep ==tumor & id_sample ==sample) %>% filter(type=="FN") %>% select(c(chrom, pos))
        jaccard = 2*nrow(intersect(caller1matrix, caller2matrix))/nrow(union_all(caller1matrix, caller2matrix))
        jaccardvector = c(jaccardvector, jaccard)
      }
    }
    
    temp_row = data.frame(caller1 = caller1, caller2 = caller2, jaccard = mean(jaccardvector, na.rm = TRUE))
    if(class(jaccard_fns)!="data.frame"){
      jaccard_fns = temp_row
    }else{
      jaccard_fns = rbind(jaccard_fns, temp_row)
    }
    
    jaccardvector = c()
    for (tumor in 1:30){
      for (sample in samples){
        caller1matrix = caller1matrixt %>% filter(id_rep ==tumor & id_sample ==sample) %>% filter(type=="TP") %>% select(c(chrom, pos))
        caller2matrix = caller2matrixt %>% filter(id_rep ==tumor & id_sample ==sample) %>% filter(type=="TP") %>% select(c(chrom, pos))
        jaccard = 2*nrow(intersect(caller1matrix, caller2matrix))/nrow(union_all(caller1matrix, caller2matrix))
        jaccardvector = c(jaccardvector, jaccard)
      }
    }
    
    temp_row = data.frame(caller1 = caller1, caller2 = caller2, jaccard = mean(jaccardvector,na.rm = TRUE))
    if(class(jaccard_tps)!="data.frame"){
      jaccard_tps = temp_row
    }else{
      jaccard_tps = rbind(jaccard_tps, temp_row)
    }
  
  }}









p1 = ggplot(jaccard_fps) +
  geom_tile(aes(y=caller1, x= caller2, fill = jaccard))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "None") +
  xlab("")+
  ylab("")+
  ggtitle("Jaccard index for the FPs")+
  scale_fill_distiller(palette = "Spectral", direction = -1)

p2 = ggplot(jaccard_fns) +
  geom_tile(aes(y=caller1, x = caller2, fill = jaccard))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "None") +
  xlab("")+
  ylab("")+
  ggtitle("Jaccard index for the FNs")+
  scale_fill_distiller(palette = "Spectral", direction = -1)

p3 = ggplot(jaccard_tps) +
  geom_tile(aes(y=caller1, x = caller2, fill = jaccard))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  xlab("")+
  ylab("")+
  ggtitle("Jaccard index for the TPs")+
  scale_fill_distiller(palette = "Spectral", direction = -1)

ggarrange(p1, p2, p3,ncol = 3)




}
