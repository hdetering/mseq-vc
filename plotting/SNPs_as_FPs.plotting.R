plot_SNPs_as_FPs <- function(df_vars){
library(ggrepel)
library(ggplot)
library(tidyverse)
temp = results %>% 
  inner_join(callers ) %>%
  filter(type =="FP") %>% 
    left_join(readcounts %>% select(id_rep, id_sample, chrom, pos, rc_ref, rc_alt)) %>%
      mutate(VAF = if_else(rc_ref+rc_alt ==0, 0,rc_alt/(rc_alt+rc_ref) )) %>%
  left_join(snps %>% select(-c(ref, alt)) %>% mutate(SNP=TRUE)) %>% 
    mutate(SNP = if_else(is.na(SNP), FALSE, SNP)) %>%
  group_by(name_caller) %>% summarise(SNPs = sum(SNP)/n(), Total = sum(SNP)) 
  
set.seed(42)
p3 =ggplot(temp)+
  geom_point(aes ( x = 0.3, y = SNPs, size = Total, color = name_caller), alpha = 0.7) +
  theme_light() + theme(axis.text.x=element_blank(),
                        axis.ticks.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.grid.major.x = element_blank(),
                        legend.position = "None")+  
  xlab("")+
  ylab("Proportion of SNPs in the FPs")+
  scale_size_continuous(name="Count of SNPs in the FPs",range = c(1,20))+
  scale_color_discrete(name = "Caller", guide = FALSE)+
  geom_text_repel(aes(y = SNPs, x = 0.3, label = name_caller), 
                  point.padding = unit(0.5, "lines"),
                  nudge_y = 0, 
                  nudge_x=3,
                  direction = "y", 
                  angle = 0, 
                  hjust = 3, 
                  force = 5)+
  xlim(c(0,1.5))


temp3 = results %>%
  filter(type == "FP") %>%
  left_join(snps %>% select(-c(ref, alt)) %>% mutate(SNP=TRUE)) %>% 
    mutate(SNP = if_else(is.na(SNP), FALSE, SNP)) %>% 
  filter(SNP) %>%
  left_join(bothhealthy_readcounts)  %>%
      mutate(healthyVAF = if_else(rc_ref_healthy+rc_alt_healthy ==0, 0,rc_alt_healthy/(rc_alt_healthy+rc_ref_healthy) )) %>%
      mutate(superhealthyVAF = if_else(rc_ref_superhealthy+rc_alt_superhealthy ==0, 0,rc_alt_superhealthy/(rc_alt_superhealthy+rc_ref_superhealthy))) 


 
p1 = ggplot(temp3)+
  geom_point(aes(x = healthyVAF,y = superhealthyVAF, color =rc_ref_superhealthy + rc_alt_superhealthy), size = 0.1)+
   scale_color_distiller(palette = "Spectral", name = "Depth at 300x")+
   xlab("VAF in the 50x healthy sample")+
   ylab("VAF in the original 300x sample")+
   theme_light()+
   theme(legend.key.width = unit(0.5, "cm"))




temp2 = results %>%
  filter(type == "FP") %>%
  inner_join(callers) %>%
  left_join(snps %>% select(-c(ref, alt)) %>% 
              mutate(SNP=TRUE)) %>%
  mutate(SNP = if_else(is.na(SNP), FALSE, SNP)) %>% 
  filter(SNP) %>%
  left_join(bothhealthy_readcounts)  %>%
  mutate(healthyVAF = if_else(rc_ref_healthy+rc_alt_healthy ==0, 0,rc_alt_healthy/(rc_alt_healthy+rc_ref_healthy) )) %>%
      mutate(superhealthyVAF = if_else(rc_ref_superhealthy+rc_alt_superhealthy ==0, 0,rc_alt_superhealthy/(rc_alt_superhealthy+rc_ref_superhealthy))) 




 
 p2 = ggplot(temp2)+
  geom_point(aes(x = healthyVAF,y = superhealthyVAF, color =name_caller), size = 0.1)+
   stat_ellipse(data = temp2 %>% filter(name_caller=="SNooPer"), aes(x = healthyVAF, y = superhealthyVAF, color = name_caller),type = "norm")+
   scale_color_discrete(name = "Caller")+
   xlab("VAF in the 50x healthy sample")+
 ylab("VAF in the original 300x sample")  +
   annotate("text", x=0.35, y = 0.60, label ="SNooPer", color = "#47A8F8", angle = 30, size = 6)+
   theme_light()+
   theme(legend.key.width = unit(0.5, "cm"))


ggarrange(p3, ggarrange(p1, p2, nrow = 2, labels = c("B","C")), ncol = 2, widths = c(1,3.5), labels = "A")
}
