gdir_denovo="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/de-novo/data/"
gdir_spikein="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/RRSV/Results/"
gdir_plots="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/Figures for the manuscript/PLOS_CompBiol/"

########
# SRSV #
########

AF_Results <- readRDS(paste(gdir_denovo,"df_AFdistances.rds",sep=""))
CallersInfo <- readRDS(paste(gdir_denovo,"df_caller.rds",sep=""))
CallersInfo <- CallersInfo[which(CallersInfo$name_caller!="MuClone_perf"),]
colnames(CallersInfo)[2] <- "name_caller_pub"


AF_Results <- AF_Results %>%
  dplyr::left_join(CallersInfo, by = "name_caller_pub") %>%
  dplyr::select (name_caller_pub, Replicate, Distance, class.x) %>%
  dplyr::rename(type=class.x)

AF_Results$name_caller_sorted <- with(AF_Results, reorder(name_caller_pub , Distance, mean , na.rm=T))
# ('hs': highly structured/low admixture, 'ms': moderately structured, 'us': unstructured/high admixture)

# Select the same two colors than for the rest of plots even with a levels missing
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 3
cols = gg_color_hue(n)
mycols=c(cols[1],cols[3])

AF_Results$type <- factor(AF_Results$type,levels=c("marginal","joint"))
AF_Results_SRSV <- AF_Results

SRSV_plot <- ggplot(na.omit(AF_Results_SRSV),aes(name_caller_sorted, Distance)) +
  geom_boxplot(aes(col=type)) +
  stat_summary(fun=mean, geom="point", 
               shape=22,fill="black") +
  facet_grid(. ~ name_caller_sorted, scales = "free") +
  labs(x="", y="Euclidean distance between\nsimulated and estimated VAFs", title = "de novo") +
  scale_color_manual(values = mycols) +
  theme_minimal() +
  ylim(0,55) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

########
# RRSV #
########

AF_Results <- readRDS(paste(gdir_spikein,"RRSV.AFdistances.rds",sep=""))
CallersInfo <- readRDS(paste(gdir_spikein,"RRSV.callers.rds",sep=""))

AF_Results <- AF_Results %>%
  dplyr::left_join(CallersInfo, by = "name_caller") %>%
  dplyr::select (name_caller, Replicate, Distance, class) %>%
  dplyr::rename(type=class)

AF_Results$name_caller_sorted <- with(AF_Results, reorder(name_caller , Distance, mean , na.rm=T))
AF_Results$type <- factor(AF_Results$type, levels = c("marginal","joint"))
AF_Results_RRSV <- AF_Results

RRSV_plot <- ggplot(na.omit(AF_Results_RRSV),aes(name_caller, Distance,)) +
  geom_boxplot(aes(col=type)) + 
  stat_summary(fun=mean, geom="point", 
               shape=22,fill="black") +
  facet_grid(. ~ name_caller_sorted, scales = "free") +
  labs(x="", y="",col="", title = "spike-in") +
  scale_color_manual(values = mycols) +
  scale_fill_manual(values = mycols) +
  theme_minimal() +
  ylim(0,55) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

####################
# Combined results #
####################

combined_plot <- ggarrange(plotlist = list(SRSV_plot,RRSV_plot),
                           labels = c("a","b"),
                           common.legend = T,
                           legend = "bottom")

ggsave(plot = combined_plot,filename = paste(gdir_plots, "DistanceAF.combined.pdf",sep=""), width = 7, height = 5, dpi = 300)


###########################
# Correlation F1-distance #
###########################

#----- SRSV --------------#

df_perf_sample <- readRDS( file.path(gdir_denovo, 'df_perf.rds') )
F1_SRSV<- df_perf_sample %>%
  dplyr::select(name_caller,name_rep,F1) %>%
  dplyr::rename(Replicate=name_rep) %>%
  dplyr::rename(name_caller_sorted=name_caller)
dist_F1 <- AF_Results_SRSV %>%
  dplyr::select(-type,-name_caller_pub) %>%
  dplyr::inner_join(F1_SRSV, by=c("name_caller_sorted", "Replicate")) %>%
  reshape2::melt(id.vars=c("Replicate","name_caller_sorted","Distance"))

dist_F1_forcor <- dist_F1[which(dist_F1$variable=="F1"),]
correlation <- cor.test(dist_F1_forcor$Distance, dist_F1_forcor$value)

# Define the number of colors you want
library(RColorBrewer)
#nb.cols <- 13
#mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
mycolors <- c(brewer.pal(12, "Paired"),"grey")


cor_F1_dist_SRSV <- ggplot(dist_F1 %>% dplyr::filter(variable %in% c("F1")) %>% dplyr::filter(!name_caller_sorted %in% c("Shimmer","SNV-PPILP","MuClone")),aes(x=Distance, y=as.numeric(value))) +
  geom_point(aes(color=name_caller_sorted)) +
  geom_smooth(method='lm',col="black") +
  scale_colour_manual(values = mycolors) +
  ylim(0,1) +
  labs(x="Euclidean distances between\nestimated and simulated VAFs", y="F1",color="caller", title = "de novo") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x=8,y=0.8, label=paste("p.value = ",formatC(correlation$p.value, format = "e", digits = 2),"\nr = ",round(correlation$estimate,2),sep=""), size=3)

#----- RRSV --------------#

perf_RRSV <- readRDS( file.path(gdir_spikein, 'df_perf.rds') )
F1_RRSV<- perf_RRSV %>%
  #select(name_caller,name_rep,F1) %>%
  dplyr::select(-id_caller,-id_rep,-class,-nclones,-nsamples,-caller,-cvg,-ttype) %>%
  dplyr::rename(Replicate=name_rep) %>%
  dplyr::rename(name_caller_sorted=name_caller)
dist_F1 <- AF_Results_RRSV %>%
  dplyr::select(-type,-name_caller) %>%
  dplyr::inner_join(F1_RRSV, by=c("name_caller_sorted", "Replicate")) %>%
  reshape2::melt(id.vars=c("Replicate","name_caller_sorted","Distance"))

dist_F1_forcor <- dist_F1[which(dist_F1$variable=="F1"),]
correlation <- cor.test(dist_F1_forcor$Distance, dist_F1_forcor$value)

cor_F1_dist_RRSV <- ggplot(dist_F1 %>% dplyr::filter(variable %in% c("F1")) %>% dplyr::filter(!name_caller_sorted %in% c("Shimmer","SNV-PPILP","MuClone")),aes(x=Distance, y=as.numeric(value))) +
  geom_point(aes(color=name_caller_sorted)) +
  geom_smooth(method='lm',col="black") +
  scale_colour_manual(values = mycolors) +
  ylim(0,1) +
  labs(x="Euclidean distances between\nestimated and simulated VAFs", y="", color="caller", title = "spike-in") +
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x=40,y=0.8, label=paste("p.value = ",formatC(correlation$p.value, format = "e", digits = 2),"\nr  = ",round(correlation$estimate,2),sep=""), size=3)



combined_cor_plot <- ggarrange(plotlist = list(cor_F1_dist_SRSV,cor_F1_dist_RRSV),
                           labels = c("a","b"),
                           common.legend = T,
                           legend = "right")

ggsave(plot = combined_cor_plot,filename = paste(gdir_plots, "DistanceAF-F1.correlation.pdf",sep=""), width = 7, height = 5, dpi = 300)


#ggplot(dist_F1 %>% filter(variable %in% c("F1","precision","recall")),aes(x=Distance, y=as.numeric(value))) +
#  geom_point(aes(color=name_caller_sorted,shape=name_caller_sorted)) +
#  #geom_smooth() +  # method='lm') +
#  labs(x="Euclidean distance between\nsimulated and estimated VAFs", y="F1",color="caller") +
#  theme_linedraw() +
#  facet_grid(. ~ variable, scales = "free")
#
#ggplot(dist_F1 %>% filter(variable %in% c("FN","FP","TP")) %>% filter(!name_caller_sorted %in% c("Shimmer","SNV-PPILP","MuClone")),
#       aes(x=name_caller_sorted, y=as.numeric(value))) +
#  geom_col(aes(fill=name_caller_sorted)) +
#  #geom_smooth() +  # method='lm') +
#  labs(x="caller", y="Counts",color="caller") +
#  theme_linedraw() +
#  facet_grid(variable ~ ., scales = "free") +
#  theme(axis.text.x = element_text(angle = 45, hjust=1),
#        legend.position = "none")
  


