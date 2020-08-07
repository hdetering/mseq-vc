gdir_denovo="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/de-novo/data/"
gdir_spikein="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/RRSV/Results/"
gdir_plots="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/Figures for the manuscript/PLOS_CompBiol/"

########
# RRSV #
########

AF_Results <- readRDS(paste(gdir_denovo,"df_AFdistances.rds",sep=""))
CallersInfo <- readRDS(paste(gdir_denovo,"df_caller.rds",sep=""))
AF_Results$type <- sapply(AF_Results$name_caller_pub, function(x) CallersInfo$`class.x`[match(x, CallersInfo$name_caller)])
#AF_Results$name_caller_pub <- sub("_","\n",AF_Results$name_caller_pub)
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

SRSV_plot <- ggplot(na.omit(AF_Results_SRSV)) +
  geom_boxplot(aes(name_caller_sorted, Distance, col=type, middle = mean(Distance))) +
  facet_grid(. ~ name_caller_sorted, scales = "free") +
  labs(x="", y="Euclidean distance between\ncalled and simualted VAFs", title = "de novo") +
  scale_color_manual(values = mycols) +
  theme_minimal() +
  ylim(0,55) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

########
# SRSV #
########

AF_Results <- readRDS(paste(gdir_spikein,"RRSV.AFdistances.rds",sep=""))
CallersInfo <- readRDS(paste(gdir_spikein,"RRSV.callers.rds",sep=""))
AF_Results$type <- sapply(AF_Results$name_caller, function(x) CallersInfo$class[match(x, CallersInfo$name_caller)])

AF_Results$name_caller_sorted <- with(AF_Results, reorder(name_caller , Distance, mean , na.rm=T))
AF_Results$type <- factor(AF_Results$type, levels = c("marginal","joint"))
AF_Results_RRSV <- AF_Results

RRSV_plot <- ggplot(na.omit(AF_Results_RRSV)) +
  geom_boxplot(aes(name_caller, Distance, col=type, middle = mean(Distance))) + # middle = mean(Distance) plots mean instead of median I think
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

ggsave(plot = combined_plot,filename = paste(gdir_plots, "DistanceAF.combined.png",sep=""), height = 6, width = 12)
