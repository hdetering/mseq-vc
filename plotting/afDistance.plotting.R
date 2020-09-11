gdir_denovo="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/de-novo/data/"
gdir_spikein="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/RRSV/Results/"
gdir_plots="/Users/tama/Google Drive/PHYLOGENOMICS/M-seq Variant Calling Benchmarking/Figures for the manuscript/PLOS_CompBiol/"

gdir_denovo <- file.path( 'data', 'de-novo' )

########
# RRSV #
########

AF_Results <- readRDS( file.path(gdir_denovo, "df_AFdistances.rds") )
CallersInfo <- readRDS( file.path(gdir_denovo, "df_caller.rds") )
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

data_dir <- file.path( 'data', 'de-novo' )

True_VAFs_muts <- readRDS( file.path(data_dir, "df_mut_sample.rds") )
mutations <- readRDS( file.path(data_dir, "df_mut.rds") )
Replicates_Info <- readRDS( file.path(data_dir, "df_rep.rds") )

df_vaf_exp <- True_VAFs_muts %>%
  inner_join( mutations, by = c('id_rep', 'id_mut') ) %>%
  inner_join( Replicates_Info, by = c('id_rep') ) %>%
  unite( 'chrom_pos', chrom, pos) %>%
  select( name_rep, id_sample, chrom_pos, vaf_exp, rc_tot, rc_alt )

# True_VAFs_muts$chrom <- sapply(True_VAFs_muts$id_mut, function(x)  mutations$chrom[match(x, mutations$id_mut)])
# True_VAFs_muts$pos <- sapply(True_VAFs_muts$id_mut, function(x)  mutations$pos[match(x, mutations$id_mut)])
# True_VAFs_muts$ref <- sapply(True_VAFs_muts$id_mut, function(x)  mutations$ref[match(x, mutations$id_mut)])
# True_VAFs_muts$alt <- sapply(True_VAFs_muts$id_mut, function(x)  mutations$alt[match(x, mutations$id_mut)])
# True_VAFs_muts$name_rep <- sapply(True_VAFs_muts$id_rep, function(x)  Replicates_Info$name_rep[match(x, Replicates_Info$id_rep)])
# # create a new column mutinfo with the four columns collapsed together
# cols <- c( 'chrom' , 'pos' , 'ref', 'alt' )
# True_VAFs_muts$mut_info <- apply( True_VAFs_muts[ , cols ] , 1 , paste , collapse = "_" )
# # create a new column mutinfo with the four columns collapsed together
# cols <- c( 'chrom' , 'pos')
# True_VAFs_muts$chrom_pos <- apply( True_VAFs_muts[ , cols ] , 1 , paste , collapse = "_" )

AFs <- readRDS( file.path(data_dir, "SRSV.AFs.rds") )
AFs$name_rep <- as.character( AFs$replicate )
AFs <- AFs %>% rename( R1 = T1, R2 = T2, R3 = T3, R4 = T4, R5 = T5 )
df_vaf_obs <- AFs %>% pivot_longer( cols = paste0('R', 1:5), names_to = 'id_sample', values_to = 'vaf' ) %>%
  mutate( vaf_obs = as.numeric(vaf) ) %>%
  select( -mut_info, -replicate, -vaf )
# df_vaf_exp <- True_VAFs_muts %>%
#   unite( 'chrom_pos', chrom, pos ) %>%
#   dplyr::filter( id_sample != 'RN' ) %>%
#   dplyr::select( name_rep, id_sample, chrom_pos, vaf_exp )
df_vaf <- df_vaf_obs %>% 
  dplyr::left_join( df_vaf_exp, by = c('name_rep', 'id_sample', 'chrom_pos') )
df_vaf_dist <- df_vaf %>%
  drop_na( c(name_rep, software, vaf_obs) ) %>%
  mutate( sq_err = (vaf_obs - vaf_exp)^2 ) %>%
  group_by( name_rep, software ) %>%
  dplyr::summarise( mean_sq_err = mean(sq_err, na.rm = T) )

SRSV_plot <- ggplot(na.omit(df_vaf_dist)) +
  geom_boxplot(aes(name_rep, mean_sq_err, middle = mean(mean_sq_err))) +
  labs(x="", y="Mean squared error between\ncalled and simulated VAFs", title = "de novo") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

#-------------------------------------------------------------------------------

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
