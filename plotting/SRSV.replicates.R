require(tidyverse)
require(dbplyr)
require(RSQLite)
require(cowplot)
require(RColorBrewer)
require(gridExtra)
require(grid)
require(ggtree)

sims_root = "/home/harry/projects/m-seq_varcall/sim_prep"

#------------------------------------------------------------------------------
# Read Simulation Meta Data
#------------------------------------------------------------------------------

# connection to analysis database
db <- './data/analysis.db'
con <- DBI::dbConnect(RSQLite::SQLite(), db)
df_rep <- tbl( con, 'replicates' ) %>% collect()
df_mut <- tbl( con, 'mutations' ) %>% collect()
df_mut_sample <- tbl( con, 'mutations_samples' ) %>% collect()
df_mut_clone <- tbl( con, 'clones_mutations' ) %>% collect()
df_prev <- tbl( con, 'clones_prev' ) %>% collect()
df_caller <- tbl( con, 'callers' ) %>% collect()
df_varcall <- tbl( con, 'varcalls' ) %>% collect()
df_rc <- tbl( con, 'readcounts' ) %>% collect()


df.trees <- tibble(
  id            = character(),
  node          = integer(),
  parent        = integer(),
  branch.length = numeric(),
  x = numeric(), y = numeric(),
  label         = character(),
  isTip         = logical(),
  branch        = numeric(),
  angle         = numeric()
)

# read clone trees
repdirs <- dir( sims_root, include.dirs = F )
for (rep.dir in repdirs) {
  fn_clone_tree <- file.path(sims_root, rep.dir, 'clone_tree.nwk')
  tree <- phytools::read.newick(fn_clone_tree)
  treedata <- fortify(tree)
  treedata$rep <- rep.dir
  df.trees <- df.trees %>% rbind(treedata)
}
df_tree <- df.trees %>% inner_join( df_rep, by = c('rep'='name_rep') )


#------------------------------------------------------------------------------
# Read Variant Caller Results
#------------------------------------------------------------------------------

df.vars <- readRDS("df_loc_perf.rds") %>%
  replace_na(list(precision = 0, F1 = 0))
# define order of variant callers (will affect plots)
df.vars$caller = factor(df.vars$id_caller, levels=c(
  "Bcftools",
  "HaplotypeCaller",
  "CaVEMan",
  "MuTect1",
  "Mutect2",
  "NeuSomatic",
  "Shimmer",
  "SNooPer",
  "SomaticSniper",
  "Strelka1",
  "Strelka2",
  "VarDict",
  "VarScan",
  "MuClone",
  "MultiSNV",
  "Mutect2_mseq",
  "SNV-PPILP"))

#------------------------------------------------------------------------------
# Plot Scenarios and Caller Results
#------------------------------------------------------------------------------

# define callers for which results will be plotted
callers = c(
  "Bcftools",
  "HaplotypeCaller",
  "CaVEMan",
  "MuTect1",
  "Mutect2",
  "NeuSomatic",
  "Shimmer",
  "SNooPer",
  "SomaticSniper",
  "Strelka2",
  "VarDict",
  "VarScan",
  "MuClone",
  "MultiSNV",
  "Mutect2_mseq",
  "SNV-PPILP")

# open PDF document to receive all plots
# NOTE: will produce empty first page, remove manually.
pdf('plots/SRSV.replicates.pdf', onefile=T, paper='A4r', width = 8, height = 4.5)

# create plots for each replicate
for ( rid in df_rep$id_rep ) {

rep_name <- df_rep %>% dplyr::filter( id_rep == rid ) %>% select( name_rep )
# prevalence matrix
df <- df_prev %>% dplyr::filter( id_rep == id_rep ) #%>% inner_join(df.reps %>% select(id, idx_rep, ttype, cvg), by='id')
p.prev <- ggplot( df, aes(x = id_clone, y = id_sample) ) +
  geom_tile( aes(fill = prev) ) +
  scale_fill_gradient( low = 'white', high = 'red', name = 'prevalence', limits = c(0, 1)) +
  ggtitle( 'Clone Prevalence' ) + theme_bw()
#ggsave("sim_prep.prev.pdf", plot=p, device=pdf, width=8, height=40, units='in')

# clone tree
df <- df.trees %>% dplyr::filter(rep == rep_name) #%>% inner_join(df.reps %>% select(id, idx_rep, ttype, cvg), by='id')
df <- df_tree %>% dplyr::filter(id_rep == rid)
x.max <- max(df$x)
y.max <- max(df$y)
p.tree <- ggtree(df) + geom_label(aes(label=label), size=2) +
  xlim(0, (x.max+.1)) + ylim(.9, (y.max+.1)) +
  ggtitle("Clone Tree") + theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
#ggsave("sim_prep.tree.pdf", plot=p, device=pdf, width=8, height=40, units='in')

# variant caller performance
df.vars.rep <- df.vars %>% dplyr::filter(id_rep == rid & name_caller %in% callers)
p_rec <- ggplot(df.vars.rep, aes(x=name_caller, y=recall)) + 
  geom_bar(aes(fill=name_caller), stat = "identity") + ylim(0, 1) +
  labs(x = "caller") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_gray() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  guides(fill = "none")

p_pre <- ggplot(df.vars.rep, aes(x = name_caller, y = precision)) + 
  geom_bar(aes(fill = name_caller), stat = "identity") + ylim(0, 1) +
  labs(x = "caller") +
  theme_gray() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  guides(fill = "none")

p_f1 <- ggplot(df.vars.rep, aes(x = name_caller, y = F1)) + 
  geom_bar(aes(fill = name_caller), stat = "identity") + ylim(0, 1) +
  labs(x = "caller") +
  theme_gray() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_title <- sprintf("Rep_id: %s, Cvg: %sx, Ttype: %s", rep_name, 
                   df_rep %>% dplyr::filter(id_rep == rid) %>% select(cvg),
                   df_rep %>% dplyr::filter(id_rep == rid) %>% select(ttype))
p.perf <- arrangeGrob(grobs = list(p_rec, p_pre, p_f1), 
                      layout_matrix = rbind(c(1,2), c(3,3), c(3,3)),
                      top = "Variant Calling Performance")

p_rep <- grid.arrange(grobs = list(p.tree, p.prev, p.perf),
                      layout_matrix = rbind(c(1,3,3), c(2,3,3)),
                      top = textGrob(p_title, gp = gpar(fontsize=16,font=3)))
#ggsave(sprintf("reps/%s.performance.local.pdf", rep), plot = p_rep, width = 8, height = 4.5)
#ggsave(sprintf('plots/reps/%s.png', rep_name), plot = p_rep, device = 'png', width = 8, height = 4.5)
}

# save PDF doc
dev.off()
