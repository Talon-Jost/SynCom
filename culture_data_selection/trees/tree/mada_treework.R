rm(list = ls())
library(ape)
library(dplyr)
library(tidyr)

tree.original <- ape::read.tree('tree/nwk/madaDB.rooted.tree.nwk')
plot(tree.original)
print(tree.original$tip.label)

clusters <- read.delim('MadaCultureDB_seqs_clustered_per100.txt', sep = '\t')
clusters.edit <- clusters %>% 
  filter(
    Idmatch != '*'
  )


names <- clusters.edit$Isolate.ID
keep <- intersect(clusters$Isolate.ID, tree.original$tip.label)

tree.pruned <- keep.tip(tree.original, keep)
plot(tree.pruned)

#png("tree_plot.png", width = 1200, height = 10000, res = 150)

ape::plot.phylo(tree.pruned,
                cex = 0.4,
                y.lim = c(0, -5),
                no.margin = TRUE)

#dev.off()

