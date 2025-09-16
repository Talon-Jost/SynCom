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

library(ape)
ape::plot.phylo(tree.pruned,
                cex = 0.2,
                y.lim = c(0, 20))




ape::plot.phylo(tree.pruned,
                cex = 0.3,
                y.lim = c(0, -5),
                no.margin = TRUE)

ape::plot.phylo(tree.original,
                cex = 0.3,
                y.lim = c(0, -5),
                no.margin = TRUE)
