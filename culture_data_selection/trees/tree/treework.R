#####

library(ape)
library(dplyr)
library(tidyr)

tree.og <- ape::read.tree('rooted.tree.export/tree.nwk')
plot(tree.og)
ape::plot.phylo(tree.og, type = 'fan', cex = 0.4, no.margin = TRUE)


#trim the tree to just the cluster matches
clusters <- read.csv('rooted.tree.export/cluster_results.csv')
names <- trimws(as.character(clusters$AmphibacReference))

print(tree.og$tip.label)

### not sure why the names are different. changing them so they're the same
clusters$AmphibacReference <- gsub("_", ".", clusters$AmphibacReference)
clusters$AmphibacReference <- gsub("", ".", clusters$AmphibacReference)
print(clusters$AmphibacReference)

keep <- intersect(clusters$AmphibacReference, tree.og$tip.label)
tree.pruned <- ape::keep.tip(tree.og, keep)
plot(tree.pruned)
