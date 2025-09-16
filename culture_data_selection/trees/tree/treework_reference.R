rm(list = ls())
library(ape)
library(dplyr)
library(tidyr)

tree.original <- ape::read.tree('Desktop/tree.nwk')
plot(tree.original)

clusters <- read.csv('Desktop/cluster_results.csv')
names <- clusters$AmphibacReference

print(tree.original$tip.label)
clusters$quoted_name <- paste0("'", clusters$AmphibacReference, "'")
keep <- intersect(clusters$quoted_name, tree.original$tip.label)
tree.pruned <- keep.tip(tree.original, keep)
plot(tree.pruned)

#this is pretty good, but we can make it better
ape::plot.phylo(tree.pruned, type = 'fan', cex = 0.4, no.margin = TRUE)

#add color coding, but we need the inhibitory metadata
meta.data <- read.csv('Desktop/AmphibBac_Metadata_2023.2.txt', sep = '\t')
meta.data <- meta.data %>% 
  dplyr::select(c('FastaID', 'AntiFungal.Function.with.AP.function')) %>% 
  dplyr::rename(AmphibacReference = 'FastaID')
meta.data <- meta.data[meta.data$AmphibacReference %in% clusters$AmphibacReference,]

cluster.refined <- clusters %>% 
  left_join(meta.data, by = 'AmphibacReference')
cluster.refined <- cluster.refined[match(tree.pruned$tip.label, 
                                         paste0("'", cluster.refined$AmphibacReference, "'")), ]
group_factor <- as.factor(cluster.refined$AntiFungal.Function.with.AP.function)
group_colors <- rainbow(length(levels(group_factor)))[group_factor]


color.map <- c(
  'Inhibitory' = 'red',
  'NonInhibitory' = 'green',
  'PossibleContamination' = 'blue',
  'NoEffect' = 'violet',
  'Facilitating' = 'yellow',
  'NotTested' = 'black'
)

group_colors <- color.map[as.character(group_factor)]

#png("Desktop/circular_tree.png", width = 3000, height = 3000, res = 300)
#
ape::plot.phylo(tree.pruned, 
                type = 'fan', 
                cex = 0.2, 
                tip.color = group_colors,
                no.margin = TRUE)

legend("topright",                          
       legend = names(color.map),           
       fill = color.map,                    
       border = NA,
       bty = "n",                           
       cex = 1) 
#dev.off()


## barchart work
meta.cleaned <- cluster.refined %>%
  mutate(NonZero_Columns = strsplit(as.character(NonZero_Columns), ",\\s*")) %>%
  mutate(num_funcs = lengths(NonZero_Columns)) %>%
  unnest_wider(NonZero_Columns, names_sep = "_")

meta.cleaned <- meta.cleaned %>%
  mutate(obs.count = rowSums(!is.na(.[, 3:63])))

meta.cleaned$obs.count <- as.numeric(meta.cleaned$obs.count)

## make into barplot tree
tree.df <- data.frame(tip = tree.pruned$tip.label) 
tree.df <- left_join(tree.df, meta.cleaned[, c("AmphibacReference", "obs.count")], by = c("tip" = "AmphibacReference"))

# Now create the bar chart data
color_scale <- colorRampPalette(c("blue", "green", "yellow", "red"))(100)  # You can adjust the color range here
norm_count <- (tree.df$obs.count - min(tree.df$obs.count)) / (max(tree.df$obs.count) - min(tree.df$obs.count))  # Normalize obs.count to [0, 1]

# Map normalized counts to colors
tree.df$color <- color_scale[round(norm_count * 99) + 1]

# Set up the plotting device
#png("Desktop/tree_with_barchart.png", width = 3000, height = 3000, res = 300)

# Plot the tree
new <- plot.phylo(tree.pruned, 
                type = "fan",
                fill = group_colors,
                cex = 0.3, 
                no.margin = TRUE)
new

# Add bar chart
for (i in 1:nrow(tree.df)) {
  rect(xleft = i - 0.3, ybottom = 0, xright = i + 0.3, ytop = tree.df$obs.count[i],
       col = tree.df$color[i], border = NA)
}

# Save the plot
#dev.off()

write.csv(meta.cleaned, 'meta_cleaned.csv')


tree <- ggtree::ggtree(tree.pruned,
               branch.length = 'none',
               layout = 'circular') #+ geom_tiplab()
amp.meta <- read.csv('Desktop/MadaCultureDB_joined_culture_with_metadata_2025.csv') %>% 
  dplyr::select(c('FastaID', 'DirTax_Order'))

amp.meta.2 <- amp.meta[amp.meta$FastaID %in% cluster.refined$AmphibacReference,]

amp.meta.2 <- amp.meta.2 %>% 
  filter(!duplicated(FastaID)) %>% 
  dplyr::rename(AmphibacReference = 'FastaID')

amp.meta.3 <- cluster.refined %>% 
  left_join(amp.meta.2, by = 'AmphibacReference')

tree$data <- tree$data %>%
  left_join(amp.meta.3, by = c("label" = "AmphibacReference"))
table(is.na(tree))
tree + 
  aes(color = DirTax_Order) + 
  geom_tree()
