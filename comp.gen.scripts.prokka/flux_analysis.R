---
title: "flux_analysis"
author: "Talon Jost"
date: "2025-08-29"
output: pdf_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```



```{r cars}
library(dplyr)
library(tidyr)

```


The files that are here (and uploaded in onedrive)
```{r pressure, echo=FALSE}
fixed106.unique <- read.csv('unique_genes_mada106.csv')
fixed85.unique <- read.csv('unique_genes_mada85.csv')

int <- intersect(fixed106.unique$gene, fixed85.unique$gene)
int <- as.data.frame(int)
int <- int %>% 
  dplyr::rename(gene = "int")
int.merged <- int %>% 
  dplyr::left_join(fixed106.unique, by = "gene")

unique.int106 = fixed106.unique %>% 
  dplyr::anti_join(int.merged, by = "gene")

unique.int85 <- fixed85.unique %>% 
  dplyr::anti_join(int.merged, by = "gene")

genes <- sort(unique(c(fixed106.unique$gene, fixed85.unique$gene)))

gene.matrix <- data.frame(
  genes = genes,
  Mada106 = ifelse(genes %in% fixed106.unique$gene, 1, 0),
  Mada85 = ifelse(genes %in% fixed85.unique$gene, 1, 0)
)

gene.matrix.num <- as.matrix(gene.matrix[,-1])
rownames(gene.matrix.num) <- gene.matrix$genes
annot.column <- data.frame(
  inhibition = c('Inhibitory', 'Noninhibitory')
)

rownames(annot.column) <- colnames(gene.matrix.num)

# Unique to Mada106
unique106 <- (gene.matrix.num[, "Mada106"] == 1) & (gene.matrix.num[, "Mada85"] == 0)

# Unique to Mada85
unique85 <- (gene.matrix.num[, "Mada106"] == 0) & (gene.matrix.num[, "Mada85"] == 1)

# Combine into a factor for row annotation
row_annot <- data.frame(
  Unique = factor(ifelse(unique106, "Mada106",
                         ifelse(unique85, "Mada85", "Shared")))
)
rownames(row_annot) <- rownames(gene.matrix.num)

library(pheatmap)
heatmap1 <- pheatmap(
  gene.matrix.num,
  cluster_rows = FALSE,      # cluster genes
  cluster_cols = FALSE,      # cluster organisms
  show_rownames = FALSE,    # hide row names if too many
  color = c("white", "steelblue"),# 0 = white, 1 = blue
  annotation_col = annot.column, 
  annotation_row = row_annot,
  filename = "mada85106_heatmap.pdf"
)
heatmap1
ggsave(heatmap1, plot = 'mada85106_heatmap.pdf', limitsize = FALSE, width = 30, height = 50)
#this is actually a very interesting figure to look at because I didn't expect so many unique genes
#I find this curious as there should be a number of conserved genes given this is a closely related species
#perhaps trying to make a circular genome with it might help?
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

```{r}
remotes::install_github("jokergoo/circlize") #might have to force this install
library(circlize)

genome <- fixed106.unique

genes$start <- cumsum(c(0, head(genome$length_bp, -1))) + 1
genes$end <- cumsum(genome$length_bp)

genome.length <- sum(genome$length_bp)

genome$type <- ifelse(genome$gene %in% unique.int106$gene, "Mada106_unique",
                      ifelse(genome$gene %in% int$gene, "Shared", "Other")

circos.clear()
circos.initialize(factors = rep("chr1", nrow(genome)), xlim = c(0, genome.length))

for(i in 1:nrow(genome)) {
  circos.rect(
    xleft = genome$start[i],
    xright = genome$end[i],
    ybottom = 0,
    ytop = 1,
    col = ifelse(genome$type[i] == "Shared", "gray",
                 ifelse(genome$type[i] == "Mada106_unique", "red", "steelblue")),
    border = "black"
  )
}

circos.initialize("chr1", xlim = c(0, genome.length))

col_vec <- ifelse(genome$type == "Shared", "gray",
                  ifelse(genome$type == "Mada106_unique", "red", "steelblue"))

# Single-track plotting without a panel loop
circos.trackPlotRegion(
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    for(i in 1:nrow(genome)) {
      circos.rect(
        xleft = genome$start[i],
        xright = genome$end[i],
        ybottom = 0,
        ytop = 1,
        col = ifelse(genome$type[i] == "Shared", "gray",
                     ifelse(genome$type[i] == "Mada106_unique", "red", "steelblue")),
        border = "black"
      )
    }
  }
)

title("Circular Genome Plot - Mada106")                    
#doesn't work

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("RCircos")
library(RCircos)



genome <- fixed106.unique

# Calculate start/end positions based on length_bp
genome$start <- cumsum(c(0, head(genome$length_bp, -1))) + 1
genome$end <- cumsum(genome$length_bp)

# Label genes as shared or unique
genome$type <- ifelse(genome$gene %in% unique.int106$gene, "Mada106_unique",
                      ifelse(genome$gene %in% intersect(fixed106.unique$gene, fixed85.unique$gene), "Shared", "Other"))

# Create karyotype table
karyotype <- data.frame(Chromosome="chr1",
                        chromStart=genome$start,
                        chromEnd=genome$end,
                        Gene=genome$gene,
                        Type=genome$type)

# Create a fake chromosome length table (only one chromosome)
cytoband <- data.frame(Chromosome="chr1", chromStart=0, chromEnd=sum(genome$length_bp), Band="p", Stain="gpos100")

# Initialize RCircos
RCircos.Set.Core.Components(cytoband, chr.exclude=NULL, tracks.inside=5, tracks.outside=0)

# Start the plot
RCircos.Set.Plot.Area()
RCircos.Draw.Chromosome.Ideogram()

# Highlight genes by type
colors <- c(Mada106_unique="red", Shared="gray", Other="steelblue")
Rcircos(karyotype, track.num=1, side="in")
RCircos.Draw.Gene.Heatmap(karyotype, data.col=5, track.num=2, side="in", col=colors)

```




