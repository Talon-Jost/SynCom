#mada 106
#remotes::install_github("Waschina/cobrar")
library(cobrar)
library(tidyverse)
library(dplyr)

df <- readRDS("Documents/Penn_State/Genomes/Mada106.RDS")
str(df)
slotNames(df)

summary(df)
names(df)

attributes(df)
as(df, "mada106s4")
df@metadata
View(as.data.frame(df@metadata))

reactions_df <- data.frame(
  id = df@react_id,
  name = df@react_name,
  lower_bound = df@lowbnd,
  upper_bound = df@uppbnd
)
View(reactions_df)

met_df <- data.frame(
  id = df@met_id,
  name = df@met_name,
  compartment = df@met_comp
)
View(met_df)


##analysing the data some more and making figures?
library(sybil)
library(sybilSBML)
library(glpkAPI)
library(ggplot2)
model <- sybilSBML::readSBMLmod("Downloads/Mada106.xml")
result <- optimizeProb(model)


#check flux balance analysis
# Objective value (e.g., growth rate)
print(paste("Objective value:", result@lp_obj))

# Convert sparse matrix to numeric vector
flux_vector_num <- as.numeric(flux_vector)

# Now create data frame
flux_df <- data.frame(
  reaction = react_id(model),
  flux = flux_vector_num
)

head(flux_df)

ggplot(flux_df, aes(x = flux)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Flux Distribution",
       x = "Flux value",
       y = "Number of reactions")


#single gene knockouts
library(sybil)
getReactionsForGene <- function(model, gene) {
  # Extract GPR rules
  gpr_rules <- gprRules(model)
  
  # Logical vector: does the gene appear in each GPR rule?
  present <- grepl(paste0("\\b", gene, "\\b"), gpr_rules)
  
  # Return indices of reactions with this gene in their GPR
  which(present)
}



knockoutGene <- function(model, gene) {
  gene_rxns <- getReactionsForGene(model, gene)
  
  ko_model <- model
  
  for (rxn in gene_rxns) {
    lowbnd(ko_model)[rxn] <- 0
    uppbnd(ko_model)[rxn] <- 0
  }
  
  return(ko_model)
}


getReactionsForGene <- function(model, gene) {
  # Extract GPR rules
  gpr_rules <- gprRules(model)
  
  # Logical vector: does the gene appear in each GPR rule?
  present <- grepl(paste0("\\b", gene, "\\b"), gpr_rules)
  
  # Return indices of reactions with this gene in their GPR
  which(present)
}


gene_list <- allGenes(model)
ko_results <- numeric(length(gene_list))
names(ko_results) <- gene_list

for (g in gene_list) {
  ko_model <- knockoutGene(model, g)
  ko_result <- optimizeProb(ko_model)
  ko_results[g] <- ko_result@lp_obj
}



# Summarize
ko_df <- data.frame(
  gene = gene_list,
  growth_after_KO = ko_results
)

# Look at top essential genes (growth drops a lot)
ko_df <- ko_df[order(ko_df$growth_after_KO), ]
head(ko_df, 10)


library(ggplot2)

# 1. Distribution of knockout objective values
wildtype_obj <- optimizeProb(model)@lp_obj

df_ko <- data.frame(
  gene = names(ko_results),
  obj_val = ko_results
)

ggplot(df_ko, aes(x = obj_val)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
  geom_vline(xintercept = wildtype_obj, color = "red", linetype = "dashed") +
  labs(
    title = "Distribution of Objective Values After Gene Knockout",
    x = "Objective Value (Growth Rate or Other)",
    y = "Number of Genes"
  ) +
  theme_minimal()

# 2. Identify essential genes
# Define essential genes as those whose knockout objective < 10% of wildtype
essential_genes <- df_ko$gene[df_ko$obj_val < 0.1 * wildtype_obj]

cat("Number of essential genes:", length(essential_genes), "\n")
cat("Essential genes:\n")
print(essential_genes)











##### mada85
mada85.df <- readRDS("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada85/Mada85.RDS")
str(mada85.df)
slotNames(mada85.df)

summary(mada85.df)
names(mada85.df)

attributes(mada85.df)

reactions_df.mada85 <- data.frame(
  id = mada85.df@react_id,
  name = mada85.df@react_name,
  lower_bound = mada85.df@lowbnd,
  upper_bound = mada85.df@uppbnd
)
View(reactions_df.mada85)

met_dfmada85 <- data.frame(
  id = mada85.df@met_id,
  name = mada85.df@met_name,
  compartment = mada85.df@met_comp
)
View(met_dfmada85)


##analysing the data some more and making figures?
library(sybil)
library(sybilSBML)
library(glpkAPI)
library(ggplot2)
model.85 <- sybilSBML::readSBMLmod("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada85/Mada85.xml")
result.85 <- optimizeProb(model.85)


#check flux balance analysis
# Objective value (e.g., growth rate)
print(paste("Objective value:", result.85@lp_obj))

# Convert sparse matrix to numeric vector
flux_vector <- fluxes(result.85)
flux_vector_num.85 <- as.numeric(flux_vector)

# Now create data frame
flux_df.85 <- data.frame(
  reaction = react_id(model.85),
  flux = flux_vector_num.85
)

head(flux_df.85)

ggplot(flux_df.85, aes(x = flux)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Flux Distribution",
       x = "Flux value",
       y = "Number of reactions")


#single gene knockouts
library(sybil)
getReactionsForGene <- function(model, gene) {
  gpr_rules <- gprRules(model)
  present <- grepl(paste0("\\b", gene, "\\b"), gpr_rules)
  which(present)
}

knockoutGene <- function(model, gene) {
  gene_rxns <- getReactionsForGene(model, gene)
  ko_model <- model
  for (rxn in gene_rxns) {
    lowbnd(ko_model)[rxn] <- 0
    uppbnd(ko_model)[rxn] <- 0
  }
  return(ko_model)
}

gene_list.85 <- allGenes(model.85)
ko_results.85 <- numeric(length(gene_list.85))
names(ko_results.85) <- gene_list.85

for (g in gene_list.85) {
  ko_model <- knockoutGene(model.85, g)
  ko_result <- optimizeProb(ko_model)
  ko_results.85[g] <- ko_result@lp_obj
}



# Summarize
ko_df.85 <- data.frame(
  gene = gene_list.85,
  growth_after_KO = ko_results.85
)

# Look at top essential genes (growth drops a lot)
ko_df <- ko_df[order(ko_df$growth_after_KO), ]
head(ko_df, 10)


library(ggplot2)

# 1. Distribution of knockout objective values
wildtype_obj <- optimizeProb(model.85)@lp_obj

df_ko.85 <- data.frame(
  gene = names(ko_results.85),
  obj_val = ko_results.85
)

ggplot(df_ko.85, aes(x = obj_val)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
  geom_vline(xintercept = wildtype_obj, color = "red", linetype = "dashed") +
  labs(
    title = "Distribution of Objective Values After Gene Knockout",
    x = "Objective Value (Growth Rate or Other)",
    y = "Number of Genes"
  ) +
  theme_minimal()

# 2. Identify essential genes
# Define essential genes as those whose knockout objective < 10% of wildtype
essential_genes <- df_ko.85$gene[df_ko.85$obj_val < 0.1 * wildtype_obj]

cat("Number of essential genes:", length(essential_genes), "\n")
cat("Essential genes:\n")
print(essential_genes)






#### comparative gneomics
#full genomes don't return any matches between the two which is very unusual. there should be some overlap between the two genomes themselves in terms of overall genes, especially core genomes
genes.85 <- sybil::allGenes(model.85)
genes.106 <- sybil::allGenes(model)

head(sybil::allGenes(model.85))
head(sybil::allGenes(model))

shared.genes.mada106.85 <- intersect(genes.85, genes.106)
unique.85 <- setdiff(genes.85, genes.106)
unique.106 <- setdiff(genes.106, genes.85)

cat("Genes in both:", length(shared.genes.mada106.85), "\n")
cat("Unique to Mada85:", length(unique.85), "\n")
cat("Unique to XX:", length(unique.106), "\n")


#reaction-level comparison
reacts.85 <- react_id(model.85)
reacts.106 <- react_id(model)

shared_rxns <- intersect(reacts.85, reacts.106)
cat("Shared reactions:", length(shared_rxns), "\n")

unique_85 <- setdiff(reacts.85, reacts.106)
unique_106 <- setdiff(reacts.106, reacts.85)

cat("Unique to Mada85:", length(unique_85), "\n")
cat("Unique to Mada106:", length(unique_106), "\n")

sybsys.85 <- model.85@subSys
subsys.106 <- model@subSys

library(pathview)
all.rxn <- union(reacts.85, reacts.106)
pres.85 <- as.integer(all.rxn %in% reacts.85)
pres.106 <- as.integer(all.rxn %in% reacts.106)
react.matrix <- data.frame(
  reaction = all.rxn,
  mada85 = pres.85,
  mada106 = pres.106
)




pathview::pathview(
  gene.data = react.matrix$mada85,
  pathway.id = "map00010",
  species = "ko",
  gene.idtype = "KEGG",
  limit = list(gene = c(0,1))
)



#####KEGG stuff
mada.106ko <- read.table("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada106/Mada.106_ko.txt",
                         sep = '\t', 
                         fill = TRUE, 
                         quote = "", 
                         stringsAsFactors = FALSE,
                         comment.char = "",
                         header = T)

ko_file <- readLines("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada106/Mada.106_ko.txt")
ko_file <- ko_file[grepl("\t", ko_file)]
ko_df <- read.table(text = ko_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ko_df <- ko_df[grepl("^K\\d{5}$", ko_df$gene), ]
ko.106 <- unique(ko_df$gene)

mada.85ko <- read.table("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada85/Mada.85_ko.txt",
                        sep = '\t',
                        fill = TRUE,
                        quote = "",
                        stringsAsFactors = FALSE,
                        comment.char = "",
                        header = TRUE)

ko_file.85 <- readLines("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada85/Mada.85_ko.txt")
ko_file.85 <- ko_file.85[grepl("\t", ko_file.85)]
ko_df.85 <- read.table(text = ko_file.85, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ko_df.85 <- ko_df.85[grepl("^K\\d{5}$", ko_df.85$gene), ]
ko.85 <- unique(ko_df$gene)

shared.kos <- intersect(ko.85, ko.106)
unique.85 <- setdiff(ko.85, ko.106)
unique.106 <- setdiff(ko.106, ko.85)

cat("Shared KO terms:", length(shared.kos), "\n")
cat("Unique to Mada85:", length(unique.85), "\n")
cat("Unique to Mada106:", length(unique.106), "\n")


ko.85 <- scan("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada85/Mada.85_ko.txt", what = "", sep = "\n")
ko.106 <- scan("Documents/Penn_State/Genomes/Vences_Genomes/Janibacter_Mada106/Mada.106_ko.txt", what = "", sep = "\n")


library(KEGGREST)

library(KEGGREST)
library(pbapply)  # for progress bar

# Fast, safe version
get_pathways_fast <- function(kos) {
  ko_path_map <- pbapply::pblapply(kos, function(ko) {
    tryCatch({
      res <- keggLink("pathway", ko)
      unname(res)
    }, error = function(e) NULL)
  })
  
  # Flatten and keep only valid results
  unique(unlist(ko_path_map))
}

pathways.85 <- get_pathways_fast(ko.85)
pathways.106 <- get_pathways_fast(ko.106)





####### redownload

ko_file <- readLines("Downloads/user_ko.txt")
ko_file <- ko_file[grepl("\t", ko_file)]
ko_df <- read.table(text = ko_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ko_df <- ko_df[grepl("^K\\d{5}$", ko_df$gene), ]
ko.106 <- unique(ko_df$gene)



ko_file.85 <- readLines("Downloads/user_ko(1).txt")
ko_file.85 <- ko_file.85[grepl("\t", ko_file.85)]
ko_df.85 <- read.table(text = ko_file.85, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ko_df.85 <- ko_df.85[grepl("^K\\d{5}$", ko_df.85$gene), ]
ko.85 <- unique(ko_df$gene)


####overlap
essential.85 <- df_ko.85$gene[df_ko.85$obj_val < 0.1 * wildtype_obj]
essential.106 <- df_ko$gene[df_ko$obj_val < 0.1 * wildtype_obj]

essential.85 <- na.omit(essential.85)
essential.106 <- na.omit(essential.106)

# Also ensure they are character vectors, not factors or NULL
essential.85 <- as.character(essential.85)
essential.106 <- as.character(essential.106)

#install.packages("VennDiagram")
library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(Mada85 = essential.85, Mada106 = essential.106),
  filename = NULL,
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  main = "Overlap of Essential Genes"
)
grid.draw(venn.plot)


#####prokka work
library(rtracklayer)

gff106 <- "Documents/Penn_State/Genomes/Vences_Genomes/prokka_annot106/PROKKA_07172025.gff"
read.gff106 <- rtracklayer::import(gff106)

head(read.gff106)

#search for genes such as tRNA
genes <- read.gff106[read.gff106$type == "gene"]
cds <- read.gff106[read.gff106$type == "CDS"]

# View attributes
mcols(cds)

#mada.85
read.gff85 <- "Documents/Penn_State/Genomes/Vences_Genomes/prokka_annot85/PROKKA_07172025.gff"

#read the tsv's in
tsv85 <- read.table("prokka_annot85/PROKKA_07172025.tsv", 
                    sep = '\t',
                    header = T)
tsv106 <- read.table("prokka_annot106/PROKKA_07172025.tsv", 
                     sep = '\t',
                     header = T)

genes_85 <- na.omit(tsv85$gene)
genes_106 <- na.omit(tsv106$gene)

shared_genes <- intersect(genes_85, genes_106)
unique_85 <- setdiff(genes_85, genes_106)
unique_106 <- setdiff(genes_106, genes_85)

cat("Shared genes:", length(shared_genes), "\n")
cat("Unique to Mada85:", length(unique_85), "\n")
cat("Unique to Mada106:", length(unique_106), "\n")

pseudotsv <- read.delim("prokka_Pseudomonas/pseudo_annot/PROKKA_01302025.tsv", sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
genes_pseudo <- na.omit(pseudotsv$gene)

shared.pseudo.85 <- intersect(genes_85, genes_pseudo)
unique_85.a <- setdiff(genes_85, genes_pseudo)
unique_pseudo <- setdiff(genes_pseudo, genes_85)

cat("Shared genes:", length(shared.pseudo.85), "\n")
cat("Unique to Mada85:", length(unique_85.a), "\n")
cat("Unique to Mada106:", length(unique_pseudo), "\n")


###filter and save
unique.pseudotsv <- unique_pseudo
pseudo.filtered.tsv <- pseudotsv[pseudotsv$gene %in% unique.pseudotsv, ]
write.csv(pseudo.filtered.tsv, file = "unique_genes_pseudomonas.csv", col.names = F)

unique.mada85 <- unique_85
mada85.filtered.tsv <- tsv85[tsv85$gene %in% unique.mada85, ]
write.csv(mada85.filtered.tsv, file = "unique_genes_mada85.csv", col.names = F)

unique.mada106 <- unique_106
mada106.filtered.tsv <- tsv106[tsv106$gene %in% unique_106, ]
write.csv(mada106.filtered.tsv, file = "unique_genes_mada106.csv", col.names = F)


####reading in the fixed files that I did by hand
setwd('Documents/Penn_State/Genomes/Vences_Genomes/')
fixed106.unique <- read.csv('Janibacter_Mada106/unique_genes_mada106.csv')
fixed85.unique <- read.csv('Janibacter_Mada85/unique_genes_mada85.csv')

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

#remotes::install_github("jokergoo/circlize")
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


###FROM FASTA FILE  
library(Biostrings)
library(circlize)
mada106fasta <- readDNAStringSet(filepath = "Janibacter_Mada106/Mada106_6148_6699.fna")
mada85fasta <- readDNAStringSet(filepath = 'Janibacter_Mada85/Mada85.fna')


len1 <- width(mada106fasta)
len2 <- width(mada85fasta)

# Standardize length to max genome for alignment
maxlen <- max(len1, len2)

# Initialize circle
circos.clear()
circos.par(start.degree = 90, gap.degree = 5, track.height = 0.1)
circos.initialize(factors = c("Genome1", "Genome2"),
                  xlim = rbind(c(0, len1), c(0, len2)))

# Genome 1 backbone
circos.track(ylim = c(0, 1), factors = "Genome1", bg.border = "black",
             panel.fun = function(x, y) {
               circos.axis(h = "bottom", labels.cex = 0.5)
             })

# Genome 2 backbone
circos.track(ylim = c(0, 1), factors = "Genome2", bg.border = "black",
             panel.fun = function(x, y) {
               circos.axis(h = "bottom", labels.cex = 0.5)
             })

circos.clear()
