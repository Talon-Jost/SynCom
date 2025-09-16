####snp work
```{r}
dist <- read.csv("C:/Users/tpj5244/Desktop/Genomes/Genomes/Vences_Genomes/snp_matrix.csv") %>% 
  janitor::clean_names()

colnames(dist) <- gsub("^X\\_", "", colnames(dist), ignore.case = TRUE)

ggplot(dist, aes(x = dist)) +
  geom_dotplot(binwidth = 1000, fill = "blue", dotsize = 0.5) +
  theme_minimal() +
  labs(title = "Dot Plot of SNP Positions",
       x = "Genomic Position",
       y = "Count (binned)")


ggplot(dist, aes(x = p1, y = p2)) +
  geom_point(alpha = 0.6) +
  labs(x = "Position in Genome 1", y = "Position in Genome 2", title = "SNP Position Dot Plot") +
  theme_minimal()


dist$diff <- abs(dist$p1 - dist$p2)
threshold <- 100000
dist$outlier <- dist$diff > threshold

ggplot(dist, aes(x = p1, y = p2, color = outlier)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(
    title = "SNP Position Dot Plot (Outliers Highlighted)",
    x = "Position in Genome 1",
    y = "Position in Genome 2",
    color = "Outlier"
  ) +
  theme_minimal()



```