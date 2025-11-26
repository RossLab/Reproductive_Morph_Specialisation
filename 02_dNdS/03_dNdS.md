# Calculate dNdS per gene #
The last step is to combine the measures of the number of synonymous sites and nonsynonymous sites for each gene, with the number of synonymous and nonsynoymous mutations for each gene. I have done this in R (v4.3.2).

Necessary packages:
```
library(ggplot2) # v 3.5.1
library(tidyr) # v 1.3.1
library(dplyr) # v 1.1.4
```
```
### Input files ###
# Synonymous and nonsynonymous sites (degen_site_counts.tsv), and also synonymous and nonsynonymous mutations (Bcop.filterstats.genes.relevant.txt)
# also:
# List of single-copy orthologs between species pairs -> Bcop_SCO_genes_list.txt
# Which genes are on which scaffold for each focal species -> Bcop_v3_genes_and_chromosomes.tsv
# List of longest transcripts -> As_and_X_longest.headers.txt

Bcop_SCO <- read.table("r_input/Bcop_SCO_genes_list.txt") # 8687 genes
Bcop_SCO <- Bcop_SCO[order(Bcop_SCO$V1),]
# If later I want to remove the transcript names
# Bcop_SCO$V1 <- sub("\\.t1$", "", Bcop_SCO$V1)

genes_and_chromosomes <- read.table("r_input/Bcop_v3_genes_and_chromosomes.tsv")
colnames(genes_and_chromosomes) <- c("chrom", "gene")

longest_transcripts <- read.table("r_input/all_longest.headers.txt")

site_counts <- read.table("r_input/degen_site_counts.tsv", header = T)
site_counts <- site_counts[,c(1,6,7)]
colnames(site_counts) <- c('transcript', 'N_nonsyn_sites', 'N_syn_sites')
nrow(site_counts) # 12414 genes
head(site_counts)
site_counts_longest <- site_counts[site_counts$transcript %in% longest_transcripts$V1, ]
nrow(site_counts_longest) # 11440 

# Synonymous and non-synonymous variants for each gene
var_counts_X <- read.table("r_input/Bcop_X.filterstats.genes.relevant.txt", header = T)
# First column is genes and second column is transcript 
nrow(var_counts_X) # 9882 genes
# subset for longest isoforms
var_counts_X_longest <- var_counts_X[var_counts_X$TranscriptId %in% longest_transcripts$V1, ]
nrow(var_counts_X_longest) # 9621 genes
colnames(var_counts_X_longest) <- c("gene", "transcript", "missense", "stop_retained", "synonymous")

var_counts_Inv <- read.table("r_input/Bcop_Inv.filterstats.genes.relevant.txt", header = T)
# First column is genes and second column is transcript 
nrow(var_counts_Inv) # 9940 genes
# subset for longest isoforms
var_counts_Inv_longest <- var_counts_Inv[var_counts_Inv$TranscriptId %in% longest_transcripts$V1, ]
nrow(var_counts_Inv_longest) # 9682 genes
colnames(var_counts_Inv_longest) <- c("gene", "transcript", "missense", "stop_retained", "synonymous")
```
\
Now that I have the variant counts from the separate versions of the genome (autosomes + X vs autosomes + Inversion), I wanted to compare the autosomal genes first as a santify check. They should be pretty much exactly equal unless the presence of the X or the Inversion causes a higher rate of mismapping and thus problems in variant calling.
```
var_counts_X_longest$N_syn_var <- var_counts_X_longest$stop_retained+var_counts_X_longest$synonymous
var_counts_X_longest$N_nonsyn_var <- var_counts_X_longest$missense
var_counts_X_longest <- var_counts_X_longest[, c(1,2,6,7)]

var_counts_Inv_longest$N_syn_var <- var_counts_Inv_longest$stop_retained+var_counts_Inv_longest$synonymous
var_counts_Inv_longest$N_nonsyn_var <- var_counts_Inv_longest$missense
var_counts_Inv_longest <- var_counts_Inv_longest[, c(1,2,6,7)]

## Merge all these data frames together
sites_and_counts_X <- merge(site_counts_longest, var_counts_X_longest, by=c("transcript"))
sites_and_counts_X_chrom <- merge(sites_and_counts_X, genes_and_chromosomes, by=c("gene"))
sites_and_counts_X_chrom <- na.omit(sites_and_counts_X_chrom)

sites_and_counts_Inv <- merge(site_counts_longest, var_counts_Inv_longest, by=c("transcript"))
sites_and_counts_Inv_chrom <- merge(sites_and_counts_Inv, genes_and_chromosomes, by=c("gene"))
sites_and_counts_Inv_chrom <- na.omit(sites_and_counts_Inv_chrom)

# calc dN/dS per gene i.e. (N_nonsyn_var/N_nonsyn_sites)/(N_syn_var/N_syn_sites)
sites_and_counts_X_chrom$dN <- sites_and_counts_X_chrom$N_nonsyn_var/sites_and_counts_X_chrom$N_nonsyn_sites
sites_and_counts_X_chrom$dS <- sites_and_counts_X_chrom$N_syn_var/sites_and_counts_X_chrom$N_syn_sites
sites_and_counts_X_chrom$dNdS <- sites_and_counts_X_chrom$dN/sites_and_counts_X_chrom$dS
# clean up values of NaN (0 variants) and Inf, which happens where dS = 0 and dN > 0
sites_and_counts_X_chrom$dNdS[is.nan(sites_and_counts_X_chrom$dNdS)] <- 0
sites_and_counts_X_chrom <- sites_and_counts_X_chrom[which(sites_and_counts_X_chrom$dNdS != Inf),]

# Now do the same for the dataframe with the inversion 
sites_and_counts_Inv_chrom$dN <- sites_and_counts_Inv_chrom$N_nonsyn_var/sites_and_counts_Inv_chrom$N_nonsyn_sites
sites_and_counts_Inv_chrom$dS <- sites_and_counts_Inv_chrom$N_syn_var/sites_and_counts_Inv_chrom$N_syn_sites
sites_and_counts_Inv_chrom$dNdS <- sites_and_counts_Inv_chrom$dN/sites_and_counts_Inv_chrom$dS
# clean up values of NaN (0 variants) and Inf, which happens where dS = 0 and dN > 0
sites_and_counts_Inv_chrom$dNdS[is.nan(sites_and_counts_Inv_chrom$dNdS)] <- 0
sites_and_counts_Inv_chrom <- sites_and_counts_Inv_chrom[which(sites_and_counts_Inv_chrom$dNdS != Inf),]

## Check that the autosome values are similar 
sites_and_counts_X_chrom_autosomes <- sites_and_counts_X_chrom[sites_and_counts_X_chrom$chrom != "X", ]
sites_and_counts_Inv_chrom_autosomes <- sites_and_counts_Inv_chrom[sites_and_counts_Inv_chrom$chrom != "Inversion", ]
both_methods_autosomes <- merge(sites_and_counts_X_chrom_autosomes, sites_and_counts_Inv_chrom_autosomes, by = "gene")
both_methods_autosomes_subset <- both_methods_autosomes[, c("gene", "dNdS.x", "dNdS.y")]
ggplot(data = both_methods_autosomes_subset, aes(x = dNdS.x, y = dNdS.y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
```
\
Most points are bang on the 1 to 1 line, which is reassuring. 
\
\
Now I want to compare the dNdS value of the X'X homologs. 
```
## Overall comparison
dNdS_inversion <- sites_and_counts_Inv_chrom[sites_and_counts_Inv_chrom$chrom == "Inversion", ]
dNdS_X <- sites_and_counts_X_chrom[sites_and_counts_X_chrom$chrom == "X", ]
boxplot(dNdS_inversion$dNdS, dNdS_X$dNdS, ylim = c(0,1)) # no big overall differences
mean(dNdS_inversion$dNdS) #0.1392735
mean(dNdS_X$dNdS) #0.1554073

## Now read in orthofinder output for a homolog-by-homolog comparison
X_Inv_orthogroups <- read.table("r_input/X_vs_Inv_HOG.tsv", header = T, sep="\t")
X_Inv_SCO <- read.table("r_input/X_vs_Inv_Orthogroups_SCOs.tsv")

X_Inv_SCO_genes <- X_Inv_orthogroups[X_Inv_orthogroups$HOG %in% X_Inv_SCO$V1, ]
nrow(X_Inv_SCO_genes) # 2152
X_Inv_SCO_genes <- X_Inv_SCO_genes[, c("Bcop_Inv_genes", "Bcop_X_genes")]
colnames(X_Inv_SCO_genes) <- c("Inv_transcript", "X_transcript")

nrow(dNdS_inversion) #464
dNdS_inversion_subset <- dNdS_inversion[, c("transcript", "dNdS")]
colnames(dNdS_inversion_subset) <- c("Inv_transcript", "Inv_dNdS")
X_Inv_SCO_genes_dNdS <- merge(X_Inv_SCO_genes, dNdS_inversion_subset, by = "Inv_transcript")
nrow(X_Inv_SCO_genes_dNdS) # 272

nrow(dNdS_X) #524
dNdS_X_subset <- dNdS_X[, c("transcript", "dNdS")]
colnames(dNdS_X_subset) <- c("X_transcript", "X_dNdS")
X_Inv_SCO_genes_dNdS <- merge(X_Inv_SCO_genes_dNdS, dNdS_X_subset, by = "X_transcript")
nrow(X_Inv_SCO_genes_dNdS) # 171
head(X_Inv_SCO_genes_dNdS) 
X_Inv_SCO_genes_dNdS$num <- c(1:171)

X_Inv_SCO_genes_dNdS_long <- X_Inv_SCO_genes_dNdS %>%
  pivot_longer(cols = c(X_dNdS, Inv_dNdS), 
               names_to = "Chrom", 
               values_to = "dNdS") %>%
  mutate(Chrom = ifelse(Chrom == "X_dNdS", "X", "Inversion"))

ggplot(X_Inv_SCO_genes_dNdS_long, aes(x = Chrom, y = dNdS, group = num)) +
  geom_boxplot(aes(group = Chrom), width = 0.5, alpha = 0.3) +  # Boxplots
  geom_line(aes(group = num), color = "gray", alpha = 0.5) +  # Connect orthologs
  geom_point(aes(color = Chrom), size = 2) +  # Points for each value
  theme_minimal() +
  labs(y = "dN/dS", x = "Chrom", title = "dN/dS Comparison of X and Inversion") +
  theme(legend.position = "none")

X_Inv_SCO_genes_dNdS$dNdS_diff <- X_Inv_SCO_genes_dNdS$X_dNdS - X_Inv_SCO_genes_dNdS$Inv_dNdS
plot(X_Inv_SCO_genes_dNdS$dNdS_diff)

X_Inv_SCO_genes_dNdS <- X_Inv_SCO_genes_dNdS |>
  mutate(
    X_transcript = str_remove(X_transcript, "\\.t\\d+$"),
    Inv_transcript = str_remove(Inv_transcript, "\\.t\\d+$")
  )
write.csv(X_Inv_SCO_genes_dNdS, "output/X_Inv_ortholog_dNdS.csv", row.names = FALSE)

X_Inv_diff_genes <- X_Inv_SCO_genes_dNdS[abs(X_Inv_SCO_genes_dNdS$dNdS_diff) != 0, ]
nrow(X_Inv_diff_genes) #123 not equal to 0
```
This will be combined later with homolog-specific expression analysis. 
