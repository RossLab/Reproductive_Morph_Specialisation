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

X_Inv_orthogroups <- read.table("r_input/X_vs_Inv_longest_OG_HOG.tsv", header = T, sep="\t")
X_Inv_SCO <- read.table("r_input/X_vs_Inv_longest_Orthogroups_SCOs.tsv")
head(X_Inv_orthogroups)
nrow(X_Inv_orthogroups) #2875
head(X_Inv_SCO)
nrow(X_Inv_SCO) #2703

X_Inv_SCO_genes <- X_Inv_orthogroups[X_Inv_orthogroups$HOG %in% X_Inv_SCO$V1, ]
nrow(X_Inv_SCO_genes) # 2703
X_Inv_SCO_genes <- X_Inv_SCO_genes[, c("Bcop_Inv_genes", "Bcop_X_genes")]
X_Inv_SCO_genes[] <- lapply(X_Inv_SCO_genes, function(x) sub("\\.t[0-9]+$", "", x))

# there is one annoying gene name that does have .t1 in it, so I'm adding it back
X_Inv_SCO_genes[which(X_Inv_SCO_genes$Bcop_X_genes == "g13114"), 2] <- "g13114.t1"
X_Inv_SCO_genes

genes_and_chromosomes <- read.table("r_input/Bcop_v3_genes_and_chromosomes.tsv")
colnames(genes_and_chromosomes) <- c("chrom", "gene")

#### annotated SNPs after mapping to outgroup (separately for X and Inv)  ####
site_counts <- read.table("r_input/degen_site_counts.tsv", header = T)
site_counts <- site_counts[,c(1,6,7)]
colnames(site_counts) <- c('transcript', 'N_nonsyn_sites', 'N_syn_sites')
nrow(site_counts) # 22737 genes
head(site_counts)

# Synonymous and non-synonymous variants for each gene
var_counts_X <- read.table("r_input/Bcop_X.filterstats.genes.relevant.txt", header = T)
# First column is genes and second column is transcript 
nrow(var_counts_X) # 9652
colnames(var_counts_X) <- c("gene", "transcript", "missense", "stop_retained", "synonymous")

var_counts_Inv <- read.table("r_input/Bcop_Inv.filterstats.genes.relevant.txt", header = T)
# First column is genes and second column is transcript 
nrow(var_counts_Inv) # 9709
colnames(var_counts_Inv) <- c("gene", "transcript", "missense", "stop_retained", "synonymous")
```
\
Now that I have the variant counts from the separate versions of the genome (autosomes + X vs autosomes + Inversion), I wanted to compare the autosomal genes first as a santify check. They should be pretty much exactly equal unless the presence of the X or the Inversion causes a higher rate of mismapping and thus problems in variant calling.
```
var_counts_X$N_syn_var <- var_counts_X$stop_retained+var_counts_X$synonymous
var_counts_X$N_nonsyn_var <- var_counts_X$missense
var_counts_X <- var_counts_X[, c(1,2,6,7)]

var_counts_Inv$N_syn_var <- var_counts_Inv$stop_retained+var_counts_Inv$synonymous
var_counts_Inv$N_nonsyn_var <- var_counts_Inv$missense
var_counts_Inv <- var_counts_Inv[, c(1,2,6,7)]

## Merge all together
site_counts
nrow(site_counts) #22737
var_counts_X
var_counts_Inv

sites_and_counts_X <- merge(site_counts, var_counts_X, by=c("transcript"))
sites_and_counts_X_chrom <- merge(sites_and_counts_X, genes_and_chromosomes, by=c("gene"))
sites_and_counts_X_chrom <- na.omit(sites_and_counts_X_chrom)
head(sites_and_counts_X_chrom)
nrow(sites_and_counts_X) #9150

sites_and_counts_Inv <- merge(site_counts, var_counts_Inv, by=c("transcript"))
sites_and_counts_Inv_chrom <- merge(sites_and_counts_Inv, genes_and_chromosomes, by=c("gene"))
sites_and_counts_Inv_chrom <- na.omit(sites_and_counts_Inv_chrom)
head(sites_and_counts_Inv_chrom)
nrow(sites_and_counts_Inv_chrom) #9118

##### calc dN/dS per gene #####
# dN/dS is the ratio of the number of nonsynonymous substitutions per non-synonymous site to the number of synonymous substitutions per synonymous site
# i.e. (N_nonsyn_var/N_nonsyn_sites)/(N_syn_var/N_syn_sites)
sites_and_counts_X_chrom$dN <- sites_and_counts_X_chrom$N_nonsyn_var/sites_and_counts_X_chrom$N_nonsyn_sites
sites_and_counts_X_chrom$dS <- sites_and_counts_X_chrom$N_syn_var/sites_and_counts_X_chrom$N_syn_sites
sites_and_counts_X_chrom$dNdS <- sites_and_counts_X_chrom$dN/sites_and_counts_X_chrom$dS
# clean up values of NaN (0 variants) and Inf, which happens where dS = 0 and dN > 0
sites_and_counts_X_chrom <- sites_and_counts_X_chrom[complete.cases(sites_and_counts_X_chrom), ]
sites_and_counts_X_chrom <- sites_and_counts_X_chrom[which(sites_and_counts_X_chrom$dNdS != Inf),]
nrow(sites_and_counts_X_chrom) #5456

sites_and_counts_Inv_chrom$dN <- sites_and_counts_Inv_chrom$N_nonsyn_var/sites_and_counts_Inv_chrom$N_nonsyn_sites
sites_and_counts_Inv_chrom$dS <- sites_and_counts_Inv_chrom$N_syn_var/sites_and_counts_Inv_chrom$N_syn_sites
sites_and_counts_Inv_chrom$dNdS <- sites_and_counts_Inv_chrom$dN/sites_and_counts_Inv_chrom$dS
# clean up values of NaN (0 variants) and Inf, which happens where dS = 0 and dN > 0
# sites_and_counts_Inv_chrom$dNdS[is.nan(sites_and_counts_Inv_chrom$dNdS)] <- 0
sites_and_counts_Inv_chrom <- sites_and_counts_Inv_chrom[complete.cases(sites_and_counts_Inv_chrom),]
sites_and_counts_Inv_chrom <- sites_and_counts_Inv_chrom[which(sites_and_counts_Inv_chrom$dNdS != Inf),]
nrow(sites_and_counts_Inv_chrom) #5453

## Check that the autosome values are similar 
sites_and_counts_X_chrom_autosomes <- sites_and_counts_X_chrom[sites_and_counts_X_chrom$chrom != "X", ]
sites_and_counts_Inv_chrom_autosomes <- sites_and_counts_Inv_chrom[sites_and_counts_Inv_chrom$chrom != "Inversion", ]
both_methods_autosomes <- merge(sites_and_counts_X_chrom_autosomes, sites_and_counts_Inv_chrom_autosomes, by = "gene")
both_methods_autosomes_subset <- both_methods_autosomes[, c("gene", "dNdS.x", "dNdS.y")]
ggplot(data = both_methods_autosomes_subset, aes(x = dNdS.x, y = dNdS.y)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1)
```
\
Most points are bang on the 1 to 1 line, which is reassuring. 
\
\
Now I want to compare the dNdS value of the X'X homologs. 
```
## Overall comparison
X_Inv_orthogroups <- read.table("r_input/X_vs_Inv_longest_OG_HOG.tsv", header = T, sep="\t")
X_Inv_SCO <- read.table("r_input/X_vs_Inv_longest_Orthogroups_SCOs.tsv")

X_Inv_SCO_genes <- X_Inv_orthogroups[X_Inv_orthogroups$HOG %in% X_Inv_SCO$V1, ]
nrow(X_Inv_SCO_genes) # 2703
X_Inv_SCO_genes <- X_Inv_SCO_genes[, c("Bcop_Inv_genes", "Bcop_X_genes")]
colnames(X_Inv_SCO_genes) <- c("Inv_transcript", "X_transcript")

nrow(dNdS_inversion) #698
dNdS_inversion_subset <- dNdS_inversion[, c("transcript", "dNdS")]
colnames(dNdS_inversion_subset) <- c("Inv_transcript", "Inv_dNdS")
X_Inv_SCO_genes_dNdS <- merge(X_Inv_SCO_genes, dNdS_inversion_subset, by = "Inv_transcript")
nrow(X_Inv_SCO_genes_dNdS) # 553

nrow(dNdS_X) #725
dNdS_X_subset <- dNdS_X[, c("transcript", "dNdS")]
colnames(dNdS_X_subset) <- c("X_transcript", "X_dNdS")
X_Inv_SCO_genes_dNdS <- merge(X_Inv_SCO_genes_dNdS, dNdS_X_subset, by = "X_transcript")
nrow(X_Inv_SCO_genes_dNdS) # 516
head(X_Inv_SCO_genes_dNdS) 
X_Inv_SCO_genes_dNdS$num <- c(1:516)

X_Inv_SCO_genes_dNdS$dNdS_diff <- X_Inv_SCO_genes_dNdS$X_dNdS - X_Inv_SCO_genes_dNdS$Inv_dNdS
plot(X_Inv_SCO_genes_dNdS$dNdS_diff)

X_Inv_SCO_genes_dNdS <- X_Inv_SCO_genes_dNdS |>
  mutate(
    X_transcript = str_remove(X_transcript, "\\.t\\d+$"),
    Inv_transcript = str_remove(Inv_transcript, "\\.t\\d+$")
  )
write.csv(X_Inv_SCO_genes_dNdS, "output/X_Inv_ortholog_dNdS.csv", row.names = FALSE)

X_Inv_diff_genes <- X_Inv_SCO_genes_dNdS[abs(X_Inv_SCO_genes_dNdS$dNdS_diff) != 0, ]
nrow(X_Inv_diff_genes) #486 not equal to 0
```
This will be combined later with allele-biased expression analysis. 
