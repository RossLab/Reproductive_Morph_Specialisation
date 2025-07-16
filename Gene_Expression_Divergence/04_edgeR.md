# Differential gene expression analysis #
Differential gene expression (DGE) analysis was carried out in R (v4.3.2), using the edgeR package.
\
The necessary packages are: 
```
library(edgeR) # v4.0.9
library(PCAtools) # v2.14.0
library(pheatmap) # v1.0.12
library(tidyverse) # v2.0.0
library(ggpubr) # v0.6.0
```

## Comparison of autosomal genes between androgenic and gynogenic females ##

The analysis was done separate for each tissue type. 
\
\
A design matrix is needed (for example, the table in 00_Data.md) to specify which group (sex/tissue) the sample belongs to. For the purpose of reading in the featureCount files that are in separate folders, I am also adding the file path to the design matrix. 
```
## Design matrix ##
design <- read.csv("inputs/gene_expression_divergence_gyno_andro_male_design.csv")
design$path <- paste0("inputs/", design$sample_id, ".txt",sep="") # Allows featureCount files to be read in more easily.
group <- factor(paste(design$sex,design$tissue,sep=".")) # Makes a separate group that combines information about the sex and tissue. Mostly for PCA plots. 
design <- cbind(design,group=group)
```
\
Because I want to subset for autosomal genes, I read in a file that specifies the contig that each gene belongs to. In this case the autosomes are II, III and IV, and the sex chromosomes are X and Inversion (aka the X'). 

```
## Genes and the chromosome they belong to ##
genes_and_chromosomes <- read.table("inputs/Bcop_v3_genes_and_chromosomes.tsv")
autosomal_genes <- genes_and_chromosomes[genes_and_chromosomes$V1 == "II" | genes_and_chromosomes$V1 == "III" | genes_and_chromosomes$V1 == "IV", ]
```
\
I then summarise the separate featureCounts txt files into one big dataframe, with one gene per row and one sample per column, populated by the read count of said gene in said sample. In the process I also subset for genes that are on autosomal contigs.

```
## Summarising featureCount read files into one dataframe ##
autosomal_featureCounts <- lapply(design$path, function(x)read.table(x, header=T)) # Applies the read.table() function to each file path, reading each file into a data frame
names(autosomal_featureCounts) <- design$sample_id # Names each dataframe in the list
autosomal_featureCounts <- lapply(autosomal_featureCounts, function(x) x[, c(1, 7)]) # For each data frame in the list, only columns 1 and 7 are kept, which is the gene name and the read count.
autosomal_featureCounts <- lapply(names(autosomal_featureCounts), function(name) {
  x <- autosomal_featureCounts[[name]]
  colnames(x)[2] <- name
  x 
}) # Assigns the corresponding sample ID (name) as the name of the second column (the count data column)
names(autosomal_featureCounts) <- design$sample_id # Names each dataframe in the list because the name got lost
autosomal_counts <- data.frame(GeneID = autosomal_genes$V2) # here I am subsetting for only autosomal genes 
colnames(autosomal_counts) <- "Geneid"
for (name in names(autosomal_featureCounts)) {
  # Extract the current table
  x <- autosomal_featureCounts[[name]]
  # Extract column 2, setting it as a new column in `autosomal_counts` with the sample ID as the column name
  autosomal_counts <- merge(autosomal_counts, x, by = "Geneid", all.x = TRUE)
}
rownames(autosomal_counts) <- autosomal_counts$Geneid
autosomal_counts <- autosomal_counts[, 2:28]
nrow(autosomal_counts) # I have a total of 16679 rows, which is the correct number of autosomal genes
```
\
Now that I have a dataframe with all my read counts, I can start to set up a DGEList, which is the edgeR object that includes both read counts and design matrix. I first do it with all the tissues and all sexes together. 
I filter lowly expressed genes, as they will not have enough information to statistically determine whether they are differentially expressed anyways. I normalise for different library sizes and composition between different samples, which is necessary for DGE analysis. 

```
#### DGE analysis setup ####

## Read the counts into a DGEList, specifying the design ##
all_DGEList <- DGEList(autosomal_counts, group = design$group)
nrow(all_DGEList$counts)

## Filter lowly expressed genes ##
keep <- filterByExpr(all_DGEList)
table(keep) # Keep 12889 genes, get rid of 3790 genes 
all_DGEList <- all_DGEList[keep, , keep.lib.sizes=FALSE]
nrow(all_DGEList$counts) # After filtering, I am left with 12889 genes, which is the correct number

## Normalisation ##
all_DGEList <- normLibSizes(all_DGEList)
head(all_DGEList)
print(all_DGEList$samples) # Check that library sizes have changed 

## Visualise sample distances to check that replicates cluster as you expect ##
plotMDS(all_DGEList)
```
\
Plot and save the PCA plot.


```
####################### PCA plots ###########################
counts_for_pca<-cpm(all_DGEList,log=TRUE,prior.count=1) 

# Set group with custom labels 
all_DGEList$samples$group <- factor(all_DGEList$samples$group, levels = c(
  "androgenic.germline",
  "androgenic.somatic_reproductive",
  "androgenic.somatic_non_reproductive",
  "gynogenic.germline",
  "gynogenic.somatic_reproductive",
  "gynogenic.somatic_non_reproductive",
  "male.germline",
  "male.somatic_reproductive",
  "male.somatic_non_reproductive"
),
labels = c(
  "Androgenic - Germline",
  "Androgenic - Somatic reproductive",
  "Androgenic - Somatic non-reproductive",
  "Gynogenic - Germline",
  "Gynogenic - Somatic reproductive",
  "Gynogenic - Somatic non-reproductive",
  "Male - Germline",
  "Male - Somatic reproductive",
  "Male - Somatic non-reproductive"
))

pca_output <- pca(counts_for_pca, metadata = all_DGEList$samples)
screeplot(pca_output) # Use a screeplot to check the amount of variation that can be explained by each PC

png(file="./plots/allsexalltissue_pca.png", height = 800, width = 1000)
all_tissue_pca <- biplot(pca_output, colby = "group", colkey = c(
  "Androgenic - Germline" = "#a1d99b",
  "Androgenic - Somatic reproductive" = "#41ab5d",
  "Androgenic - Somatic non-reproductive" = "#006d2c",
  "Gynogenic - Germline" = "#d4b9da",
  "Gynogenic - Somatic reproductive" = "#af8dc3",
  "Gynogenic - Somatic non-reproductive" = "#762a83",
  "Male - Germline" = "#9ecae1",
  "Male - Somatic reproductive" = "#4292c6",
  "Male - Somatic non-reproductive" = "#084594"), legendPosition = "bottom", lab = NULL) + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 18)) +
  guides(
    colour = guide_legend(
      nrow = 3, ncol = 3,
      byrow = TRUE
    )
  )
all_tissue_pca
dev.off()
```
\
Having checked that the samples are behaving as they should, I move on to analysing differential gene expression between androgenic and gynogenic females at each tissue level. The first step is to repeat the setup, but only with samples from one tissue. I am first analysing the somatic non-reproductive tissue. 


```
##### SOMATIC NON-REPRODUCTIVE TISSUE #####
## Making a DGEList object ##
design_nonrepro <- design[c(7:9, 16:18, 25:27), ]
counts_nonrepro <- autosomal_counts[, c(7:9, 16:18, 25:27)]
DGEList_nonrepro <- DGEList(counts = counts_nonrepro, group = design_nonrepro$sex)
head(DGEList_nonrepro)
nrow(DGEList_nonrepro$counts) #16679

# Constructing a model matrix for glm 
design_nonrepro_matrix <- model.matrix(~0 + design_nonrepro$sex)
rownames(design_nonrepro_matrix) <- design_nonrepro$sample_id
colnames(design_nonrepro_matrix) <- sort(unique(design_nonrepro$sex))
design_nonrepro_matrix

## Filter lowly expressed genes ##
keep <- filterByExpr(DGEList_nonrepro)
table(keep) # Keep 11656 genes, get rid of 5023 genes 
DGEList_nonrepro <- DGEList_nonrepro[keep, , keep.lib.sizes=FALSE]
nrow(DGEList_nonrepro$counts) #11656

## Normalisation ##
DGEList_nonrepro <- normLibSizes(DGEList_nonrepro)
head(DGEList_nonrepro)

## Visualise sample distances ##
plotMDS(DGEList_nonrepro)
counts_for_pca<-cpm(DGEList_nonrepro,log=TRUE,prior.count=1) 
pca_output <- pca(counts_for_pca, metadata = DGEList_nonrepro$samples)

## PCA plot for paper ##
png(file="./plots/somaticnonrepro_pca.png", height = 800, width = 800)
non_repro_pca <- biplot(pca_output, colby = "group", legendPosition = "none", lab = NULL, colkey = c(
  "androgenic" = "#006d2c",
  "gynogenic" = "#762a83",
  "male" = "#084594"), title = "C", titleLabSize = 25) + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust = -0.1, vjust = -1.5))
dev.off()
screeplot(pca_output)
```
\
Common dispersion and tagwise dispersion is estimated, then analysed for differentially expressed genes (DEGs) by a glm QLFT test. 

```
# Estimate common dispersion and tagwise dispersion in one go 
DGEList_nonrepro <- estimateDisp(DGEList_nonrepro, design = design_nonrepro_matrix)
# This is a plot of the negative binomial dispersion for the dataset, to check that it is to expectations. 
plotBCV(DGEList_nonrepro)

# Given raw counts, NB dispersion(s) and a design matrix, glmQLFit() fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
fit_nonrepro <- glmQLFit(DGEList_nonrepro, design_nonrepro_matrix)
# Compares androgenics (as the baseline) to gynogenics 
qlf.androvsgyno_nonrepro <- glmQLFTest(fit_nonrepro, contrast=c(-1, 1, 0))
qlf.androvsgyno_nonrepro_FDR <- topTags(qlf.androvsgyno_nonrepro, 11656) # Saves the output in a dataframe. 11656 is because that is the number of genes that is being considered after filtering for sufficient expression.
qlf.androvsgyno_nonrepro_FDR_sig <- qlf.androvsgyno_nonrepro_FDR$table[qlf.androvsgyno_nonrepro_FDR$table$FDR < 0.05,] # Filter for significantly DEG with a false discovery rate (FDR) of less than 5%
summary(decideTests(qlf.androvsgyno_nonrepro)) # Summarises the number of genes upregulated, downregulated or unchanged. 
write.csv(qlf.androvsgyno_nonrepro_FDR_sig, "outputs/nonrepro_androvsgyno_sig.csv")
```
\
Visualise results with heatmaps/volcano plots. 
```
## Heatmap ## 
logcpm_nonrepro <- cpm(DGEList_nonrepro, log=TRUE) # Extract cpm values (normalised gene counts)
logcpm_nonrepro_females <- logcpm_nonrepro[, 4:9] # Extract only values for androgenic and gynogenic females
heatmap_all_nonrepro_50 <- pheatmap(logcpm_nonrepro_females[rownames(head(qlf.androvsgyno_nonrepro_FDR_sig, 50)), ],
                                    cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T,
                                    main = "Top 50 significant DEGs between female morphs in nonreproductive tissue")
 
### Volcano plots ###
# Assign expression groups for legend 
qlf.androvsgyno_nonrepro_FDR_df <- as.data.frame(qlf.androvsgyno_nonrepro_FDR)
qlf.androvsgyno_nonrepro_FDR_df <- qlf.androvsgyno_nonrepro_FDR_df %>% 
  mutate(
    Expression = case_when(FDR < 0.05 & logFC <= 0 ~ "Androgenic-biased",
                           FDR < 0.05 & logFC >= 0 ~ "Gynogenic-biased",
                           TRUE ~ "Unchanged")
  )
table(qlf.androvsgyno_nonrepro_FDR_df$Expression) # Should match summary(decideTests(qlf.androvsgyno_nonrepro))

nonrepro_volcano <- ggplot(qlf.androvsgyno_nonrepro_FDR_df, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 1) +
  scale_color_manual(values = c("chartreuse3", "purple", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept = -log(0.05,10), colour = "gray") +
  # geom_vline(xintercept = c(1, -1), colour="gray") +
  ggtitle("F") + theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = -0.1, vjust = -1.5))
```
