Necessary packages include: 

```
library(edgeR) # v4.0.9
library(PCAtools) # v2.14.0
library(pheatmap) # v1.0.12
library(tidyverse) # v2.0.0
library(ggpubr) # v0.6.0
library(cowplot)
```
\
Differential deposition of maternal transcripts are analysed in the same way as in 01_Gene_Expression_Divergence/04_edgeR.md
```
#### Read in files, and wrangle them into useable formats ####
## Design matrix ##
design <- read.csv("r_input/maternal_deposit_design.csv")
design$path <- paste0("r_input/", design$sample_id, ".multi.featureCounts.txt",sep="")
group <- factor(paste(design$sex,design$time,sep="."))
design <- cbind(design,group=group)
design

## Genes and the chromosome they belong to ##
genes_and_chromosomes <- read.table("r_input/Bcop_v3_genes_and_chromosomes.tsv")
autosomal_genes <- genes_and_chromosomes[genes_and_chromosomes$V1 == "II" | genes_and_chromosomes$V1 == "III" | genes_and_chromosomes$V1 == "IV", ]

## Summarising featureCount read files into one dataframe ##
autosomal_featureCounts <- lapply(design$path, function(x)read.table(x, header=T)) # applies the read.table() function to each file path, reading each file into a data frame
names(autosomal_featureCounts) <- design$sample_id # Names each dataframe in the list
autosomal_featureCounts <- lapply(autosomal_featureCounts, function(x) x[, c(1, 7)]) # For each data frame in the list, only columns 1 and 7 are kept.
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
autosomal_counts <- autosomal_counts[, 2:13]
nrow(autosomal_counts) # 16679

# The dataset originally had eggs collected at 4-6 hours as well, but I am only interested in the eggs from 0-4 hours, which is when I expect there to be more maternally deposited transcript.
####################### 0-4 hours #######################
## Making a DGEList object ##
design_04 <- design[c(1:3, 7:9), ]
counts_04 <- autosomal_counts[, c(1:3, 7:9)]
DGEList_04 <- DGEList(counts = counts_04, group = design_04$sex)
head(DGEList_04)
nrow(DGEList_04$counts) #16679

# Constructing a model matrix for glm 
design_04_matrix <- model.matrix(~0 + design_04$sex)
rownames(design_04_matrix) <- design_04$sample_id
colnames(design_04_matrix) <- sort(unique(design_04$sex))
design_04_matrix

## Filter lowly expressed genes ##
keep <- filterByExpr(DGEList_04)
table(keep) # Keep 7452 genes, get rid of 9227 genes 
DGEList_04 <- DGEList_04[keep, , keep.lib.sizes=FALSE]
nrow(DGEList_04$counts) #7452

## Normalisation ##
DGEList_04 <- normLibSizes(DGEList_04)
head(DGEList_04)

## Visualise sample distances ##
plotMDS(DGEList_04)
counts_for_pca<-cpm(DGEList_04,log=TRUE,prior.count=1) 
pca_output <- pca(counts_for_pca, metadata = DGEList_04$samples)
pca_output$metadata$group<- factor(pca_output$metadata$group, levels = c("gynogenic", "androgenic"))

# png(file="./outputs/maternaldeposit_0-4_pca.png", height = 500, width = 500)
pca_04 <- biplot(pca_output, colby = "group", legendPosition = "bottom", gridlines.major = FALSE, lab = NULL, title = "A - Autosomal PCA", titleLabSize = 25) +
  coord_fixed(ratio = 1, xlim = c(-25, 25), ylim = c(-25, 25)) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 23, hjust = 0, margin = margin(b = 5), face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.margin = margin(2.5, 2.5, 2.5, 2.5),
    legend.position = c(0.75, 0.87),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=14),
    legend.background = element_rect(fill = NULL, color = "black", linewidth = 0.5, linetype = "solid")) + 
  scale_color_manual(labels = c("Gynogenic", "Androgenic"), values = c("purple","chartreuse3")) +
  labs(col = "Reproductive morph") +
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2)) 
pca_04
# dev.off()
screeplot(pca_output)

# Estimate common dispersion and tagwise dispersion in one go 
DGEList_04 <- estimateDisp(DGEList_04, design = design_04_matrix)
# This is a plot of the NB dispersion for the dataset 
plotBCV(DGEList_04)

# Given raw counts, NB dispersion(s) and a design matrix, glmQLFit() fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
fit_04 <- glmQLFit(DGEList_04, design_04_matrix)
head(fit_04$coefficients)
# compares androgenics (as the baseline) to gynogenics 
qlf.androvsgyno_04 <- glmQLFTest(fit_04, contrast=c(-1, 1))
qlf.androvsgyno_04_FDR <- topTags(qlf.androvsgyno_04, 7452)
qlf.androvsgyno_04_FDR_sig <- qlf.androvsgyno_04_FDR$table[qlf.androvsgyno_04_FDR$table$FDR < 0.05,]
summary(decideTests(qlf.androvsgyno_04))
# write.csv(qlf.androvsgyno_04_FDR_sig, "outputs/04_androvsgyno_sig.csv")

## Plotting heatmaps etc 
logcpm_04 <- cpm(DGEList_04, log=TRUE)
logcpm_04
heatmap_all_04_50 <- pheatmap(logcpm_04[rownames(head(qlf.androvsgyno_04_FDR_sig,43)), ],
                                    cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T,
                                    main = "Significant autosomal DEGs between female morphs in eggs from 0-4 hrs")

#### Volcano plots (which requires differential gene analysis) ####
# compares androgenics (as the baseline) to gynogenics 
summary(decideTests(qlf.androvsgyno_04))
plotMD(qlf.androvsgyno_04)
abline(h=c(-1,1), col="blue") 


qlf.androvsgyno_04_FDR_df <- as.data.frame(qlf.androvsgyno_04_FDR)

qlf.androvsgyno_04_FDR_df <- qlf.androvsgyno_04_FDR_df %>% 
  mutate(
    Expression = case_when(FDR < 0.05 & logFC <= 0 ~ "Androgenic-biased",
                           FDR < 0.05 & logFC >= 0 ~ "Gynogenic-biased",
                           TRUE ~ "Unbiased")
  )
table(qlf.androvsgyno_04_FDR_df$Expression)

volcano_04 <- ggplot(qlf.androvsgyno_04_FDR_df, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 1.5) +
  scale_color_manual(values = c("chartreuse3", "purple", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept = -log(0.05,10), colour = "gray") +
  # geom_vline(xintercept = c(1, -1), colour="gray") +
  ggtitle("B - Autosomal DGE") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title.position = "plot",
    plot.title = element_text(size = 23, hjust = 0, margin = margin(b = 10), face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title=element_text(size=16), 
    legend.text=element_text(size=16)
  )+
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2))

ggarrange(pca_04, volcano_04, align = "hv")
# dev.off()
```
