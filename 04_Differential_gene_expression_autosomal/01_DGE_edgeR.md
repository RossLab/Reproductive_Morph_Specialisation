# Differential gene expression analysis of autosomal genes, in edgeR #
To look at divergent expression in autosomal genes, I can use the same files generated for differential gene expression analysis in X-linked genes. Instead of filtering for X-linked genes, I am filtering for autosomal genes this time. I also don't need to add up the read counts from separate gametologs since I am mapping to a reference genome which has just one "copy" of autosomal genes. The analysis below is largely similar otherwise, and also carried out in R (v4.3.2), using the edgeR package.

\
The necessary packages are: 
```
library(edgeR) # v4.0.9
library(PCAtools) # v2.14.0
library(pheatmap) # v1.0.12
library(tidyverse) # v2.0.0
library(ggpubr) # v0.6.0
```
I first wrangle the data, which is a separate file for each sample, into a single table. In the process, I also subset for only autosomal genes.
```
## Design matrix ##
design <- read.csv("inputs/gene_expression_divergence_gyno_andro_male_design.csv")
design$path <- paste0("inputs/", design$sample_id, ".txt",sep="")
group <- factor(paste(design$sex,design$tissue,sep="."))
design <- cbind(design,group=group)
design

## Genes and the chromosome they belong to ##
genes_and_chromosomes <- read.table("inputs/Bcop_v3_genes_and_chromosomes.tsv")
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
autosomal_counts <- autosomal_counts[, 2:28]
nrow(autosomal_counts) # 16679
```
I look first at the sample distances for all tissue types together. 
```
####################### DGE analysis setup #######################

## Read the counts into a DGEList, specifying the design ##
all_DGEList <- DGEList(autosomal_counts, group = design$group)
nrow(all_DGEList$counts)

## Filter lowly expressed genes ##
keep <- filterByExpr(all_DGEList)
table(keep) # Keep 12889 genes, get rid of 3790 genes 
all_DGEList <- all_DGEList[keep, , keep.lib.sizes=FALSE]
nrow(all_DGEList$counts) #12889

## Normalisation ##
all_DGEList <- normLibSizes(all_DGEList)
head(all_DGEList)
print(all_DGEList$samples)

## Visualise sample distances ##
plotMDS(all_DGEList)

####################### PCA plots ###########################
counts_for_pca<-cpm(all_DGEList,log=TRUE,prior.count=1) 

# Set group with desired order and custom labels
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
screeplot(pca_output)

# png(file="./plots/allsexalltissue_pca.png", height = 800, width = 1000)
all_tissue_pca <- biplot(pca_output, colby = "group", colkey = c(
  "Androgenic - Germline" = "#a1d99b",
  "Androgenic - Somatic reproductive" = "#41ab5d",
  "Androgenic - Somatic non-reproductive" = "#006d2c",
  "Gynogenic - Germline" = "#d4b9da",
  "Gynogenic - Somatic reproductive" = "#af8dc3",
  "Gynogenic - Somatic non-reproductive" = "#762a83",
  "Male - Germline" = "#9ecae1",
  "Male - Somatic reproductive" = "#4292c6",
  "Male - Somatic non-reproductive" = "#084594"), legendPosition = "bottom", lab = NULL, title = "A", titleLabSize = 25) + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 18, margin = margin(t = 0, b = 0)), plot.title = element_text(hjust = -0.05, vjust = -1.5)) +
  guides(colour = guide_legend(nrow = 3, ncol = 3, byrow = TRUE))
all_tissue_pca
# dev.off()
```
\
Now I compare autosomal gene expression between gynogenic and androgenic females in somatic non-reproductive tissue. 
```
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

# png(file="./plots/somaticnonrepro_pca.png", height = 800, width = 800)
non_repro_pca <- biplot(pca_output, colby = "group", legendPosition = "none", gridlines.major = FALSE, lab = NULL, 
                        title = "C - Somatic non-reproductive", titleLabSize = 20)  +
  scale_colour_manual(values = c(
      "androgenic" = "chartreuse3",
      "gynogenic" = "purple",
      "male" = "blue2"), labels = c(
      "androgenic" = "Androgenic female",
      "gynogenic" = "Gynogenic female",
      "male" = "Male")) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme(plot.title.position = "plot",
    plot.title = element_text(size = 22, hjust = 0, vjust = -4),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none",
    legend.title = element_text(size = 18), legend.text = element_text(size = 18, margin = margin(t = 0, b = 0))
  ) 
non_repro_pca

non_repro_pca_legend <- biplot(pca_output, colby = "group", legendPosition = "none", gridlines.major = FALSE, lab = NULL, 
                        title = "C - Somatic non-reproductive", titleLabSize = 20)  +
  scale_colour_manual(name = "Reproductive morph", values = c(
    "androgenic" = "chartreuse3",
    "gynogenic" = "purple",
    "male" = "blue2"), labels = c(
      "androgenic" = "Androgenic female",
      "gynogenic" = "Gynogenic female",
      "male" = "Male")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(plot.title.position = "plot",
        plot.title = element_text(size = 22, hjust = 0, vjust = -4),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_text(size = 18), legend.text = element_text(size = 18, margin = margin(t = 0, b = 0))
  ) 

# dev.off()
screeplot(pca_output)

# Estimate common dispersion and tagwise dispersion in one go 
DGEList_nonrepro <- estimateDisp(DGEList_nonrepro, design = design_nonrepro_matrix)
# This is a plot of the NB dispersion for the dataset 
plotBCV(DGEList_nonrepro)

# Given raw counts, NB dispersion(s) and a design matrix, glmQLFit() fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
fit_nonrepro <- glmQLFit(DGEList_nonrepro, design_nonrepro_matrix)
head(fit_nonrepro$coefficients)
# compares androgenics (as the baseline) to gynogenics 
qlf.androvsgyno_nonrepro <- glmQLFTest(fit_nonrepro, contrast=c(-1, 1, 0))
qlf.androvsgyno_nonrepro_FDR <- topTags(qlf.androvsgyno_nonrepro, 11656)
qlf.androvsgyno_nonrepro_FDR_sig <- qlf.androvsgyno_nonrepro_FDR$table[qlf.androvsgyno_nonrepro_FDR$table$FDR < 0.05,]
summary(decideTests(qlf.androvsgyno_nonrepro))
# write.csv(qlf.androvsgyno_nonrepro_FDR_sig, "outputs/nonrepro_androvsgyno_sig.csv")

## Plotting heatmaps etc 
logcpm_nonrepro <- cpm(DGEList_nonrepro, log=TRUE)
logcpm_nonrepro
logcpm_nonrepro_females <- logcpm_nonrepro[, 4:9]
heatmap_all_nonrepro_50 <- pheatmap(logcpm_nonrepro_females[rownames(head(qlf.androvsgyno_nonrepro_FDR_sig, 50)), ],
                                    cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T,
                                    main = "Top 50 significant DEGs between female morphs in nonreproductive tissue")
 
### Volcano plots (which requires differential gene analysis) ###
# compares androgenics (as the baseline) to gynogenics 
summary(decideTests(qlf.androvsgyno_nonrepro))
plotMD(qlf.androvsgyno_nonrepro)
abline(h=c(-1,1), col="blue") 


qlf.androvsgyno_nonrepro_FDR_df <- as.data.frame(qlf.androvsgyno_nonrepro_FDR)

qlf.androvsgyno_nonrepro_FDR_df <- qlf.androvsgyno_nonrepro_FDR_df %>% 
  mutate(
    Expression = case_when(FDR < 0.05 & logFC <= 0 ~ "Androgenic-biased",
                           FDR < 0.05 & logFC >= 0 ~ "Gynogenic-biased",
                           TRUE ~ "Unchanged")
  )
table(qlf.androvsgyno_nonrepro_FDR_df$Expression)

nonrepro_volcano <- ggplot(qlf.androvsgyno_nonrepro_FDR_df, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 2, alpha = 0.6) +
  scale_color_manual(values = c("chartreuse3", "purple", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  geom_hline(yintercept = -log(0.05,10), colour = "gray") +
  # geom_vline(xintercept = c(1, -1), colour="gray") +
  ggtitle("F - Somatic non-reproductive") + theme_bw() +
  theme(plot.title.position = "plot",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 22, hjust = 0.05, vjust = 2, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))
nonrepro_volcano

nonrepro_volcano_legend <- ggplot(qlf.androvsgyno_nonrepro_FDR_df, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 1.5) +
  scale_color_manual(values = c("chartreuse3", "purple", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  geom_hline(yintercept = -log(0.05,10), colour = "gray") +
  # geom_vline(xintercept = c(1, -1), colour="gray") +
  ggtitle("F - Somatic non-reproductive") + theme_bw() +
  theme(plot.title.position = "plot",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 22, hjust = -0.1, vjust = 2, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
    legend.position = "bottom", legend.title = element_text(size = 18), legend.text = element_text(size = 18, margin = margin(t = 0, b = 0)))
```
\
I repeat this for germline and somatic reproductive tissue. Finally, I put all the plots together with ggarrange for Figure 5.
```
pca_legend <- ggpubr::get_legend(non_repro_pca_legend) 
pca_legend_gg <- as_ggplot(pca_legend) + theme(plot.margin = margin(t = -40, r = 0, b = 0, l = 0))
vol_legend <- ggpubr::get_legend(nonrepro_volcano_legend)
png(file="./plots/gene_exp_div.png", height = 850, width = 1000)
ggarrange(ggarrange(germline_pca,repro_pca, non_repro_pca, ncol = 3, common.legend = FALSE), pca_legend_gg, ggarrange(germline_volcano, repro_volcano, nonrepro_volcano, legend = "none", ncol = 3), vol_legend, nrow = 4, heights = c(1, 0.05, 0.9, 0.1))
dev.off()
```