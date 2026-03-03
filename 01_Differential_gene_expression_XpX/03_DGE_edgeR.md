# Differential gene expression analysis of X-linked genes, in edge R #
Differential gene expression (DGE) analysis was carried out in R (v4.3.2), using the edgeR package.The analysis was done separate for each tissue type. 
\
The necessary packages are: 
```
library(edgeR) #v4.0.9
library(PCAtools) #v2.14.0
library(pheatmap) #v1.0.12
library(tidyverse) #v2.0.0
library(ggpubr) #v0.6.0
library(dplyr) #v1.1.4
library(tidyr) #v1.3.1
library(purrr) #v1.0.2
library(data.table) #v1.14.10

```
I first read in the gametolog pairs which I obtained from Orthofinder. 
```
X_Inv_orthogroups <- read.table("inputs/X_vs_Inv_longest_OG_HOG.tsv", header = T, sep="\t")
X_Inv_SCO <- read.table("inputs/X_vs_Inv_longest_Orthogroups_SCOs.tsv")
head(X_Inv_orthogroups)
nrow(X_Inv_orthogroups) #2875
head(X_Inv_SCO)
nrow(X_Inv_SCO) #2703

X_Inv_SCO_genes <- X_Inv_orthogroups[X_Inv_orthogroups$HOG %in% X_Inv_SCO$V1, ]
nrow(X_Inv_SCO_genes) # 2703
X_Inv_SCO_genes <- X_Inv_SCO_genes[, c("Bcop_Inv_genes", "Bcop_X_genes")]
X_Inv_SCO_genes[] <- lapply(
  X_Inv_SCO_genes,
  function(x) sub("\\.t[0-9]+$", "", x)
)
# there is one annoying gene name that does have .t1 in it, so I'm adding it back
X_Inv_SCO_genes[which(X_Inv_SCO_genes$Bcop_X_genes == "g13114"), 2] <- "g13114.t1"
```
\
A design matrix is needed (for example, the table in 00_Data.md) to specify which group (sex/tissue) the sample belongs to. For the purpose of reading in the featureCount files that are in separate folders, I am also adding the file path to the design matrix. 
I also need a gene-chromosome file to locate which gene maps onto which chromosome, so that I can subset for X-linked genes only.
```
## Design matrix ##
design_XpX <- read.csv("inputs/gene_expression_divergence_gyno_andro_male_design.csv")
design_XpX$path_XpX <- paste0("inputs/XpX_mapping/", design_XpX$sample_id, ".txt",sep="")
group <- factor(paste(design_XpX$sex,design_XpX$tissue,sep="."))
design_XpX <- cbind(design_XpX,group=group)
design_XpX

## Genes and the chromosome they belong to ##
genes_and_chromosomes <- read.table("inputs/Bcop_v3_genes_and_chromosomes.tsv")
```
\
I wrangle the data into a format where I can sum the expression read counts of the X and X' gametologs, per sample. 

```
## Make a long featureCounts dataframe where one sample is put after the next 
featureCounts_long <- map2_dfr(
  design_XpX$path_XpX,
  design_XpX$sample_id,
  ~ read.table(.x, header = TRUE) %>%
    dplyr::select(Geneid, count = 7) %>%   # keep Geneid + counts
    dplyr::mutate(sample_id = .y)
)
## Make the ortholog dataframe longer to facilitate joining, and assign each pair a pair number 
orthologs_long <- X_Inv_SCO_genes %>%
  mutate(ortholog_id = row_number()) %>%
  pivot_longer(
    cols = c(Bcop_X_genes, Bcop_Inv_genes),
    names_to = "chrom",
    values_to = "Geneid"
  ) 
## Left join (to keep all the orthologs ordered and get rid of the genes that are not in the list)
expr_ortho <- left_join(orthologs_long, featureCounts_long, by = "Geneid")
# Sum all the values per group per sample 
expr_summed <- expr_ortho %>%
  group_by(ortholog_id, sample_id) %>%
  summarise(count = sum(count), .groups = "drop")
# Pivot wide 
expr_master <- expr_summed %>%
  pivot_wider(
    names_from  = sample_id,
    values_from = count,
    values_fill = 0
  )
expr_master<- as.data.frame(left_join(expr_master, orthologs_long[orthologs_long$chrom == "Bcop_X_genes",], by = "ortholog_id"))
rownames(expr_master) <- expr_master$Geneid
expr_master <- expr_master[, c(2, 13, 22:28, 3:12, 14:21)]
nrow(expr_master) #2703
head(expr_master)
```
I first do this with all tissue to get an idea of the sample distances, before repeating the comparisons by tissue type. 
```
## Read the counts into a DGEList, specifying the design ##
DGEList_XpX <- DGEList(expr_master, group = design_XpX$group)
nrow(DGEList_XpX$counts) #2703

## Filter lowly expressed genes ##
keep <- filterByExpr(DGEList_XpX)
table(keep) # Keep 2487 genes, get rid of 216 genes 
DGEList_XpX <- DGEList_XpX[keep, , keep.lib.sizes=FALSE]
nrow(DGEList_XpX$counts) #2487

## Normalisation ##
DGEList_XpX <- normLibSizes(DGEList_XpX)
head(DGEList_XpX)
print(DGEList_XpX$samples)

## Visualise sample distances ##
plotMDS(DGEList_XpX)

####################### PCA plots ###########################
counts_for_pca<-cpm(DGEList_XpX,log=TRUE,prior.count=1) 

# Set group with desired order and custom labels
DGEList_XpX$samples$group <- factor(DGEList_XpX$samples$group, levels = c(
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

pca_output <- pca(counts_for_pca, metadata = DGEList_XpX$samples)
screeplot(pca_output)


all_tissue_pca_XpX <- biplot(pca_output, colby = "group", colkey = c(
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
all_tissue_pca_XpX
```
\
Now I compare X-linked gene expression between androgenic and gynogenic females in somatic non-reproductive tissue. 
```
######################################################## SOMATIC NON-REPRODUCTIVE TISSUE (XpX mapping) ########################################################
## Making a DGEList object ##
design_XpX_nonrepro <- design_XpX[c(7:9, 16:18, 25:27), ]
counts_XpX_nonrepro <- expr_master[, c(7:9, 16:18, 25:27)]
DGEList_XpX_nonrepro <- DGEList(counts = counts_XpX_nonrepro, group = design_XpX_nonrepro$sex)
head(DGEList_XpX_nonrepro)
nrow(DGEList_XpX_nonrepro$counts) #2703

# Constructing a model matrix for glm 
design_XpX_nonrepro_matrix <- model.matrix(~0 + design_XpX_nonrepro$sex)
rownames(design_XpX_nonrepro_matrix) <- design_XpX_nonrepro$sample_id
colnames(design_XpX_nonrepro_matrix) <- sort(unique(design_XpX_nonrepro$sex))
design_XpX_nonrepro_matrix

## Filter lowly expressed genes ##
keep <- filterByExpr(DGEList_XpX_nonrepro)
table(keep) # Keep 2336 genes, get rid of 367 genes 
DGEList_XpX_nonrepro <- DGEList_XpX_nonrepro[keep, , keep.lib.sizes=FALSE]
nrow(DGEList_XpX_nonrepro$counts) #2336

## Normalisation ##
DGEList_XpX_nonrepro <- normLibSizes(DGEList_XpX_nonrepro)
head(DGEList_XpX_nonrepro)

## Visualise sample distances ##
plotMDS(DGEList_XpX_nonrepro)
counts_for_pca<-cpm(DGEList_XpX_nonrepro,log=TRUE,prior.count=1) 
pca_output <- pca(counts_for_pca, metadata = DGEList_XpX_nonrepro$samples)

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

## Here I am making a separate ggplot just for extracting the legend in downstream analysis. I only need to do this once and not repeat it with the other tissues.
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
DGEList_XpX_nonrepro <- estimateDisp(DGEList_XpX_nonrepro, design = design_XpX_nonrepro_matrix)
# This is a plot of the NB dispersion for the dataset 
plotBCV(DGEList_XpX_nonrepro)

# Given raw counts, NB dispersion(s) and a design matrix, glmQLFit() fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
fit_XpX_nonrepro <- glmQLFit(DGEList_XpX_nonrepro, design_XpX_nonrepro_matrix)
head(fit_XpX_nonrepro$coefficients)
# compares androgenics (as the baseline) to gynogenics 
qlf.androvsgyno_XpX_nonrepro <- glmQLFTest(fit_XpX_nonrepro, contrast=c(-1, 1, 0))
qlf.androvsgyno_XpX_nonrepro_FDR <- topTags(qlf.androvsgyno_XpX_nonrepro, 2336)
# write.csv(qlf.androvsgyno_XpX_nonrepro_FDR, "outputs/XpX_nonrepro_androvsgyno_all.csv")
qlf.androvsgyno_XpX_nonrepro_FDR_sig <- qlf.androvsgyno_XpX_nonrepro_FDR$table[qlf.androvsgyno_XpX_nonrepro_FDR$table$FDR < 0.05,]
summary(decideTests(qlf.androvsgyno_XpX_nonrepro))
# write.csv(qlf.androvsgyno_XpX_nonrepro_FDR_sig, "outputs/XpX_nonrepro_androvsgyno_sig.csv")

## Plotting heatmaps etc 
logcpm_XpX_nonrepro <- cpm(DGEList_XpX_nonrepro, log=TRUE)
logcpm_XpX_nonrepro
logcpm_XpX_nonrepro_females <- logcpm_XpX_nonrepro[, 4:9]
heatmap_XpX_nonrepro_50 <- pheatmap(logcpm_XpX_nonrepro_females[rownames(head(qlf.androvsgyno_XpX_nonrepro_FDR_sig, 50)), ],
                                    cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = T,
                                    main = "Top 50 significant DEGs between female morphs in nonreproductive tissue (XpX mapping)")

### Volcano plots (which requires differential gene analysis) ###
# compares androgenics (as the baseline) to gynogenics 
summary(decideTests(qlf.androvsgyno_XpX_nonrepro))
plotMD(qlf.androvsgyno_XpX_nonrepro)
abline(h=c(-1,1), col="blue") 


qlf.androvsgyno_XpX_nonrepro_FDR_df <- as.data.frame(qlf.androvsgyno_XpX_nonrepro_FDR)

qlf.androvsgyno_XpX_nonrepro_FDR_df <- qlf.androvsgyno_XpX_nonrepro_FDR_df %>% 
  mutate(
    Expression = case_when(FDR < 0.05 & logFC <= 0 ~ "Androgenic-biased",
                           FDR < 0.05 & logFC >= 0 ~ "Gynogenic-biased",
                           TRUE ~ "Unchanged")
  )
table(qlf.androvsgyno_XpX_nonrepro_FDR_df$Expression)

nonrepro_volcano <- ggplot(qlf.androvsgyno_XpX_nonrepro_FDR_df, aes(logFC, -log(FDR,10))) +
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

## As above with the PCA, here I am making a separate ggplot just for extracting the legend in downstream analysis. I only need to do this once and not repeat it with the other tissues.
nonrepro_volcano_legend <- ggplot(qlf.androvsgyno_XpX_nonrepro_FDR_df, aes(logFC, -log(FDR,10))) +
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
        plot.title = element_text(size = 22, hjust = -0.1, vjust = 2, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.position = "bottom", legend.title = element_text(size = 18), legend.text = element_text(size = 18, margin = margin(t = 0, b = 0)))
        
```
\
I then repeat the same process with germline and somatic reproductive tissue. Finally I use ggarrange to put all the plots together for Figure 3.
```
pca_legend <- ggpubr::get_legend(non_repro_pca_legend) 
pca_legend_gg <- as_ggplot(pca_legend) + theme(plot.margin = margin(t = -40, r = 0, b = 0, l = 0))
vol_legend <- ggpubr::get_legend(nonrepro_volcano_legend)
png(file="./outputs/gene_exp_div_XpX_mapping.png", height = 850, width = 1000)
ggarrange(ggarrange(germline_pca,repro_pca, non_repro_pca, ncol = 3, common.legend = FALSE), pca_legend_gg, ggarrange(germline_volcano, repro_volcano, nonrepro_volcano, legend = "none", ncol = 3), vol_legend, nrow = 4, heights = c(1, 0.05, 0.9, 0.1))
dev.off()
```