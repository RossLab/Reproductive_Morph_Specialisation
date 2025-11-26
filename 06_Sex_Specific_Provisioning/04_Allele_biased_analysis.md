Necessary packages:
```
library(dplyr) # v1.1.4
library(GenomicRanges) # v1.54.4
library(tidyr) # v1.3.1
library(stringr) # v1.5.1
library(lme4) # v1.1-35.1
library(lmerTest) # v3.1-3
library(emmeans) # v1.10.0
```

Allele-biased deposition of X'X transcripts within the eggs of gynogenic females (ake eggs determined to become females) was done in the same way as in 03_Allele_Biased_Expression

```
## Files I need: 
# RNA seq design 
design_ase <- read.csv("r_input/maternal_deposit_design.csv")
design_ase <- design_ase[1:6,]
design_ase$path <- paste0("r_input/", design_ase$sample_id, ".ase.featureCounts.txt",sep="")
design_ase


# All genes 
genes_and_chromosomes <- read.table("r_input/Bcop_v3_genes_and_chromosomes.tsv")
X_Inv <- genes_and_chromosomes[(genes_and_chromosomes$V1 == "X") | (genes_and_chromosomes$V1 == "Inversion"), ]
nrow(X_Inv) # 7371 (3898 X genes and 3473 Inv genes)

# Which X and X' genes are orthologs, output of orthofinder. 
X_Inv_orthogroups <- read.table("r_input/X_vs_Inv_HOG.tsv", header = T, sep="\t")
X_Inv_SCO <- read.table("r_input/X_vs_Inv_Orthogroups_SCOs.tsv")

X_Inv_SCO_genes <- X_Inv_orthogroups[X_Inv_orthogroups$HOG %in% X_Inv_SCO$V1, ]
nrow(X_Inv_SCO_genes) # 2152
X_Inv_SCO_genes <- X_Inv_SCO_genes[, c("Bcop_Inv_genes", "Bcop_X_genes")]
colnames(X_Inv_SCO_genes) <- c("Inv_transcript", "X_transcript")
X_Inv_SCO_genes[] <- lapply(X_Inv_SCO_genes, function(x) sub("\\.t.*", "", x)) # Removes the transcript number
X_Inv_SCO_genes

#### I need to normalise the gene expression by length. Because the X and X' genes are homologs I can't assume that they are the same length
# Get gene lengths for each gene 
# Build a TxDb object from  GTF
txdb <- GenomicFeatures::makeTxDbFromGFF("r_input/Bcop_v3_geneid.augustus.gtf")
# Get exons by gene
exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
# Calculate total non-overlapping exon length per gene
gene_lengths <- sum(width(reduce(exons_by_gene)))
length(gene_lengths) # 24050

X_Inv_lengths <- gene_lengths[names(gene_lengths) %in% X_Inv$V2]
length(X_Inv_lengths) # 7371


### Read in the featureCounts counts, each column is one sample
ase_featureCounts <- lapply(design_ase$path, function(x)read.table(x, header=T)) # applies the read.table() function to each file path, reading each file into a data frame
names(ase_featureCounts) <- design_ase$sample_id # Names each dataframe in the list
ase_featureCounts <- lapply(ase_featureCounts, function(x) x[, c(1, 7)]) # For each data frame in the list, only columns 1 and 7 are kept.
ase_featureCounts <- lapply(names(ase_featureCounts), function(name) {
  x <- ase_featureCounts[[name]]
  colnames(x)[2] <- name
  x 
}) # Assigns the corresponding sample ID (name) as the name of the second column (the count data column)
names(ase_featureCounts) <- design_ase$sample_id # Names each dataframe in the list because the name got lost

ase_raw_gene_counts <- data.frame("Geneid" = X_Inv$V2) # Subset for only the genes on the X or Inversion
for (name in names(ase_featureCounts)) {
  # Extract the current table
  x <- ase_featureCounts[[name]]
  # Extract column 2, setting it as a new column in `autosomal_counts` with the sample ID as the column name
  ase_raw_gene_counts <- merge(ase_raw_gene_counts, x, by = "Geneid", all.x = TRUE)
}
rownames(ase_raw_gene_counts) <- ase_raw_gene_counts$Geneid
ase_raw_gene_counts <- ase_raw_gene_counts[ , 2:7]
nrow(ase_raw_gene_counts) # 7371


### Now I have the gene length of each gene, and the RNA read count, I can normalise the counts by length
ase_raw_gene_counts$Lengths <- X_Inv_lengths[match(rownames(ase_raw_gene_counts), names(X_Inv_lengths))]
head(ase_raw_gene_counts)
# Store the lengths
lengths_kb <- ase_raw_gene_counts$Lengths / 1000
# Compute reads per kilobase for all samples
rpk <- sweep(ase_raw_gene_counts[, -ncol(ase_raw_gene_counts)], 1, lengths_kb, FUN = "/")
# Compute the scaling factors (per column)
scaling_factors <- colSums(rpk) / 1e6
# Compute transcripts per million
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/")
# Add gene names back
rownames(tpm) <- rownames(ase_raw_gene_counts)
head(tpm)

### Put all the length normalised gene counts in a table 
X_Inv_SCO_genes_count <- data.frame(Number = c(1:nrow(X_Inv_SCO_genes)), X_Inv_SCO_genes)
# Adds a number ordering the table. This is to recognise the orthogroups later and also so merge doesn't change the order of the genes
head(X_Inv_SCO_genes_count)

tpm_X <- tpm
colnames(tpm_X) <- c(paste0(design_ase$sample_id, "_X",sep=""))
tpm_X$X_transcript <- rownames(tpm_X)
tpm_Inv <- tpm
colnames(tpm_Inv)<- c(paste0(design_ase$sample_id, "_Inv",sep=""))
tpm_Inv$Inv_transcript <- rownames(tpm_Inv)

X_Inv_SCO_genes_count
X_Inv_SCO_genes_count <- merge(X_Inv_SCO_genes_count, tpm_X, by = "X_transcript")
X_Inv_SCO_genes_count <- merge(X_Inv_SCO_genes_count, tpm_Inv, by = "Inv_transcript")
X_Inv_SCO_genes_count <- X_Inv_SCO_genes_count[order(X_Inv_SCO_genes_count$Number),]

# write.csv(X_Inv_SCO_genes_count, "outputs/allele_specific_gene_count_maternal_deposit.csv", row.names = FALSE)

######################## Allele-specific expression ######################## 
#### Now that I have all the length normalised gene counts for X genes and their ortholog on the Inv, I can carry out stats to see if they are differentially expressed
######## 0-4 hrs ########
design_ase
# X_Inv_SCO_genes_count <- read.csv("output/allele_specific_gene_count_all_tissues.csv")
gene_count_04hr <- X_Inv_SCO_genes_count[, c(1:3, 4:6, 10:12)]
nrow(gene_count_04hr) # 2152, same number as single orthologs between X and Inv


gene_count_04hr_long <- gene_count_04hr %>%
  pivot_longer(
    cols = starts_with("female"),
    names_to = "Sample_Group",
    values_to = "Expression"
  ) %>%
  mutate(
    Sample = str_extract(Sample_Group, "rep\\d+"),
    Chrom = ifelse(str_detect(Sample_Group, "_X$"), "X", "Inv")
  ) %>%
  select(Number, Sample, Chrom, Expression, Inv_transcript, X_transcript)

gene_count_04hr_long$Number <- as.factor(gene_count_04hr_long$Number)
gene_count_04hr_long$logTPM <- log1p(gene_count_04hr_long$Expression)



model <- lmer(logTPM ~ -1 + Chrom:Number + (1 | Sample), data = gene_count_04hr_long)
model_summary <- summary(model)
# This estimates the average expression for each homolog on different genes. Significance in this case means that they are significantly different from 0 

# Now I am extracting the estimated mean for each homolog on each chromosome, to do pairwise comparisons
emm <- emmeans(model, ~ Chrom:Number)
contrast <- contrast(emm, method = "pairwise", by = "Number") |> summary(infer = TRUE, adjust = "BH")

contrast_summary <- contrast |> 
  as.data.frame() |> 
  filter(str_detect(contrast, "Inv - X")) |> 
  mutate(Number = str_extract(as.character(Number), "\\d+"),
         logFC = estimate,
         pval = p.value,
         sig = case_when(
           pval < 0.001 ~ "***",
           pval < 0.01  ~ "**",
           pval < 0.05  ~ "*",
           TRUE         ~ "ns")) |> 
  select(Number, logFC, pval, sig)

X_Inv_SCO_genes_count_summary <- X_Inv_SCO_genes_count[, c(1:3)]

allele_specific_summary <- merge(X_Inv_SCO_genes_count_summary, contrast_summary, by = "Number")
allele_specific_summary[order(allele_specific_summary$logFC, decreasing = T),]

# write.csv(allele_specific_summary, "outputs/allele_specific_exp_summary_04hrs.csv", row.names = FALSE)


allele_specific_summary_readin <- read.csv("outputs/allele_specific_exp_summary_04hrs.csv")
#OR
# allele_specific_summary_readin <- allele_specific_summary
dNdS_readin <- read.csv("C:/Users/s2556496/Desktop/All/Morph_Specialisation_Gene_Expression/B_coprophila_morph_gene_divergence_FINAL/04_dnds/output/X_Inv_ortholog_dNdS.csv")

nrow(allele_specific_summary_readin) # 2152
nrow(allele_specific_summary_readin[allele_specific_summary_readin$pval < 0.05, ]) # 1097

allele_specific_summary_sig <- allele_specific_summary_readin |>
  mutate(
    bias = case_when(
      logFC > 0 & pval <= 0.05 ~ "Inv_biased",
      logFC < 0 & pval <= 0.05 ~ "X_biased",
      TRUE ~ "Unbiased"
    )
  )
table(allele_specific_summary_sig$bias)

## MA plot

mean_expr <- gene_count_04hr_long %>%
  group_by(Number, Chrom) %>%
  summarise(mean_TPM = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Chrom, values_from = mean_TPM) %>%
  mutate(Number = as.integer(as.character(Number)))

allele_specific_summary_MA <- allele_specific_summary_sig %>%
  left_join(mean_expr, by = "Number") %>%
  mutate(
    A = log2((Inv + X) / 2 + 1),  # average expression, +1 to avoid log(0)
    M = logFC                     # same as log2(Inv) - log2(X)
  )

MAplot_ASE_04hr <- ggplot(allele_specific_summary_MA, aes(x = A, y = M, color = bias)) +
  geom_point(size = 1.5) +
  scale_color_manual(labels = c("Inversion allele", "Unbiased", "X allele"), values = c("purple", "gray50", "chartreuse3")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Average expression (log2 TPM)",
    y = "Allele-specific expression (log2 Inv - X)"
  ) +
  labs(title = "C - Allele-specific expression") + theme_bw() +
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
  guides(color = guide_legend(override.aes = list(size = 5))) + labs(col = "Expression bias")
MAplot_ASE_04hr

### Adding dNdS information to expression 
exp_dNdS_summary <- merge(allele_specific_summary_sig, dNdS_readin, by = "X_transcript")
exp_dNdS_summary <- exp_dNdS_summary |>
  select(X_transcript, Inv_transcript.x, logFC, pval, bias, Inv_dNdS, X_dNdS, dNdS_diff)

## There are two annoying out of bound points
# Define plot limits
xlim_vals <- c(0, 1.25)
ylim_vals <- c(0, 1.25)
# Subset the points that are out of bounds
out_of_bounds <- exp_dNdS_summary |>
  filter(X_dNdS > xlim_vals[2] | Inv_dNdS > ylim_vals[2])

set.seed(123)  # Set seed once
out_of_bounds$jittered_X <- pmin(out_of_bounds$X_dNdS, xlim_vals[2]) + runif(nrow(out_of_bounds), -0.05, 0.05)
out_of_bounds$jittered_Y <- pmin(out_of_bounds$Inv_dNdS, ylim_vals[2])

exp_dNdS_plot_04hr <- ggplot(data= exp_dNdS_summary, aes(x = X_dNdS, y = Inv_dNdS, col = bias)) +
  geom_point(alpha = 0.7) +
  geom_point(data = out_of_bounds,
             aes(x = jittered_X, y = jittered_Y),
             shape = 2,
             size = 2,
             inherit.aes = TRUE,
             show.legend = FALSE) +
  coord_cartesian(xlim=xlim_vals, ylim=ylim_vals) +
  geom_abline(intercept = 0, slope = 1, col = "grey") + 
  scale_color_manual("Expression",values=c("purple","grey", "chartreuse2")) +
  ggtitle("D - Rate of evolution of X'X allele") + theme_bw() +
  xlab("dNdS of X allele") + ylab("dNdS of Inversion allele") +
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
  )
```
\
\
All plots were arranged together for Figure 6.
```
png(file="./outputs/maternal_deposit_new.png", height = 900, width = 800)
ggarrange(ggarrange(pca_04, volcano_04), ggarrange(MAplot_ASE_04hr, exp_dNdS_plot_04hr, align = "hv", common.legend = TRUE, legend = "bottom"), nrow = 2)
dev.off()
```
