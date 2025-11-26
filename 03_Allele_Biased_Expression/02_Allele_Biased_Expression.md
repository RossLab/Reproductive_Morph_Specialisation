The big difference between conventional allele-specific expression and my analysis, is that I cannot assume the gametologs on the X and X' are the same length (whether due to degeneration or incomplete assembly on the X'). Therefore I have to length-normalise before being able to carry out any expression comparisons between the gametologs. 

This is done in R (v4.3.2).

Required packages:
```
Packages:
library(tidyverse) # v2.0.0
library(ggpubr) # v0.6.0
library(dplyr) # v1.1.4
library(GenomicRanges) # v1.54.1
library(tidyr) # v1.3.1
library(stringr) # v1.5.1
library(lme4) # v1.1-35.1
library(lmerTest) #v1.1-3
library(emmeans) #v1.10.0
```
```
## Files I need: 
# RNA seq design 
design <- read.csv("r_input/allele_specific_expression_gyno_design.csv", header = T)

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

# Read in featureCount files
design$path <- paste0("r_input/", design$sample_id, "_featureCounts.txt",sep="")
design

#### Read count normalisation according to sequence length of homolog 
# Get gene lengths for each gene 
txdb <- GenomicFeatures::makeTxDbFromGFF("r_input/Bcop_v3_geneid.augustus.gtf")
# Get exons by gene
exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
# Calculate total non-overlapping exon length per gene
gene_lengths <- sum(width(reduce(exons_by_gene)))
length(gene_lengths) # 24050 genes, which is the correct number 

X_Inv_lengths <- gene_lengths[names(gene_lengths) %in% X_Inv$V2]
length(X_Inv_lengths) # 7371 genes, again the correct number 

### Read in the featureCounts counts, each column is one sample
featureCounts <- lapply(design$path, function(x)read.table(x, header=T)) # applies the read.table() function to each file path, reading each file into a data frame
names(featureCounts) <- design$sample_id # Names each dataframe in the list
featureCounts <- lapply(featureCounts, function(x) x[, c(1, 7)]) # For each data frame in the list, only columns 1 and 7 are kept.
featureCounts <- lapply(names(featureCounts), function(name) {
  x <- featureCounts[[name]]
  colnames(x)[2] <- name
  x 
}) # Assigns the corresponding sample ID (name) as the name of the second column (the count data column)
names(featureCounts) <- design$sample_id # Names each dataframe in the list because the name got lost

raw_gene_counts <- data.frame("Geneid" = X_Inv$V2) # Subset for only the genes on the X or Inversion
for (name in names(featureCounts)) {
  # Extract the current table
  x <- featureCounts[[name]]
  # Extract column 2, setting it as a new column in `autosomal_counts` with the sample ID as the column name
  raw_gene_counts <- merge(raw_gene_counts, x, by = "Geneid", all.x = TRUE)
}
rownames(raw_gene_counts) <- raw_gene_counts$Geneid
raw_gene_counts <- raw_gene_counts[ , 2:10]
nrow(raw_gene_counts) # 7371

# Now I have the gene length of each gene, and the RNA read count, I can normalise the counts by length. I am transforming the values in transcript per million (tpm) which takes into account gene length.
raw_gene_counts$Lengths <- X_Inv_lengths[match(rownames(raw_gene_counts), names(X_Inv_lengths))]
head(raw_gene_counts)
# Store the lengths
lengths_kb <- raw_gene_counts$Lengths / 1000
# Compute reads per kilobase for all samples
rpk <- sweep(raw_gene_counts[, -ncol(raw_gene_counts)], 1, lengths_kb, FUN = "/")
# Compute the scaling factors (per column)
scaling_factors <- colSums(rpk) / 1e6
# Compute transcripts per million
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/")
# Add gene names back
rownames(tpm) <- rownames(raw_gene_counts)

### Put all the length normalised gene counts in a table 
X_Inv_SCO_genes_count <- data.frame(Number = c(1:nrow(X_Inv_SCO_genes)), X_Inv_SCO_genes)
# Adds a number ordering the table. This is to recognise the orthogroups later and also so merge doesn't change the order of the genes
head(X_Inv_SCO_genes_count)

tpm_X <- tpm
colnames(tpm_X) <- c(paste0(design$sample_id, "_X",sep=""))
tpm_X$X_transcript <- rownames(tpm_X)
tpm_Inv <- tpm
colnames(tpm_Inv)<- c(paste0(design$sample_id, "_Inv",sep=""))
tpm_Inv$Inv_transcript <- rownames(tpm_Inv)

X_Inv_SCO_genes_count
X_Inv_SCO_genes_count <- merge(X_Inv_SCO_genes_count, tpm_X, by = "X_transcript")
X_Inv_SCO_genes_count <- merge(X_Inv_SCO_genes_count, tpm_Inv, by = "Inv_transcript")
X_Inv_SCO_genes_count <- X_Inv_SCO_genes_count[order(X_Inv_SCO_genes_count$Number),]

write.csv(X_Inv_SCO_genes_count, "output/allele_specific_gene_count_all_tissues.csv", row.names = FALSE)
# Throughout the script I write some of the tables into csvs so when I restart an R session I don't have to rerun the code, but just read in the csv.
```
\
Now that I have all the length-normalised read counts for X and X' homologs, I use a linear model to assess differential expression between pairs of homologs. As before, this is done by tissue. 
```
#### Somatic non-reproductive tissue ####
X_Inv_SCO_genes_count <- read.csv("output/allele_specific_gene_count_all_tissues.csv")
gene_count_somatic_nonrepro <- X_Inv_SCO_genes_count[, c(1:3, 10:12, 19:21)]
nrow(gene_count_somatic_nonrepro) # 2152, same number as single orthologs between X and Inv, which is good 

gene_count_somatic_nonrepro_long <- gene_count_somatic_nonrepro %>%
  pivot_longer(
    cols = starts_with("B"),
    names_to = "Sample_Group",
    values_to = "Expression"
  ) %>%
  mutate(
    Sample = str_extract(Sample_Group, "B\\d+"),
    Chrom = ifelse(str_detect(Sample_Group, "_X$"), "X", "Inv")
  ) %>%
  select(Number, Sample, Chrom, Expression, Inv_transcript, X_transcript)

gene_count_somatic_nonrepro_long$Number <- as.factor(gene_count_somatic_nonrepro_long$Number)
gene_count_somatic_nonrepro_long$logTPM <- log1p(gene_count_somatic_nonrepro_long$Expression)

model <- lmer(logTPM ~ -1 + Chrom:Number + (1 | Sample), data = gene_count_somatic_nonrepro_long)
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

write.csv(allele_specific_summary, "output/allele_specific_exp_summary_nonrepro.csv", row.names = FALSE)
# The model took ages to run so I'm writing the output into a csv instead of having to run it every time I restart an R session.
```
\
Now that I have the statistical comparison of differential expression between X and X' homologs in this tissue, I can combine that with dNdS information. 
```
allele_specific_summary_readin <- read.csv("output/allele_specific_exp_summary_nonrepro.csv")
dNdS_readin <- read.csv("C:/Users/s2556496/Desktop/All/Morph_Specialisation_Gene_Expression/B_coprophila_morph_gene_divergence_FINAL/04_dnds/output/X_Inv_ortholog_dNdS.csv")

nrow(allele_specific_summary_readin) # 2152
nrow(allele_specific_summary_readin[allele_specific_summary_readin$pval < 0.05, ]) # 1186

allele_specific_summary_sig <- allele_specific_summary_readin |>
  mutate(
    bias = case_when(
      logFC > 0 & pval <= 0.05 ~ "Inv_biased",
      logFC < 0 & pval <= 0.05 ~ "X_biased",
      TRUE ~ "Unbiased"
    )
  )
table(allele_specific_summary_sig$bias)
```
\
To visualise the data, I am first making an MA plot to show the differential expression between homologs. I am then comparing their dNdS values, colour coded by whether there is biased expression towards one homolog.
```
mean_expr <- gene_count_somatic_nonrepro_long %>%
  group_by(Number, Chrom) %>%
  summarise(mean_TPM = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Chrom, values_from = mean_TPM)

allele_specific_summary_MA <- allele_specific_summary_sig %>%
  left_join(mean_expr, by = "Number") %>%
  mutate(
    A = log2((Inv + X) / 2 + 1),  # average expression, +1 to avoid log(0)
    M = logFC                     # same as log2(Inv) - log2(X)
  )

MAplot_ASE_nonrepro <- ggplot(allele_specific_summary_MA, aes(x = A, y = M, color = bias)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("purple", "gray50", "chartreuse3")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Average expression (log2 TPM)",
    y = "Allele-specific expression (log2 Inv - X)") +
  ggtitle("C - Somatic non-reproductive") + theme_bw() +
  theme(plot.title.position = "plot",
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

### dNdS comparison between homologs, with expression bias in non-reproductive tissue. 
# Adding dNdS information to expression 
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

exp_dNdS_plot_nonrepro <- ggplot(data= exp_dNdS_summary, aes(x = X_dNdS, y = Inv_dNdS, col = bias)) +
  geom_point(alpha = 0.7) +
  geom_point(data = out_of_bounds,
             aes(x = jittered_X, y = jittered_Y),
             shape = 2,
             size = 2,
             inherit.aes = TRUE,
             show.legend = FALSE) +
  coord_cartesian(xlim=xlim_vals, ylim=ylim_vals) +
  geom_abline(intercept = 0, slope = 1, col = "grey") + 
  scale_color_manual("Expression bias",labels = c("Inversion allele", "Unbiased", "X allele"), values=c("purple","grey", "chartreuse2")) +
  ggtitle("F - Somatic non-reproductive") + theme_bw() +
  xlab("dNdS of X allele") + ylab("dNdS of Inversion allele") +
  theme(plot.title.position = "plot",
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")
```
\
Repeating this with somatic non-reproductive and germline tissue.
```
#### Somatic reproductive tissue ####
design
gene_count_somatic_repro <- X_Inv_SCO_genes_count[, c(1:3, 7:9, 16:18)]
nrow(gene_count_somatic_repro) # 2152, same number as single orthologs between X and Inv

gene_count_somatic_repro_long <- gene_count_somatic_repro %>%
  pivot_longer(
    cols = starts_with("B"),
    names_to = "Sample_Group",
    values_to = "Expression"
  ) %>%
  mutate(
    Sample = str_extract(Sample_Group, "B\\d+"),
    Chrom = ifelse(str_detect(Sample_Group, "_X$"), "X", "Inv")
  ) %>%
  select(Number, Sample, Chrom, Expression, Inv_transcript, X_transcript)

gene_count_somatic_repro_long$Number <- as.factor(gene_count_somatic_repro_long$Number)
gene_count_somatic_repro_long$logTPM <- log1p(gene_count_somatic_repro_long$Expression)

model <- lmer(logTPM ~ -1 + Chrom:Number + (1 | Sample), data = gene_count_somatic_repro_long)
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
# write.csv(allele_specific_summary, "output/allele_specific_exp_summary_repro.csv", row.names = FALSE)

allele_specific_summary <- read.csv("output/allele_specific_exp_summary_repro.csv")
nrow(allele_specific_summary) # 2152
nrow(allele_specific_summary[allele_specific_summary$pval < 0.05, ]) # 765

allele_specific_summary_sig <- allele_specific_summary |>
  mutate(
    bias = case_when(
      logFC > 0 & pval <= 0.05 ~ "Inv_biased",
      logFC < 0 & pval <= 0.05 ~ "X_biased",
      TRUE ~ "Unbiased"
    )
  )

table(allele_specific_summary_sig$bias)

#### MA plot 
mean_expr <- gene_count_somatic_repro_long %>%
  group_by(Number, Chrom) %>%
  summarise(mean_TPM = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Chrom, values_from = mean_TPM)

allele_specific_summary_MA <- allele_specific_summary_sig %>%
  left_join(mean_expr, by = "Number") %>%
  mutate(
    A = log2((Inv + X) / 2 + 1),  # average expression, +1 to avoid log(0)
    M = logFC                     # same as log2(Inv) - log2(X)
  )
MAplot_ASE_repro <- ggplot(allele_specific_summary_MA, aes(x = A, y = M, color = bias)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("purple", "gray50", "chartreuse3")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Average expression (log2 TPM)",
    y = "Allele-specific expression (log2 Inv - X)") +
  ggtitle("B - Somatic reproductive") + theme_bw() +
  theme(plot.title.position = "plot",
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

exp_dNdS_summary <- merge(allele_specific_summary_sig, dNdS_readin, by = "X_transcript")
exp_dNdS_summary <- exp_dNdS_summary |>
  select(X_transcript, Inv_transcript.x, logFC, pval, bias, Inv_dNdS, X_dNdS, dNdS_diff)

# Define plot limits
xlim_vals <- c(0, 1.25)
ylim_vals <- c(0, 1.25)
# Subset the points that are out of bounds
out_of_bounds <- exp_dNdS_summary |>
  filter(X_dNdS > xlim_vals[2] | Inv_dNdS > ylim_vals[2])

exp_dNdS_plot_repro <- ggplot(data= exp_dNdS_summary, aes(x = X_dNdS, y = Inv_dNdS, col = bias)) +
  geom_point(alpha = 0.7) +
  geom_point(data = out_of_bounds,
             aes(x = jittered_X, y = jittered_Y),
             shape = 2,
             size = 2,
             inherit.aes = TRUE,
             show.legend = FALSE)+
  coord_cartesian(xlim=xlim_vals, ylim=ylim_vals) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  scale_color_manual("Expression",values=c("purple","grey", "chartreuse2")) +
  ggtitle("E - Somatic reproductive") + theme_bw() +
  xlab("dNdS of X allele") + ylab("dNdS of Inversion allele") +
  theme(plot.title.position = "plot",
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")


#### Germline tissue ####
design
gene_count_germline <- X_Inv_SCO_genes_count[, c(1:3, 4:6, 13:15)]
nrow(gene_count_germline) # 2152, same number as single orthologs between X and Inv

gene_count_germline_long <- gene_count_germline %>%
  pivot_longer(
    cols = starts_with("B"),
    names_to = "Sample_Group",
    values_to = "Expression"
  ) %>%
  mutate(
    Sample = str_extract(Sample_Group, "B\\d+"),
    Chrom = ifelse(str_detect(Sample_Group, "_X$"), "X", "Inv")
  ) %>%
  select(Number, Sample, Chrom, Expression, Inv_transcript, X_transcript)

gene_count_germline_long$Number <- as.factor(gene_count_germline_long$Number)
gene_count_germline_long$logTPM <- log1p(gene_count_germline_long$Expression)

model <- lmer(logTPM ~ -1 + Chrom:Number + (1 | Sample), data = gene_count_germline_long)
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
write.csv(allele_specific_summary, "output/allele_specific_exp_summary_germline.csv", row.names = FALSE)

allele_specific_summary <- read.csv("output/allele_specific_exp_summary_germline.csv")
nrow(allele_specific_summary) # 2152
nrow(allele_specific_summary[allele_specific_summary$pval < 0.05, ]) # 1078

allele_specific_summary_sig <- allele_specific_summary |>
  mutate(
    bias = case_when(
      logFC > 0 & pval <= 0.05 ~ "Inv_biased",
      logFC < 0 & pval <= 0.05 ~ "X_biased",
      TRUE ~ "Unbiased"
    )
  )

table(allele_specific_summary_sig$bias)

mean_expr <- gene_count_germline_long %>%
  group_by(Number, Chrom) %>%
  summarise(mean_TPM = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = Chrom, values_from = mean_TPM)

allele_specific_summary_MA <- allele_specific_summary_sig %>%
  left_join(mean_expr, by = "Number") %>%
  mutate(
    A = log2((Inv + X) / 2 + 1),  # average expression, +1 to avoid log(0)
    M = logFC                     # same as log2(Inv) - log2(X)
  )

MAplot_ASE_germline <- ggplot(allele_specific_summary_MA, aes(x = A, y = M, color = bias)) +
  geom_point(size = 1.5) +
  scale_color_manual(labels = c("Inversion allele", "Unbiased", "X allele"), values = c("purple", "gray50", "chartreuse3")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Average expression (log2 TPM)",
    y = "Allele-specific expression (log2 Inv - X)") +
  ggtitle("A - Germline") + theme_bw() +
  labs(col = "Expression bias") +
  theme(plot.title.position = "plot",
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16))+
  guides(color = guide_legend(override.aes = list(size = 3)))

exp_dNdS_summary <- merge(allele_specific_summary_sig, dNdS_readin, by = "X_transcript")
exp_dNdS_summary <- exp_dNdS_summary |>
  select(X_transcript, Inv_transcript.x, logFC, pval, bias, Inv_dNdS, X_dNdS, dNdS_diff)

# Define plot limits
xlim_vals <- c(0, 1.25)
ylim_vals <- c(0, 1.25)
# Subset the points that are out of bounds
out_of_bounds <- exp_dNdS_summary |>
  filter(X_dNdS > xlim_vals[2] | Inv_dNdS > ylim_vals[2])

exp_dNdS_plot_germline <- ggplot(data= exp_dNdS_summary, aes(x = X_dNdS, y = Inv_dNdS, col = bias)) +
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
  ggtitle("D - Germline") + theme_bw() +
  xlab("dNdS of X allele") + ylab("dNdS of Inversion allele") +
  theme(plot.title.position = "plot",
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")
```
Finally I use ggarrange to put the plots together for Fig 4.
```
png(file="./output/all_exp_dNdS.png", height = 750, width = 1000)
ggarrange(MAplot_ASE_germline, MAplot_ASE_repro, MAplot_ASE_nonrepro, exp_dNdS_plot_germline, exp_dNdS_plot_repro, exp_dNdS_plot_nonrepro, nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom")
dev.off()
```
