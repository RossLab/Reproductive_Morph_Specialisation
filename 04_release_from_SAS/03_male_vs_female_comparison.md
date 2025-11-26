Now that I have counts for reads of each sex (gynogenic/androgenic female, male) mapped to the N-masked genome. I move to R to look at any differences in the degree of correlation between androgenic female vs male expression, and gynogenic female vs male expression. I first filter away lowly expressed genes and also use edgeR to transform the counts into logCPM values. I do this separately for each tissue.
```
setwd("C:/Users/s2556496/Desktop/All/Morph_Specialisation_Gene_Expression/B_coprophila_morph_gene_divergence_FINAL/06_sexually_antagonistic_selection/")

library(edgeR) # v4.0.9
library(ggplot2) # v3.5.1
library(tidyr) # v1.3.1

## Design matrix ##
design <- read.csv("inputs/gene_expression_divergence_gyno_andro_male_design.csv")
design$path_nmask <- paste0("inputs/N-mask/", design$sample_id, "_nmask.txt",sep="")
group <- factor(paste(design$sex,design$tissue,sep="."))
design <- cbind(design,group=group)
design

## Genes and the chromosome they belong to ##
genes_and_chromosomes <- read.table("inputs/Bcop_v3_genes_and_chromosomes.tsv")
genes_and_chromosomes <- genes_and_chromosomes[genes_and_chromosomes$V1 == "II" | genes_and_chromosomes$V1 == "III" | genes_and_chromosomes$V1 == "IV" | genes_and_chromosomes$V1 == "X", ]
autosomal_genes <- genes_and_chromosomes[genes_and_chromosomes$V1 == "II" | genes_and_chromosomes$V1 == "III" | genes_and_chromosomes$V1 == "IV", ]
X_genes <- genes_and_chromosomes[genes_and_chromosomes$V1 == "X", ]

## Summarising featureCount read files into one dataframe ##
featureCounts_nmask <- lapply(design$path_nmask, function(x)read.table(x, header=T)) # applies the read.table() function to each file path, reading each file into a data frame
names(featureCounts_nmask) <- design$sample_id # Names each dataframe in the list
featureCounts_nmask <- lapply(featureCounts_nmask, function(x) x[, c(1, 7)]) # For each data frame in the list, only columns 1 and 7 are kept.
featureCounts_nmask <- lapply(names(featureCounts_nmask), function(name) {
  x <- featureCounts_nmask[[name]]
  colnames(x)[2] <- name
  x 
}) # Assigns the corresponding sample ID (name) as the name of the second column (the count data column)
names(featureCounts_nmask) <- design$sample_id # Names each dataframe in the list because the name got lost
counts_nmask <- data.frame(GeneID = genes_and_chromosomes$V2)  
colnames(counts_nmask) <- "Geneid"
for (name in names(featureCounts_nmask)) {
  # Extract the current table
  x <- featureCounts_nmask[[name]]
  # Extract column 2, setting it as a new column in `counts_nmask` with the sample ID as the column name
  counts_nmask <- merge(counts_nmask, x, by = "Geneid", all.x = TRUE)
}
rownames(counts_nmask) <- counts_nmask$Geneid
counts_nmask <- counts_nmask[, 2:28]
nrow(counts_nmask) # 20577
head(counts_nmask, n=50)

# Filter away lowly expressed genes 
smallestGroupSize <- 3
keep <- rowSums(counts_nmask >= 10) >= smallestGroupSize
nrow(counts_nmask) # 20577
counts_nmask <- counts_nmask[keep,]
nrow(counts_nmask) # 16082

# Use the edgeR package to produce logcpm values 
counts_nmask_edgeR <- DGEList(counts=counts_nmask)
# Calculate logCPM
counts_nmask_logcpm <- cpm(counts_nmask_edgeR, log=TRUE, prior.count=1)
```
\
Plotting androgenic/gynogenic female expression for autosomal, X'X, or all genes for somatic non-reproductive tissue.
```
##### NON REPRODUCTIVE TISSUE #####
logcpm_nonrepro <- counts_nmask_logcpm[, c(7:9, 16:18, 25:27)]
logcpm_nonrepro_means <- data.frame(
  Male = rowMeans(logcpm_nonrepro[, 1:3]),
  Gynogenic = rowMeans(logcpm_nonrepro[, 4:6]),
  Androgenic = rowMeans(logcpm_nonrepro[, 7:9])
)
logcpm_nonrepro_means

logcpm_nonrepro_means$Gene <- rownames(logcpm_nonrepro_means)
# Reshape the data to long format for easier plotting
logcpm_nonrepro_means_long <- pivot_longer(
  logcpm_nonrepro_means,
  cols = c(Gynogenic, Androgenic),
  names_to = "Group",
  values_to = "Value"
)

logcpm_nonrepro_means_long_autosomes <- logcpm_nonrepro_means_long[logcpm_nonrepro_means_long$Gene %in% autosomal_genes$V2, ]
logcpm_nonrepro_means_long_X <- logcpm_nonrepro_means_long[logcpm_nonrepro_means_long$Gene %in% X_genes$V2, ]

# Plot using ggplot
male_vs_female_nonrepro_all <- ggplot(logcpm_nonrepro_means_long, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE, size = 1) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("G") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

male_vs_female_nonrepro_autosomes <- ggplot(logcpm_nonrepro_means_long_autosomes, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("H") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

male_vs_female_nonrepro_X <- ggplot(logcpm_nonrepro_means_long_X, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("I") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

# Fit a combined linear model with an interaction term
summary(lm(Value ~ Male * Group, data = logcpm_nonrepro_means_long_autosomes))
summary(lm(Value ~ Male * Group, data = logcpm_nonrepro_means_long_X))
# Interaction term is not significant, the slope isn't different 
```
\
Now I do the same for somatic reproductive tissue and germline tissue.
```
##### REPRODUCTIVE TISSUE #####
logcpm_repro <- counts_nmask_logcpm[, c(4:6, 13:15, 22:24)]
logcpm_repro_means <- data.frame(
  Male = rowMeans(logcpm_repro[, 1:3]),
  Gynogenic = rowMeans(logcpm_repro[, 4:6]),
  Androgenic = rowMeans(logcpm_repro[, 7:9])
)
logcpm_repro_means

logcpm_repro_means$Gene <- rownames(logcpm_repro_means)
# Reshape the data to long format for easier plotting
logcpm_repro_means_long <- pivot_longer(
  logcpm_repro_means,
  cols = c(Gynogenic, Androgenic),
  names_to = "Group",
  values_to = "Value"
)

logcpm_repro_means_long_autosomes <- logcpm_repro_means_long[logcpm_repro_means_long$Gene %in% autosomal_genes$V2, ]
logcpm_repro_means_long_X <- logcpm_repro_means_long[logcpm_repro_means_long$Gene %in% X_genes$V2, ]

# Plot using ggplot
male_vs_female_repro_all <- ggplot(logcpm_repro_means_long, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("D") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

male_vs_female_repro_autosomes <- ggplot(logcpm_repro_means_long_autosomes, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("E") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

male_vs_female_repro_X <- ggplot(logcpm_repro_means_long_X, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("F") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

# Fit a combined linear model with an interaction term
summary(lm(Value ~ Male * Group, data = logcpm_repro_means_long_autosomes))
summary(lm(Value ~ Male * Group, data = logcpm_repro_means_long_X))
# Interaction term is not significant, the slope isn't different 

##### GERMLINE TISSUE #####
logcpm_germline <- counts_nmask_logcpm[, c(1:3, 10:12, 19:21)]
logcpm_germline_means <- data.frame(
  Male = rowMeans(logcpm_germline[, 1:3]),
  Gynogenic = rowMeans(logcpm_germline[, 4:6]),
  Androgenic = rowMeans(logcpm_germline[, 7:9])
)
logcpm_germline_means

logcpm_germline_means$Gene <- rownames(logcpm_germline_means)
# Reshape the data to long format for easier plotting
logcpm_germline_means_long <- pivot_longer(
  logcpm_germline_means,
  cols = c(Gynogenic, Androgenic),
  names_to = "Group",
  values_to = "Value"
)

logcpm_germline_means_long_autosomes <- logcpm_germline_means_long[logcpm_germline_means_long$Gene %in% autosomal_genes$V2, ]
logcpm_germline_means_long_X <- logcpm_germline_means_long[logcpm_germline_means_long$Gene %in% X_genes$V2, ]

# Plot using ggplot
male_vs_female_germline_all <- ggplot(logcpm_germline_means_long, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("A") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(linewidth = 2, shape = NA)))

male_vs_female_germline_autosomes <- ggplot(logcpm_germline_means_long_autosomes, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("B") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")

male_vs_female_germline_X <- ggplot(logcpm_germline_means_long_X, aes(x = Male, y = Value, color = Group)) +
  geom_point(shape = 1, alpha = 0.2) +                       # Scatter points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  labs(x = "Male expression (logCPM)", y = "Female expression (logCPM)", color = "Group") +
  scale_color_manual(values = c("chartreuse3", "purple")) +
  xlim(c(NA, 20)) +
  ylim(c(NA, 20)) +
  ggtitle("C") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 25, hjust = -0.1, vjust = 0.1, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "none")
# Fit a combined linear model with an interaction term
summary(lm(Value ~ Male * Group, data = logcpm_germline_means_long_autosomes))
summary(lm(Value ~ Male * Group, data = logcpm_germline_means_long_X))
# Interaction term is not significant, the slope isn't different
```
Put all the plots together using ggarrange for Fig. 5.
```
png(file="./outputs/SAS_comparison.png", height = 1000, width = 1000)
ggarrange(male_vs_female_germline_all, male_vs_female_germline_autosomes, male_vs_female_germline_X, male_vs_female_repro_all, male_vs_female_repro_autosomes, male_vs_female_repro_X, male_vs_female_nonrepro_all, male_vs_female_nonrepro_autosomes, male_vs_female_nonrepro_X, ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()
```
