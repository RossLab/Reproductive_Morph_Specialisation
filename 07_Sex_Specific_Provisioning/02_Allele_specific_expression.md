# Mapping embryo RNA-seq libraries to the X'X scaffolds, and summarising read counts #
I am interested in whether there is biased deposition of transcripts from the X or X' gametolog in the eggs of gynogenic females. To study this I first map the RNA reads again with STAR, but with stricter filters for read quality and mapping quality. 
```
LOC="inputs/"

# Trim RNA reads to retain only high quality reads to minimising sequencing errors
# Default fastp quality threshold is phred quality of >=15. I'm going to bump that up to 20.
for file in $(ls ${LOC}female*_1.fastq.gz)
do
        base=$(basename $file "_1.fastq.gz")
        fastp -q 20 --correction -i ${LOC}${base}_1.fastq.gz \
        -I ${LOC}${base}_2.fastq.gz \
        -o ${base}_1.ase.trimmed.fastq.gz -O ${base}_2.ase.trimmed.fastq.gz
done

# Sync files in or define file names for RNA mapping
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3.augustus.gtf'

mkdir Bcop_v3-chromosomes.STAR

 run genomeGenerate
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeSAindexNbases 13 \
--outFileNamePrefix Bcop_v3-chromosomes \
--sjdbGTFfile ${ANNO} \
--genomeDir Bcop_v3-chromosomes.STAR \
--genomeFastaFiles ${GENOME}

### Mapping RNA reads with STAR. Because I am mapping competitively between X and X' transcripts, I am only keeping primary alignments and not allowing multimapping.
echo "aligning RNAseq reads"
for file in ${SCRATCH}/*_1.ase.trimmed.fastq.gz
do
        echo "Processing file: $file"
        base=$(basename $file "_1.ase.trimmed.fastq.gz")
        STAR \
        --runThreadN 16 \
        --alignTranscriptsPerReadNmax 20000 \
        --outFilterMultimapNmax 10 \
        --alignEndsType EndToEnd \
        --outSAMattributes Standard \
        --outSAMprimaryFlag AllBestScore \
        --outSAMmultNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${base}_1.ase.trimmed.fastq.gz ${base}_2.ase.trimmed.fastq.gz \
        --readFilesCommand zcat \
        --outTmpDir ${base}.out \
        --outFileNamePrefix ${base}.STARase. \
        --genomeDir Bcop_v3-chromosomes.STAR
done

# AllBestScore outputs all alignments with the best score as primary alignments i.e. perfect multimappers, so I can filter them later
# outSAMmultiNmax 1 only prints the alignment with the highest score

echo "filtering BAM files"
for file in ${SCRATCH}/*.STARase.Aligned.sortedByCoord.out.bam
do
        base=$(basename $file ".STARase.Aligned.sortedByCoord.out.bam")
        samtools view -b -q 10 ${base}.STARase.Aligned.sortedByCoord.out.bam > ${base}.STARase.Aligned.sortedByCoord.filtered.out.bam
done
```
I then summarise the read counts with featureCounts. 
```
# Define file names
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3_geneid.augustus.gtf'

echo "running featureCounts"
for file in $(ls *.STARase.Aligned.sortedByCoord.filtered.out.bam)
do
        echo "running featureCounts for" $file
        base=$(basename $file ".STARase.Aligned.sortedByCoord.filtered.out.bam")
        featureCounts -T 5 -p --countReadPairs -M --primary -a ${ANNO} -t exon -g gene_id -o ${base}.ase.featureCounts.txt ${base}.STARase.Aligned.sortedByCoord.filtered.out.bam
done
```
\
I can now take the read count files into R for allele-specific analysis. 
\
Packages I need:
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
\
I first read in my design table, read counts, and normalise genes counts by gene length (which may differ between the gametologs).
```
## Files I need: 
# RNA seq design 
design_ase <- read.csv("r_input/maternal_deposit_design.csv")
design_ase <- design_ase[1:3,]
design_ase$path <- paste0("r_input/", design_ase$sample_id, ".ase.featureCounts.txt",sep="")
design_ase


# All genes 
genes_and_chromosomes <- read.table("r_input/Bcop_v3_genes_and_chromosomes.tsv")
X_Inv <- genes_and_chromosomes[(genes_and_chromosomes$V1 == "X") | (genes_and_chromosomes$V1 == "Inversion"), ]
nrow(X_Inv) # 7371 (3898 X genes and 3473 Inv genes)

# Which X and X' genes are orthologs, output of orthofinder. 
X_Inv_orthogroups <- read.table("r_input/X_vs_Inv_longest_OG_HOG.tsv", header = T, sep="\t")
X_Inv_SCO <- read.table("r_input/X_vs_Inv_longest_Orthogroups_SCOs.tsv") 
nrow(X_Inv_SCO) #2703

X_Inv_SCO_genes <- X_Inv_orthogroups[X_Inv_orthogroups$HOG %in% X_Inv_SCO$V1, ]
nrow(X_Inv_SCO_genes) # 2703
X_Inv_SCO_genes <- X_Inv_SCO_genes[, c("Bcop_Inv_genes", "Bcop_X_genes")]
colnames(X_Inv_SCO_genes) <- c("Inv_transcript", "X_transcript")
X_Inv_SCO_genes[] <- lapply(X_Inv_SCO_genes, function(x) sub("\\.t.*", "", x)) # Removes the transcript number
X_Inv_SCO_genes[which(X_Inv_SCO_genes$X_transcript == "g13114"), 2] <- "g13114.t1"
head(X_Inv_SCO_genes) 

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
ase_raw_gene_counts <- ase_raw_gene_counts[ , 2:4]
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

write.csv(X_Inv_SCO_genes_count, "outputs/allele_specific_gene_count_maternal_deposit.csv", row.names = FALSE)
```
\
Now that I have all the length normalised gene counts for X genes and their ortholog on the Inv, I can carry out stats to see if they are differentially expressed.
```
gene_count_04hr <- X_Inv_SCO_genes_count
nrow(gene_count_04hr) # 2703, same number as single orthologs between X and Inv

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
write.csv(allele_specific_summary, "outputs/allele_specific_exp_summary_04hrs.csv", row.names = FALSE)
```
\
I add in information about the rate of molecular evolution for each gametolog, obtained in an earlier analysis. 
```
allele_specific_summary_readin <- read.csv("outputs/allele_specific_exp_summary_04hrs.csv")

nrow(allele_specific_summary_readin) # 2703
nrow(allele_specific_summary_readin[allele_specific_summary_readin$pval < 0.05, ]) # 1433

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
  geom_point(size = 2, alpha =0.6) +
  scale_color_manual(name = "Expression", labels = c("Inversion-biased", "Unbiased", "X-biased"), values = c("maroon4", "gray50", "plum")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Average expression (log2 TPM)",
    y = "Allele-specific expression (log2 Inv - X)") +
  labs(title = "C - Allele-specific expression") + theme_bw() +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title.position = "plot",
    plot.title = element_text(size = 22, hjust = 0.05, vjust = 2, margin = margin(b = 10), face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.title=element_text(size=16), 
    legend.text=element_text(size=16),
    legend.spacing.y = unit(0.2, "cm"),
    legend.box.spacing = unit(0.2, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2))
MAplot_ASE_04hr
```
\
As with female tissue, I am  interested in whether differential deposition between gynogenic and androgenic females in X-linked genes is associated with biased deposition between X and X' gametologs in gynogenic females.
```
XpX_mapping_androvsgyno_04hrs_readin <- read.csv("outputs/XpX_04hrs_androvsgyno_all.csv")
colnames(XpX_mapping_androvsgyno_04hrs_readin)[1] <- "X_transcript"
head(XpX_mapping_androvsgyno_04hrs_readin)
nrow(XpX_mapping_androvsgyno_04hrs_readin) #1761

allele_exp_04hrs_readin <- read.csv("outputs/allele_specific_exp_summary_04hrs.csv")
nrow(allele_exp_04hrs_readin) #2703

predict <- inner_join(XpX_mapping_androvsgyno_04hrs_readin, allele_exp_04hrs_readin, by = "X_transcript")
predict <- predict[, c(1, 2, 6, 9, 10)]
colnames(predict) <- c("X_transcript", "GA_logFC", "GA_FDR", "XpX_logFC", "XpX_pval")

predict_plot <- ggplot(data = predict, aes(x = GA_logFC, y = XpX_logFC)) + 
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "purple", alpha = 0.1)  + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= "chartreuse3", alpha = 0.1) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= "white", alpha = 0.1) + 
  annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= "white", alpha = 0.1)+
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "spearman", size = 4.5, label.x.npc = 0.5, label.y.npc = 0.1) +
  ggtitle("D - X-linked genes DGE vs ASE") +
  xlab("Androgenic vs Gynogenic logFC") + ylab("X vs Inversion Allele logFC") +
  theme_minimal() +
  theme(plot.title.position = "plot",
        axis.line = element_blank(),
        # panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 22, hjust = 0.05, vjust = 2, face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.position = "none")
predict_plot
cor.test(predict$GA_logFC, predict$XpX_logFC, method = "spearman", exact = FALSE)
```