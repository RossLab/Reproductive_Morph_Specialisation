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
Packages:
library(tidyverse) # v2.0.0
library(ggpubr) # v0.6.0
library(dplyr) # v1.1.4
library(GenomicRanges) # v1.54.1
library(glmmTMB) # v1.1.11
library(rstatix) # v0.7.2
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

X_Inv_lengths <- gene_lengths[names(gene_lengths) %in% X_Inv_SCO_genes$Inv_transcript | names(gene_lengths) %in% X_Inv_SCO_genes$X_transcript]
length(X_Inv_lengths) # 5406

X_Inv_lengths_df <- data.frame(Inv_transcript = names(X_Inv_lengths), Inv_length = X_Inv_lengths)
X_Inv_SCO_genes <- left_join(X_Inv_SCO_genes, X_Inv_lengths_df, by = "Inv_transcript") 
nrow(X_Inv_SCO_genes)
X_Inv_lengths_df <- data.frame(X_transcript = names(X_Inv_lengths), X_length = X_Inv_lengths)
X_Inv_SCO_genes <- left_join(X_Inv_SCO_genes, X_Inv_lengths_df, by = "X_transcript") 
nrow(X_Inv_SCO_genes) #2703
ggplot(data = X_Inv_SCO_genes, aes(x = X_length, y = Inv_length)) +
  geom_point(size = 2, alpha = 0.2) +
  geom_smooth() + geom_abline(slope = 1, intercept = 0)
## There is a length difference between the X and Inv transcripts, especially for longer genes 
X_Inv_SCO_genes <- X_Inv_SCO_genes %>%
  mutate(Number = as.factor(row_number()))%>%
  mutate(p0       = Inv_length / (Inv_length + X_length),
         logit_p0 = qlogis(p0))
X_Inv_SCO_genes_length_long <- bind_rows(transmute(X_Inv_SCO_genes, Geneid = Inv_transcript, Number, Allele = "Inv", Len = Inv_length, logit_p0),
                                         transmute(X_Inv_SCO_genes, Geneid = X_transcript, Number, Allele = "X", Len = X_length, logit_p0)) 

X_Inv_SCO_genes_length_long <- X_Inv_SCO_genes_length_long %>%
  arrange(Number)
nrow(X_Inv_SCO_genes_length_long) # 5406
head(X_Inv_SCO_genes_length_long)

### Read in the featureCounts counts, each column is one sample
featureCounts <- lapply(design_ase$path, function(x)read.table(x, header=T)) # applies the read.table() function to each file path, reading each file into a data frame
names(featureCounts) <- design_ase$sample_id # Names each dataframe in the list
featureCounts <- lapply(featureCounts, function(x) x[, c(1, 7)]) # For each data frame in the list, only columns 1 and 7 are kept.
featureCounts <- lapply(names(featureCounts), function(name) {
  x <- featureCounts[[name]]
  colnames(x)[2] <- name
  x 
}) # Assigns the corresponding sample ID (name) as the name of the second column (the count data column)
names(featureCounts) <- design_ase$sample_id # Names each dataframe in the list because the name got lost
for (name in names(featureCounts)) {
  # Extract the current table
  x <- featureCounts[[name]]
  # Extract column 2, setting it as a new column in counts table with the sample ID as the column name
  X_Inv_SCO_genes_length_long <- merge(X_Inv_SCO_genes_length_long, x, by = "Geneid", all.x = TRUE)
}
X_Inv_SCO_genes_length_long <- X_Inv_SCO_genes_length_long %>%
  arrange(Number)
nrow(X_Inv_SCO_genes_length_long) #5406

## Allele-specific expression ##
design_ase
gene_counts_matdep <- X_Inv_SCO_genes_length_long
sample_cols <- grep("^female_0-4h_rep\\d+", names(gene_counts_matdep), value = TRUE)
gene_counts_matdep_long <- gene_counts_matdep %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "Sample", values_to = "Count") %>%
  group_by(Number, Sample) %>%
  summarise(
    Count_Inv = sum(Count[Allele == "Inv"], na.rm = TRUE),
    Count_X   = sum(Count[Allele == "X"],   na.rm = TRUE),
    logit_p0  = dplyr::first(logit_p0),     # same within a Pair
    .groups   = "drop"
  ) 
nrow(gene_counts_matdep_long) # 8109
head(gene_counts_matdep_long, n =30) 

## Get rid of zero counts -> unexpressed 

gene_counts_matdep_long_filt <- gene_counts_matdep_long %>%
  mutate(Total = Count_Inv + Count_X) %>% 
  filter(Total > 0)
nrow(gene_counts_matdep_long_filt) #6071

#### Running the model ####
# Function for fitting the model
fit_bb_tmb <- function(df) {
  
  n_samples <- nrow(df)
  
  # Cannot fit betabinomial with a single observation
  if (n_samples < 2) {
    return(tibble(
      intercept = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, fitted_p_inv = NA_real_,
      null_p_inv     = plogis(mean(df$logit_p0)),
      grad_max       = NA_real_,
      total_reads    = sum(df$Total),
      n_samples_used = n_samples,
      converged      = FALSE,
      flag           = "single_sample"
    ))
  }
  
  tryCatch({
    fit <- glmmTMB(
      cbind(Count_Inv, Count_X) ~ 1 + offset(logit_p0),
      family = betabinomial(link = "logit"),
      data   = df
    )
    
    # Convergence checks
    pdHess     <- isTRUE(fit$sdr$pdHess)
    grad       <- tryCatch(max(abs(fit$sdr$gradient.fixed)), error = function(e) NA_real_)
    grad_ok    <- isTRUE(grad < 0.001)
    
    s          <- summary(fit)$coefficients$cond
    intercept  <- s["(Intercept)", "Estimate"]
    se         <- s["(Intercept)", "Std. Error"]
    z          <- s["(Intercept)", "z value"]
    p_val      <- s["(Intercept)", "Pr(>|z|)"]
    separation <- isTRUE(abs(intercept) > 10 | se > 10)
    converged  <- pdHess & grad_ok & !separation
    
    flag <- case_when(
      separation ~ "separation",
      !pdHess    ~ "non_pd_hessian",
      !grad_ok   ~ paste0("large_gradient:", round(grad, 4)),
      TRUE       ~ "ok"
    )
    
    mean_offset  <- mean(df$logit_p0)
    fitted_p_inv <- plogis(mean_offset + intercept)
    null_p_inv   <- plogis(mean_offset)
    
    tibble(
      intercept, se, z, p_value = p_val,
      fitted_p_inv, null_p_inv,
      grad_max       = grad,
      total_reads    = sum(df$Total),
      n_samples_used = n_samples,
      converged      = converged,
      flag           = flag
    )
    
  }, error = function(e) {
    tibble(
      intercept = NA_real_, se = NA_real_, z = NA_real_,
      p_value = NA_real_, fitted_p_inv = NA_real_,
      null_p_inv     = plogis(mean(df$logit_p0)),
      grad_max       = NA_real_,
      total_reads    = sum(df$Total),
      n_samples_used = n_samples,
      converged      = FALSE,
      flag           = paste0("error: ", conditionMessage(e))
    )
  })
}

### Running in chunks for memory 
all_numbers <- unique(gene_counts_matdep_long_filt$Number)
chunk_size  <- 400
chunks      <- split(all_numbers, ceiling(seq_along(all_numbers) / chunk_size))

results_list <- vector("list", length(chunks))

for (i in seq_along(chunks)) {
  cat("Fitting chunk", i, "of", length(chunks), "\n")
  
  results_list[[i]] <- gene_counts_matdep_long_filt %>%
    filter(Number %in% chunks[[i]]) %>%
    group_by(Number) %>%
    nest() %>%
    mutate(fit = map(data, \(d) suppressMessages(suppressWarnings(fit_bb_tmb(d))))) %>%
    unnest(fit) %>%
    select(-data)
  
  gc()
}

matdep_results_expressed <- bind_rows(results_list)

# Re-attach all-zero pairs from unfiltered table
matdep_zero_pairs <- gene_counts_matdep_long %>%
  group_by(Number) %>%
  filter(all(Count_Inv == 0 & Count_X == 0)) %>%
  summarise(logit_p0 = dplyr::first(logit_p0), .groups = "drop") %>%
  transmute(
    Number,
    intercept      = NA_real_, se = NA_real_, z = NA_real_,
    p_value        = NA_real_, fitted_p_inv = NA_real_,
    null_p_inv     = plogis(logit_p0),
    grad_max       = NA_real_,
    total_reads    = 0L,
    n_samples_used = 0L,
    converged      = FALSE,
    flag           = "all_zero"
  )

matdep_results_all <- bind_rows(matdep_results_expressed, matdep_zero_pairs) %>%
  arrange(Number)

converged_pdr <- which(matdep_results_all$converged & !is.na(matdep_results_all$p_value))
matdep_results_all$q_value <- NA_real_
matdep_results_all$q_value[converged_pdr] <- p.adjust(
  matdep_results_all$p_value[converged_pdr], method = "BH"
)

matdep_results_all <- matdep_results_all %>%
  mutate(
    bias = case_when(
      flag == "all_zero"                                         ~ "unexpressed",
      flag == "separation" & total_reads >= 10 & intercept > 0  ~ "Inv_biased",
      flag == "separation" & total_reads >= 10 & intercept < 0  ~ "X_biased",
      !converged                                                 ~ "model_failed",
      q_value < 0.05 & intercept > 0                            ~ "Inv_biased",
      q_value < 0.05 & intercept < 0                            ~ "X_biased",
      q_value >= 0.05                                            ~ "unbiased"
    )
  )

print(table(matdep_results_all$bias, useNA = "ifany"))

matdep_results_all<-left_join(matdep_results_all, X_Inv_SCO_genes, by = "Number")
matdep_results_all <- matdep_results_all[, c(1,15:16,2:14) ]
nrow(matdep_results_all) #2703
write.csv(matdep_results_all, "/outputs/matdep_allele_biased_exp_new.csv")

##### MA Plot #####
# Compute M and A from the raw counts
matdep_MA_cal <- gene_counts_matdep_long_filt %>%
  group_by(Number) %>%
  summarise(
    # A: mean of log2 total counts (average expression level across samples)
    mean_log_total = mean(log2(Total + 1)),
    
    # M: log2 fold change Inv/X, averaged across samples
    # add pseudocount of 0.5 to avoid log(0)
    mean_log2FC    = mean(log2((Count_Inv + 0.5) / (Count_X + 0.5))),
    .groups = "drop"
  )
matdep_results_all <- matdep_results_all %>%
  left_join(matdep_MA_cal, by = "Number")

matdep_results_all <- matdep_results_all %>%
  filter(bias != "unexpressed") %>%           # exclude all-zero pairs
  mutate(
    bias_plot = case_when(
      bias == "Inv_biased"   ~ "Inv biased",
      bias == "X_biased"     ~ "X biased",
      bias == "unbiased"     ~ "Unbiased",
      TRUE                   ~ "Model failed"
    ),
    bias_plot = factor(bias_plot, levels = c("Inv biased", "X biased", 
                                             "Unbiased", "Model failed"))
  )

MAplot_ASE_matdep <- ggplot(data = matdep_results_all, aes(x = mean_log_total, y = mean_log2FC, colour = bias_plot)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  scale_colour_manual("Expression bias", values = c(
    "Inv biased"   = "maroon4",
    "X biased"     = "plum",
    "Unbiased"     = "grey50",
    "Model failed" = "grey85"
  )) +
  labs(
    x = "Average expression (log2 counts)",
    y = "Allele-specific expression (log2 Inv - X)"
  ) +
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
MAplot_ASE_matdep
```
\
As with female tissue, I am  interested in whether differential deposition between gynogenic and androgenic females in X-linked genes is associated with biased deposition between X and X' gametologs in gynogenic females.
```
XpX_mapping_androvsgyno_04hrs_readin <- read.csv("outputs/XpX_04hrs_androvsgyno_all.csv")
colnames(XpX_mapping_androvsgyno_04hrs_readin)[1] <- "X_transcript"
head(XpX_mapping_androvsgyno_04hrs_readin)
nrow(XpX_mapping_androvsgyno_04hrs_readin) #1761

allele_exp_04hrs_readin <- read.csv("outputs/matdep_allele_biased_exp_new.csv")
allele_exp_04hrs_readin$Number <- as.factor(allele_exp_04hrs_readin$Number)
nrow(allele_exp_04hrs_readin) #2703
allele_exp_04hrs_readin_MA_cal <- gene_counts_matdep_long_filt %>%
  group_by(Number) %>%
  summarise(
    # A: mean of log2 total counts (average expression level across samples)
    mean_log_total = mean(log2(Total + 1)),
    
    # M: log2 fold change Inv/X, averaged across samples
    # add pseudocount of 0.5 to avoid log(0)
    mean_log2FC    = mean(log2((Count_Inv + 0.5) / (Count_X + 0.5))),
    .groups = "drop"
  )
allele_exp_04hrs_readin <- allele_exp_04hrs_readin %>%
  left_join(allele_exp_04hrs_readin_MA_cal, by = "Number")
nrow(allele_exp_04hrs_readin) #2703

predict <- inner_join(XpX_mapping_androvsgyno_04hrs_readin, allele_exp_04hrs_readin, by = "X_transcript")
predict <- predict[, c(1, 2, 6, 24, 21)]
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