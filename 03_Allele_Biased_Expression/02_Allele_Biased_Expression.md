# Allele-Biased Expression Analysis in R #
The big difference between conventional allele-specific expression and my analysis, is that I cannot assume the gametologs on the X and X' are the same length (whether due to degeneration or incomplete assembly on the X'). To counter this, I incorporate a measure of the length difference between the gametologs into my model as an "offset". My model therefore asks, per pair of gametologs, whether they have different read counts beyond what can be expected from differences in length.  

This is done in R (v4.3.2).

Required packages:
```
Packages:
library(tidyverse) # v2.0.0
library(ggpubr) # v0.6.0
library(dplyr) # v1.1.4
library(GenomicRanges) # v1.54.1
library(glmmTMB) # v1.1.11
library(rstatix) # v0.7.2
```
## Files I need: 
# RNA seq design 
design <- read.csv("r_input/allele_specific_expression_gyno_design.csv", header = T)

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

# Read in featureCount files
design$path <- paste0("r_input/", design$sample_id, "_featureCounts.txt",sep="")
design

## Get gene lengths for each gene ##
# Build a TxDb object from  GTF
txdb <- GenomicFeatures::makeTxDbFromGFF("r_input/Bcop_v3_geneid.augustus.gtf")
# Get exons by gene
exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
# Calculate total non-overlapping exon length per gene
gene_lengths <- sum(width(reduce(exons_by_gene)))
length(gene_lengths) # 24050

X_Inv_lengths <- gene_lengths[names(gene_lengths) %in% X_Inv_SCO_genes$Inv_transcript | names(gene_lengths) %in% X_Inv_SCO_genes$X_transcript] #Filters for genes only in the X or inversion 
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

## Calculates p0, a measure of length difference between the gametologs
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
featureCounts <- lapply(design$path, function(x)read.table(x, header=T)) # applies the read.table() function to each file path, reading each file into a data frame
names(featureCounts) <- design$sample_id # Names each dataframe in the list
featureCounts <- lapply(featureCounts, function(x) x[, c(1, 7)]) # For each data frame in the list, only columns 1 and 7 are kept.
featureCounts <- lapply(names(featureCounts), function(name) {
  x <- featureCounts[[name]]
  colnames(x)[2] <- name
  x 
}) # Assigns the corresponding sample ID (name) as the name of the second column (the count data column)
names(featureCounts) <- design$sample_id # Names each dataframe in the list because the name got lost
for (name in names(featureCounts)) {
  # Extract the current table
  x <- featureCounts[[name]]
  # Extract column 2, setting it as a new column in counts table with the sample ID as the column name
  X_Inv_SCO_genes_length_long <- merge(X_Inv_SCO_genes_length_long, x, by = "Geneid", all.x = TRUE)
}
X_Inv_SCO_genes_length_long <- X_Inv_SCO_genes_length_long %>%
  arrange(Number)
nrow(X_Inv_SCO_genes_length_long) #5406
```
\
Now that I have a table of the X-Inv gametologs, their read count, and their length difference, I use a generalised linear model with betabinomial distribution to assess differential expression between pairs of homologs. As before, this is done by tissue. 
```
############################ SOMATIC NON-REPRODUCTIVE ALLELE-BIASED EXP ############################ 
design
gene_counts_nonrepro <- X_Inv_SCO_genes_length_long[, c(1:5, 12:14)]
sample_cols <- grep("^B\\d+", names(gene_counts_nonrepro), value = TRUE)
gene_counts_nonrepro_long <- gene_counts_nonrepro %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "Sample", values_to = "Count") %>%
  group_by(Number, Sample) %>%
  summarise(
    Count_Inv = sum(Count[Allele == "Inv"], na.rm = TRUE),
    Count_X   = sum(Count[Allele == "X"],   na.rm = TRUE),
    logit_p0  = dplyr::first(logit_p0),     # same within a Pair
    .groups   = "drop"
  ) 
nrow(gene_counts_nonrepro_long) # 8109
head(gene_counts_nonrepro_long, n =30) 

## Get rid of zero counts -> unexpressed 

gene_counts_nonrepro_long_filt <- gene_counts_nonrepro_long %>%
  mutate(Total = Count_Inv + Count_X) %>% 
  filter(Total > 0)
nrow(gene_counts_nonrepro_long_filt) #7506
write.csv(gene_counts_nonrepro_long_filt, "output/nonrepro_gene_counts_new.csv")

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
all_numbers <- unique(gene_counts_nonrepro_long_filt$Number)
chunk_size  <- 400
chunks      <- split(all_numbers, ceiling(seq_along(all_numbers) / chunk_size))

results_list <- vector("list", length(chunks))

for (i in seq_along(chunks)) {
  cat("Fitting chunk", i, "of", length(chunks), "\n")
  
  results_list[[i]] <- gene_counts_nonrepro_long_filt %>%
    filter(Number %in% chunks[[i]]) %>%
    group_by(Number) %>%
    nest() %>%
    mutate(fit = map(data, \(d) suppressMessages(suppressWarnings(fit_bb_tmb(d))))) %>%
    unnest(fit) %>%
    select(-data)
  
  gc()
}

nonrepro_results_expressed <- bind_rows(results_list)

# Re-attach all-zero pairs from unfiltered table
nonrepro_zero_pairs <- gene_counts_nonrepro_long %>%
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

nonrepro_results_all <- bind_rows(nonrepro_results_expressed, nonrepro_zero_pairs) %>%
  arrange(Number)

converged_pdr <- which(nonrepro_results_all$converged & !is.na(nonrepro_results_all$p_value))
nonrepro_results_all$q_value <- NA_real_
nonrepro_results_all$q_value[converged_pdr] <- p.adjust(
  nonrepro_results_all$p_value[converged_pdr], method = "BH"
)

nonrepro_results_all <- nonrepro_results_all %>%
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

print(table(nonrepro_results_all$bias, useNA = "ifany"))

nonrepro_results_all<-left_join(nonrepro_results_all, X_Inv_SCO_genes, by = "Number")
nonrepro_results_all <- nonrepro_results_all[, c(1,15:16,2:14) ]
nrow(nonrepro_results_all) #2703
write.csv(nonrepro_results_all, "output/nonrepro_allele_biased_exp_new.csv")

##### MA Plot #####
# Compute M and A from the raw counts
nonrepro_MA_cal <- gene_counts_nonrepro_long_filt %>%
  group_by(Number) %>%
  summarise(
    # A: mean of log2 total counts (average expression level across samples)
    mean_log_total = mean(log2(Total + 1)),
    
    # M: log2 fold change Inv/X, averaged across samples
    # add pseudocount of 0.5 to avoid log(0)
    mean_log2FC    = mean(log2((Count_Inv + 0.5) / (Count_X + 0.5))),
    .groups = "drop"
  )
nonrepro_results_all <- nonrepro_results_all %>%
  left_join(nonrepro_MA_cal, by = "Number")

nonrepro_results_all <- nonrepro_results_all %>%
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

MAplot_ASE_nonrepro <- ggplot(data = nonrepro_results_all, aes(x = mean_log_total, y = mean_log2FC, colour = bias_plot)) +
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
  ggtitle("C - Somatic non-reproductive") + theme_bw() +
  theme(plot.title.position = "plot",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.position = "none") 
MAplot_ASE_nonrepro

#### dNdS plot ####
dNdS_readin <- read.csv("C:/Users/s2556496/Desktop/All/Morph_Specialisation_Gene_Expression/B_coprophila_morph_gene_divergence_FINAL/04_dnds/output/X_Inv_ortholog_dNdS.csv")

exp_dNdS_summary <- merge(nonrepro_results_all, dNdS_readin, by = "X_transcript")
exp_dNdS_summary <- exp_dNdS_summary |>
  select(X_transcript, Inv_transcript.x, bias, Inv_dNdS, X_dNdS, dNdS_diff)
nrow(exp_dNdS_summary)
exp_dNdS_summary <- exp_dNdS_summary[exp_dNdS_summary$bias != "model_failed",]
nrow(exp_dNdS_summary)
```
\
Allele-specific expression shows that a number of inversion gametologs are upregulated, which is surprising and out of the null expectation of degradation. I am interested in using rate of molecular evolution (rate of nonsynonymous/synonymous substitution) as further evidence that this Inv-biased expression is due to purifying selection, and therefore likely to be functional.
```
set.seed(123)  # Set seed once
### Do the different bias groups have different dNdS differences between the X and Inv homolog?
exp_dNdS_summary$bias <- factor(exp_dNdS_summary$bias)
exp_dNdS_summary$bias <- relevel(exp_dNdS_summary$bias, ref = "unbiased")
welch_anova_test(exp_dNdS_summary, dNdS_diff ~ bias)
stat_gh_nonrepro <- games_howell_test(exp_dNdS_summary, dNdS_diff ~ bias) %>%
  add_xy_position(x = "bias", fun = "max", step.increase = 0.10)
stat_gh_nonrepro$label <- stat_gh_nonrepro$p.adj.signif 

exp_dNdS_plot_nonrepro <- ggplot(data = exp_dNdS_summary, aes(x = bias, y = dNdS_diff, col = bias)) +
  scale_color_manual("Expression bias",labels = c("Unbiased", "Inversion allele", "X allele"), values=c("grey50","maroon4", "plum")) +
  scale_fill_manual(values = c("grey", "plum"))+
  scale_x_discrete(labels = c(unbiased = "Unbiased", Inv_biased = "Inversion", X_biased = "X"))+
  geom_boxplot(alpha = 0) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 0.02, height = 0)) +
  ggtitle("F - Somatic non-reproductive") +
  ylab("dNdS of X allele - dNdS of Inversion allele") + xlab("Expression bias")+
  stat_pvalue_manual(stat_gh_nonrepro, label = "label", tip.length = 0, step.increase = 0.06, y.position = 0.5, bracket.size = 0.75) +  
  coord_cartesian(ylim = c(-1, 0.75)) +
  theme_minimal()+
  theme(plot.title.position = "plot",
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 1.2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.position = "none")
exp_dNdS_plot_nonrepro
```
\
It seems that genes with Inv-biased expression also have more conserved evolution in the inversion gametolog, relative to the X gametolog. I want to further test if there is a link between differential expression between female morphs and allele-specific expression within X'X gynogenic females.
```
# Read in results of X-linked differential gene expression between gynogenic and androgenic females 
XpX_mapping_androvsgyno_nonrepro_readin <- read.csv("/09_X_edgeR/outputs/XpX_nonrepro_androvsgyno_all.csv")
colnames(XpX_mapping_androvsgyno_nonrepro_readin)[1] <- "X_transcript"
head(XpX_mapping_androvsgyno_nonrepro_readin)
nrow(XpX_mapping_androvsgyno_nonrepro_readin) #2336

allele_exp_nonrepro_readin <- read.csv("output/nonrepro_allele_biased_exp_new.csv")
nrow(allele_exp_nonrepro_readin) #2703
allele_exp_nonrepro_readin$Number <- as.factor(allele_exp_nonrepro_readin$Number)
nonrepro_MA_cal <- gene_counts_nonrepro_long_filt %>%
  group_by(Number) %>%
  summarise(
    # A: mean of log2 total counts (average expression level across samples)
    mean_log_total = mean(log2(Total + 1)),
    
    # M: log2 fold change Inv/X, averaged across samples
    # add pseudocount of 0.5 to avoid log(0)
    mean_log2FC    = mean(log2((Count_Inv + 0.5) / (Count_X + 0.5))),
    .groups = "drop"
  )
allele_exp_nonrepro_readin <- allele_exp_nonrepro_readin %>%
  left_join(nonrepro_MA_cal, by = "Number")
nrow(allele_exp_nonrepro_readin) #2703

nonrepro_predict <- inner_join(XpX_mapping_androvsgyno_nonrepro_readin, allele_exp_nonrepro_readin, by = "X_transcript")
nonrepro_predict <- nonrepro_predict[, c(1, 2, 6, 24, 21)]
colnames(nonrepro_predict) <- c("X_transcript", "GA_logFC", "GA_FDR", "XpX_logFC", "XpX_pval")

nonrepro_predict_plot <- ggplot(data = nonrepro_predict, aes(x = GA_logFC, y = XpX_logFC)) + 
  annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "purple", alpha = 0.1)  + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= "chartreuse3", alpha = 0.1) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf, fill= "white", alpha = 0.1) + 
  annotate("rect", xmin = 0, xmax = -Inf, ymin = Inf, ymax = 0, fill= "white", alpha = 0.1) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", col = "black") +
  stat_cor(method = "spearman", size = 4.5, label.x.npc = 0.5, label.y.npc = 0.1) +
  ggtitle("I - Somatic non-reproductive") +
  xlab("Androgenic vs Gynogenic logFC") + ylab("X vs Inversion Allele logFC") +
  theme_minimal() +
  theme(plot.title.position = "plot",
        axis.line = element_blank(),
        # panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 20, hjust = 0.1, vjust = 1, face = "bold"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.position = "none")
nonrepro_predict_plot
cor.test(nonrepro_predict$GA_logFC, nonrepro_predict$XpX_logFC, method = "spearman", exact = FALSE)
```
I do the same for germline and somatic reproductive tissue, and use use ggarrange to put all the plots together for Fig 4.
```
png(file="./output/all_exp_dNdS.png", height = 1100, width = 1000)
ggarrange(ggarrange(MAplot_ASE_germline, MAplot_ASE_repro, MAplot_ASE_nonrepro, exp_dNdS_plot_germline, exp_dNdS_plot_repro, exp_dNdS_plot_nonrepro, nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom"), ggarrange(germline_predict_plot, repro_predict_plot, nonrepro_predict_plot, nrow = 1, ncol = 3), ncol = 1, nrow = 2, heights = c(2,1))
dev.off()

```
