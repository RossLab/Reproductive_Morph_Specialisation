# Differential gene expression analysis of X-linked genes #

To assess whether gynogenic and androgenic females deposit different RNA into their eggs as a method of sex-specific provisioning, I carried out the same analyses as I did with the separate female tissues, but with maternally-depositted RNA in male and female eggs instead.

I first look to see if there are differentially-depositted transcripts that come from the X or X' inversion. I already have the X'X gametologs from a previous Orthofinder run, so that does not need to be regenerated. 

## Mapping of RNA reads ##
RNA reads extracted from male and female pooled embryos were downloaded. FastQC and fastp were used to check the quality of and trim the reads. 

```
# fastqc pre-trim
fastqc -t 4 *.fastq.gz

# trim reads
for file in $(ls *_1.fastq.gz)
do
        base=$(basename $file "_1.fastq.gz")
        fastp -i ${base}_1.fastq.gz \
        -I ${base}_2.fastq.gz \
        -o ${base}_1.trimmed.fastq.gz -O ${base}_2.trimmed.fastq.gz
done

# fastqc post-trim
fastqc -t 4 *.trimmed.fastq.gz
```
\
STAR was used to map the trimmed reads to the _Bradysia coprophila_ v3 reference genome. 
```
# Sync files in or define file names
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3.augustus.gtf'
LOC='/07_maternal_deposit/inputs/'

mkdir Bcop_v3-chromosomes.STAR

# run genomeGenerate
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeSAindexNbases 12 \
--outFileNamePrefix Bcop_v3-chromosomes \
--sjdbGTFfile ${ANNO} \
--genomeDir Bcop_v3-chromosomes.STAR \
--genomeFastaFiles ${GENOME}

echo "aligning RNAseq reads"
for file in $(ls ${LOC}*_1.trimmed.fastq.gz)
do
        base=$(basename $file "_1.trimmed.fastq.gz")
        STAR \
        --runThreadN 16 \
        --alignTranscriptsPerReadNmax 20000 \
        --alignEndsType EndToEnd \
        --outSAMattributes Standard \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${LOC}${base}_1.trimmed.fastq.gz ${LOC}${base}_2.trimmed.fastq.gz \
        --readFilesCommand zcat \
        --outTmpDir ${base}.out \
        --outFileNamePrefix ${base}.STAR. \
        --genomeDir Bcop_v3-chromosomes.STAR
done
```
\
FeatureCounts from the Subreads package was used to for read summation up to gene level. Multimappers were counted fractionally.
```
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3_geneid.augustus.gtf'
LOC='07_maternal_deposit/outputs/'

echo "running featureCounts"
for file in $(ls ${LOC}*.STAR.Aligned.sortedByCoord.out.bam)
do
        echo "running featureCounts for" $file
        base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
        featureCounts -T 5 -p --countReadPairs -M --fraction -a ${ANNO} -t exon -g gene_id -o ${base}.multi.featureCounts.txt ${LOC}${base}.STAR.Aligned.sortedByCoord.out.bam
done
```
## Differential gene expression with edgeR ##
I read in the Orthofinder output, as well as a design table and a file mapping genes to their chromosomes. I make a table with the X and X' count per gametolog group, so that I can sum the read counts up.
```
# Which X and X' genes are orthologs, output of orthofinder. 
X_Inv_orthogroups <- read.table("r_input/X_vs_Inv_longest_OG_HOG.tsv", header = T, sep="\t")
X_Inv_SCO <- read.table("r_input/X_vs_Inv_longest_Orthogroups_SCOs.tsv")
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

design_XpX <- read.csv("r_input/maternal_deposit_design.csv")
design_XpX$path_XpX <- paste0("r_input/", design_XpX$sample_id, ".multi.featureCounts.txt",sep="")
design_XpX
design_XpX_04hrs <- design_XpX

## Genes and the chromosome they belong to ##
genes_and_chromosomes <- read.table("r_input/Bcop_v3_genes_and_chromosomes.tsv")

### Here I need to wrangle the data into a format where I sum the expression count of the X and Inv orthologs, per sample. 
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
expr_master <- expr_master[, c(2:7)]
nrow(expr_master) #2703
head(expr_master)
```
\
With this table, I can carry out DGE analysis. 
```
## Read the counts into a DGEList, specifying the design ##
DGEList_XpX <- DGEList(expr_master, group = design_XpX$sex)
nrow(DGEList_XpX$counts) #2703

## Filter lowly expressed genes ##
keep <- filterByExpr(DGEList_XpX)
table(keep) # Keep 1761 genes, get rid of 942 genes 
DGEList_XpX <- DGEList_XpX[keep, , keep.lib.sizes=FALSE]
nrow(DGEList_XpX$counts) #1761

## Normalisation ##
DGEList_XpX <- normLibSizes(DGEList_XpX)
head(DGEList_XpX)
print(DGEList_XpX$samples)

# Constructing a model matrix for glm 
design_XpX_04hrs_matrix <- model.matrix(~0 + design_XpX$sex)
rownames(design_XpX_04hrs_matrix) <- design_XpX_04hrs$sample_id
colnames(design_XpX_04hrs_matrix) <- sort(unique(design_XpX_04hrs$sex))
design_XpX_04hrs_matrix


## Visualise sample distances ##
plotMDS(DGEList_XpX)
counts_for_pca<-cpm(DGEList_XpX,log=TRUE,prior.count=1) 
pca_output <- pca(counts_for_pca, metadata = DGEList_XpX$samples)

pca_04hrs <- biplot(pca_output, colby = "group", gridlines.major = FALSE, lab = NULL, 
                        title = "A - X-linked genes PCA", titleLabSize = 20)  +
  scale_colour_manual(name = "Reproductive morph", values = c(
    "androgenic" = "chartreuse3",
    "gynogenic" = "purple"), labels = c(
      "androgenic" = "Androgenic",
      "gynogenic" = "Gynogenic")) +
  guides(colour = guide_legend(override.aes = list(size=4), nrow = 2)) +
  theme(plot.title.position = "plot",
        plot.title = element_text(size = 22, hjust = 0, vjust = 3.5),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_text(size = 16), legend.text = element_text(size = 16, margin = margin(t = 0, b = 0)),
        legend.spacing.y = unit(0.2, "cm"),
        legend.box.spacing = unit(0.2, "cm")) 
pca_04hrs


# dev.off()
screeplot(pca_output)

# Estimate common dispersion and tagwise dispersion in one go 
DGEList_XpX_04hrs <- estimateDisp(DGEList_XpX, design = design_XpX_04hrs_matrix)
# This is a plot of the NB dispersion for the dataset 
plotBCV(DGEList_XpX_04hrs)

# Given raw counts, NB dispersion(s) and a design matrix, glmQLFit() fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components.
fit_XpX_04hrs <- glmQLFit(DGEList_XpX_04hrs, design_XpX_04hrs_matrix)
head(fit_XpX_04hrs$coefficients)
# compares androgenics (as the baseline) to gynogenics 
qlf.androvsgyno_XpX_04hrs <- glmQLFTest(fit_XpX_04hrs, contrast=c(-1, 1))
qlf.androvsgyno_XpX_04hrs_FDR <- topTags(qlf.androvsgyno_XpX_04hrs, 2336)
# write.csv(qlf.androvsgyno_XpX_04hrs_FDR, "outputs/XpX_04hrs_androvsgyno_all.csv")
qlf.androvsgyno_XpX_04hrs_FDR_sig <- qlf.androvsgyno_XpX_04hrs_FDR$table[qlf.androvsgyno_XpX_04hrs_FDR$table$FDR < 0.05,]
summary(decideTests(qlf.androvsgyno_XpX_04hrs))
# write.csv(qlf.androvsgyno_XpX_04hrs_FDR_sig, "outputs/XpX_04hrs_androvsgyno_sig.csv")


### Volcano plots (which requires differential gene analysis) ###
# compares androgenics (as the baseline) to gynogenics 
summary(decideTests(qlf.androvsgyno_XpX_04hrs))
plotMD(qlf.androvsgyno_XpX_04hrs)
abline(h=c(-1,1), col="blue") 


qlf.androvsgyno_XpX_04hrs_FDR_df <- as.data.frame(qlf.androvsgyno_XpX_04hrs_FDR)

qlf.androvsgyno_XpX_04hrs_FDR_df <- qlf.androvsgyno_XpX_04hrs_FDR_df %>% 
  mutate(
    Expression = case_when(FDR < 0.05 & logFC <= 0 ~ "Androgenic-biased",
                           FDR < 0.05 & logFC >= 0 ~ "Gynogenic-biased",
                           TRUE ~ "Unbiased")
  )
table(qlf.androvsgyno_XpX_04hrs_FDR_df$Expression)

volcano_XpX_04hrs <- ggplot(qlf.androvsgyno_XpX_04hrs_FDR_df, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 2, alpha = 0.6) +
  scale_color_manual(values = c("chartreuse3", "purple", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept = -log(0.05,10), colour = "gray") +
  # geom_vline(xintercept = c(1, -1), colour="gray") +
  ggtitle("B - X-linked genes DGE") + theme_bw() +
  theme(plot.title.position = "plot",
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        axis.ticks = element_line(linewidth = 1.2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 22, hjust = 0.05, vjust = 2, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14), 
        legend.position = "bottom",
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        legend.spacing.y = unit(0.2, "cm"),
        legend.box.spacing = unit(0.2, "cm"))+
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2))
volcano_XpX_04hrs
```