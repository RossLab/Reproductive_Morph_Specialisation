# Read summarisation with featureCounts #

featureCounts from the subread package (v2.0.8) was used to summarise the mapped reads. 

```
# Specify the genome and annotation file
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
ANNO='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3_geneid.augustus.gtf'

# Run featureCuounts
for file in $(ls /mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/01_STAR/outputs/*.STAR.Aligned.sortedByCoord.out.bam)
do
        echo "running featureCounts for" $file
        base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
        featureCounts -T 5 -p --countReadPairs -M --fraction -a ${ANNO} -t exon -g gene_id -o ${base}.multi.featureCounts.txt /mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/01_STAR/outputs/${base}.STAR.Aligned.sortedByCoord.out.bam
done
```

The `--countReadPairs -M --fraction` tag allows multimapping reads to be count, with the count distributed evenly between all sites that the read maps to. This produces a txt file with read counts, summarised at the gene level for each gene. This txt file can then be taken into R for differential gene expression analysis.
