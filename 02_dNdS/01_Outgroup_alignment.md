# Outgroup alignment #
\
I want to compare the rate of evolution of homologous genes on the X vs the X'. To do that I align the _Bradysia coprophila_ genome to the closest relative for which we have a genome, _B. odoriphaga_. _B. odoriphaga_ does not have a X' inversion. I do this with minimap2 (v2.28), samtools (v1.21) and bedtools (v2.31.1).

```
# Define the file names for the focal species (B. coprophila) and the outgroup (B. odoriphaga).
FOCAL='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
OUTGROUP='/mnt/loki/ross/assemblies/flies/sciaridae/genomes_from_qm/Bradysia_odoriphaga/GCA_016920775.1_B2019C_genomic.fna'
NAME='Bcop_Bodo_minimap2'

## Align genomes using minimap2
# --secondary = no keeps only primary alignments
# -ax asm20 is optimized for aligning highly similar assemblies and whole genome comparisons
minimap2 -t 30 -c --secondary=no -ax asm20 $FOCAL $OUTGROUP > $NAME.sam
samtools view -bo $NAME.bam $NAME.sam && samtools sort -o $NAME.sorted.bam $NAME.bam && rm $NAME.bam

bedtools bamtobed -ed -i $NAME.sorted.bam > $NAME.sorted.bed && rm $NAME.sorted.bam

# Sync outputs out
rsync -av $NAME.sorted.bed /mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/output/
```

This produces a bed file which 
