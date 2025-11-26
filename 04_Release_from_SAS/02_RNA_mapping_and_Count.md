# Map RNA to N-masked genome and summarise read counts #

Once I have my genome, subsetted for only autosomes and X chromosome (aka no Inversion), and N-masked at variant sites on the X chromosome, I can align my reads to this N-masked genome.
\
\
I use the trimmed reads I have previously produced in 01_Gene_Expression_Divergence, and map them to the N-masked genome with STAR (v2.7.11b). 

```

# Define file name
NGENOME="outputs/Bcop_v3-chromosomes_X_Nmasked.fasta"
ANNO="Bcop_v3_geneid.augustus.gtf"
READLOC="outputs/"

awk '$1 == "X" || $1 == "II" || $1 == "III" || $1 == "IV"' ${ANNO} > Bcop_v3_geneid_noinv.augustus.gtf
ANNO="Bcop_v3_geneid_noinv.augustus.gtf"

mkdir Bcop_v3-chromosomes_X_Nmasked.STAR

# run genomeGenerate
echo "run genomeGenerate"
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeSAindexNbases 12 \
--outFileNamePrefix Bcop_v3-chromosomes_X_Nmasked.fasta \
--sjdbGTFfile ${ANNO} \
--genomeDir Bcop_v3-chromosomes_X_Nmasked.STAR \
--genomeFastaFiles ${NGENOME}

# map RNAseq reads to genome with X and X'.
echo "aligning RNAseq reads"
for file in ${READLOC}*_1.trimmed.fq.gz
do
        base=$(basename $file "_1.trimmed.fq.gz")
        STAR \
        --runThreadN 16 \
        --alignTranscriptsPerReadNmax 20000 \
        --alignEndsType EndToEnd \
        --outSAMattributes Standard \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${READLOC}${base}_1.trimmed.fq.gz ${READLOC}${base}_2.trimmed.fq.gz \
        --readFilesCommand zcat \
        --outTmpDir ${base}.out \
        --outFileNamePrefix ${base}.STAR. \
        --genomeDir Bcop_v3-chromosomes_X_Nmasked.STAR
done
```
\
I checked the `*.STAR.Log.final.out` for the % of uniquely mapped reads, multimapping reads, unmapped reads, etc. They look more or less comparable to mapping to the full genome with the Inversion scaffold. 
\
\
Then I summarise the reads with featureCounts.
```
# Define file names
NGENOME="outputs/Bcop_v3-chromosomes_X_Nmasked.fasta"
ANNO="Bcop_v3_geneid.augustus.gtf"

awk '$1 == "X" || $1 == "II" || $1 == "III" || $1 == "IV"' ${ANNO} > Bcop_v3_geneid_noinv.augustus.gtf
ANNO="Bcop_v3_geneid_noinv.augustus.gtf"

# STAR aligned bam files
BAMLOC="outputs/"

echo "running featureCounts"
for file in ${BAMLOC}*.out.bam
do
        echo "running featureCounts for" $file
        base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
        featureCounts -T 5 -p --countReadPairs -M --fraction -a ${ANNO} -t exon -g gene_id -o ${BAMLOC}${base}.multi.featureCounts.txt ${BAMLOC}${base}.STAR.Aligned.sortedByCoord.out.bam
done
```
