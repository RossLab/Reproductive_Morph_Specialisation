I first map the RNA-seq reads to the _B. coprophila_ genome. (Mention that the genome and the RNA seq reads came from the same inbred lab line, but that i am only using reads from gynogenic females) Since this is not exactly allele-specific expression, but more like homolog specific expression, I am using strict mapping rules to ensure that reads are mapping to the correct homolog on the X vs the X'. This involves keeping only primary alignments in the case of imperfect multimappers, and filtering away perfect multimappers. In this case I should be left with only informative reads. 

```
# Trim RNA reads to retain only high quality reads to minimising sequencing errors
# Default fastp quality threshold is phred quality of >=15. I'm going to bump that up to 20.
for file in $(cat gyno_list.txt)
do
        base=$(basename $file "_1.fq.gz")
        fastp -q 20 --correction -i /mnt/hel/ross/sequencing/raw/novogene_may_2023_gnats/X204SC23040157-Z01-F001_02/01.RawData/B*/${base}_1.fq.gz \
        -I /mnt/hel/ross/sequencing/raw/novogene_may_2023_gnats/X204SC23040157-Z01-F001_02/01.RawData/B*/${base}_2.fq.gz \
        -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz
done

# Define file names for RNA mapping
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
ANNO='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3.augustus.gtf'

#mkdir Bcop_v3-chromosomes.STAR

# run genomeGenerate
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
for file in ${SCRATCH}/*_1.trimmed.fq.gz
do
        echo "Processing file: $file"
        base=$(basename $file "_1.trimmed.fq.gz")
        STAR \
        --runThreadN 16 \
        --alignTranscriptsPerReadNmax 20000 \
        --outFilterMultimapNmax 10 \
        --alignEndsType EndToEnd \
        --outSAMattributes Standard \
        --outSAMprimaryFlag AllBestScore \
        --outSAMmultNmax 1 \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${base}_1.trimmed.fq.gz ${base}_2.trimmed.fq.gz \
        --readFilesCommand zcat \
        --outTmpDir ${base}.out \
        --outFileNamePrefix ${base}.STAR. \
        --genomeDir Bcop_v3-chromosomes.STAR
done
# AllBestScore outputs all alignments with the best score as primary alignments i.e. perfect multimappers, so I can filter them later
# outSAMmultiNmax 1 only prints the alignment with the highest score

echo "filtering BAM files"
for file in ${SCRATCH}/*.STAR.Aligned.sortedByCoord.out.bam
do
        base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
        samtools view -b -q 10 ${base}.STAR.Aligned.sortedByCoord.out.bam > ${base}.STAR.Aligned.sortedByCoord.filtered.out.bam
done
# This filters out perfect multimappers, which will have a very low mapping quality score. 
```
This leaves me with a bam file of reads that are not perfect multimappers, or the primary alignment in imperfect multimappers. Now I can summarise the read count with featureCounts. 
```
# Define file names
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
ANNO='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3_geneid.augustus.gtf'

echo "running featureCounts"
for file in $(ls *.STAR.Aligned.sortedByCoord.filtered.out.bam)
do
        echo "running featureCounts for" $file
        base=$(basename $file ".STAR.Aligned.sortedByCoord.filtered.out.bam")
        featureCounts -T 5 -p --countReadPairs -M --primary -a ${ANNO} -t exon -g gene_id -o ${base}.featureCounts.txt ${base}.STAR.Aligned.sortedByCoord.filtered.out.bam
done
```
I then move to R to carry out homolog-specific expression. 
