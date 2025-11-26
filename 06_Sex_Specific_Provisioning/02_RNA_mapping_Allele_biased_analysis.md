I am interested in whether there is biased deposition of transcripts from the X or X' gametolog in the eggs of gynogenic females. To study this I first map the RNA reads again, but with stricter filters for read quality and mapping quality. 
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
\
Now I have read counts for both autosomal genes between gynogenic and androgenic females, and X'X gametologs within gynogenic females. I can move on to differential gene expression analysis or allele-biased expression analysis. 
