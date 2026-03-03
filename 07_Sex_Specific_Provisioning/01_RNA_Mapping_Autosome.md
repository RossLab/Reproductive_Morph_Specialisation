# Map embryo RNA-seq libraries to ref genome, and summarise read counts #
To assess whether gynogenic and androgenic females deposit different RNA into their eggs as a method of sex-specific provisioning, I carried out the same analysis as in 01_Gene_Expression_Divergence and 03_Homolog_specific_expression. 
\
\
First I will look at differentially deposited autosomal genes. I download the RNA-seq libraries from Genbank using `sra-tools` (v 3.2.0) and map it to the Bcop_v3 genome which has the X' inverison scaffold. 
```
FILE="Maternal_Deposit_RNAseq_SRAruns.txt"

# Set up directories
SRA_DIR="$SCRATCH/sra"
FASTQ_DIR="$SCRATCH/sra_fastq"
mkdir -p $SRA_DIR $FASTQ_DIR logs

while read ACCESSION; do
    echo "Processing: $ACCESSION"
    # Download the SRA file
    prefetch $ACCESSION  --max-size 30GB --output-directory $SRA_DIR
    # Convert SRA to FASTQ
    fasterq-dump --split-3 --threads 4 --outdir $FASTQ_DIR $SRA_DIR/$ACCESSION/$ACCESSION.sra
    # Compress FASTQ files
    gzip $FASTQ_DIR/${ACCESSION}_*.fastq
    echo "Completed: $ACCESSION"
done < ${FILE}

echo "All accessions processed!"
```
The quality of the reads were assessed before trimming with FastQC (v0.12.1), trimmed with Fastp (v0.24.0), then QC'ed again. 
```
# fastqc pre-trim
fastqc -t 4 inputs/*.fastq.gz

# trim reads
for file in $(ls inputs/*_1.fastq.gz)
do
        base=$(basename $file "_1.fastq.gz")
        fastp -i inputs/${base}_1.fastq.gz \
        -I inputs/${base}_2.fastq.gz \
        -o ${base}_1.trimmed.fastq.gz -O ${base}_2.trimmed.fastq.gz
done

# fastqc post-trim
fastqc -t 4 *.trimmed.fastq.gz
```
Map trimmed RNA reads to the genome with STAR (v2.7.11b)
```
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3.augustus.gtf'
LOC='inputs/'

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
featureCounts from the subread package (v2.0.8) was used to summarise the mapped reads.
```
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3_geneid.augustus.gtf'
LOC='outputs/'

echo "running featureCounts"
for file in $(ls ${LOC}*.STAR.Aligned.sortedByCoord.out.bam)
do
        echo "running featureCounts for" $file
        base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
        featureCounts -T 5 -p --countReadPairs -M --fraction -a ${ANNO} -t exon -g gene_id -o ${base}.multi.featureCounts.txt ${LOC}${base}.STAR.Aligned.sortedByCoord.out.bam
done
```
