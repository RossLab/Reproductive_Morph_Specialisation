# QC and Trim RNA-seq reads #

The quality of the reads were assessed before trimming with FastQC (v0.12.1), trimmed with Fastp (v0.24.0), then QC'ed again.

```
# Assess quality of reads with fastqc
fastqc -t 4 *.fq.gz
# -t 4 to use 4 threads

# Trim reads with fastp 
for file in $(ls B*/*_1.fq.gz)
do
        base=$(basename $file "_1.fq.gz")
        fastp -i B*/${base}_1.fq.gz \
        -I B*/${base}_2.fq.gz \
        -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz
done

# Assess quality of tirmmed reads with fastqc
fastqc -t 4 *.trimmed.fq.gz
```

The FastQC reports looked normal, so I can proceed to the next stage, which is mapping the reads to the genome.

# Map RNA-seq reads to the reference genome #

Map trimmed RNA reads to the genome with STAR (v2.7.11b). The genome used here is the Bcop_v3 genome which has the X' inversion scaffold as well. Ideally, reads coming from the X or the X' in gynogenic females should map to their gametolog of origin, and reads in androgenic females should only map to the X gametolog. However, since I am adding up the read count in each gametolog pair later on to compare between female morphs, this doesn't matter hugely. 
The other advantage of mapping to a reference that includes the inversion is that, since the androgenic females (XX) and males (X0) should not map to the inversion, I can compare the mapping rate of those samples to the X' to that of mapping to a reference without the X', to get an estimate of the mapping rate error. 

```
# Define the file name for the genome and for the annotation 
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3.augustus.gtf'

# genomeGenerate requires that a directory has been made beforehand for the genome index to go into
mkdir Bcop_v3-chromosomes.STAR

# run genomeGenerate to generate a genome index 
STAR \
--runThreadN 16 \ 
--runMode genomeGenerate \
--genomeSAindexNbases 12 \  
--outFileNamePrefix Bcop_v3-chromosomes \
--sjdbGTFfile ${ANNO} \  # STAR will extract splice junctions from the annotation, which improves mapping
--genomeDir Bcop_v3-chromosomes.STAR \
--genomeFastaFiles ${GENOME}

echo "aligning RNAseq reads"
for file in $(ls outputs/*_1.trimmed.fq.gz)
do
        base=$(basename $file "_combined_1.fq.gz")
        STAR \
        --runThreadN 16 \  # specifies the use of 16 threads
        --alignTranscriptsPerReadNmax 20000 \
        --alignEndsType EndToEnd \
        --outSAMattributes Standard \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn outputs/${base}_combined_1.fq.gz outputs/${base}_combined_2.fq.gz \
        --readFilesCommand zcat \
        --outTmpDir ${base}.out \
        --outFileNamePrefix ${base}.STAR. \
        --genomeDir Bcop_v3-chromosomes.STAR
done
```
Average mapping rate was ~93%, average unique mapping rate was ~83%. Females (both androgenic and gynogenic) tended to have higher unique mapping rates than males in all tissues, which may be due to higher TE expression in males? 

Looking only at androgenic females and males, comparing the % of unmapped reads when mapping to the genome with the X' scaffold, and mapping to the same genome with the X' scaffold excluded, provides an estimate of the mismapping rate. Since they do not have the X' chromosome, the difference in % of unmapped reads should be due to reads mismapping to the X' scaffold. The average difference in % of unmapped reads in the two mapping exercises is ~1.9% for males and ~1.2% for androgenic females, showing low mismapping rate. 

Once reads have been mapped to the reference genome, the number of reads can be summarised to the gene level with featureCounts.

# Read summarisation with featureCounts #

featureCounts from the subread package (v2.0.8) was used to summarise the mapped reads. 

```
# Specify the genome and annotation file
GENOME='Bcop_v3-chromosomes.fasta'
ANNO='Bcop_v3_geneid.augustus.gtf'

# Run featureCuounts
for file in $(ls outputs/*.STAR.Aligned.sortedByCoord.out.bam)
do
        echo "running featureCounts for" $file
        base=$(basename $file ".STAR.Aligned.sortedByCoord.out.bam")
        featureCounts -T 5 -p --countReadPairs -M --fraction -a ${ANNO} -t exon -g gene_id -o ${base}.multi.featureCounts.txt outputs/${base}.STAR.Aligned.sortedByCoord.out.bam
done
```

The `--countReadPairs -M --fraction` tag allows multimapping reads to be count, with the count distributed evenly between all sites that the read maps to. Therefore, if there are reads that map equally well to the X or the X' genome, they will get fractionally distributed between the two. Again, since I am adding up read counts from both gametologs in downstream analysis, this doesn't matter hugely. 

AFter all this, I have a txt file per sample with read counts, summarised at the gene level for each gene. This txt file can then be taken into R for differential gene expression analysis.
