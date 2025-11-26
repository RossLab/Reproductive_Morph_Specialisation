## Map RNA-seq reads to the reference genome ##

Map trimmed RNA reads to the genome with STAR (v2.7.11b). The genome used here is the Bcop_v3 genome which has the X' inversion scaffold as well. Reads from androgenic (XX) females and males (X0) should not map to the X' scaffold, so the mapping rate of those samples to X' gives an estimate of mapping error rate. 

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
--genomeSAindexNbases 12 \  # Should be min(14, log2(GenomeLength)/2 - 1) for smaller genomes
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

Average mapping rate was ~93%, average unique mapping rate was ~83%. Females (both androgenic and gynogenic) tended to have higher unique mapping rates than males in all tissues. 

Looking only at androgenic females and males, comparing the % of unmapped reads when mapping to the genome with the X' scaffold, and mapping to the same genome with the X' scaffold excluded, provides an estimate of the mismapping rate. Since they do not have the X' chromosome, the difference in % of unmapped reads should be due to reads mismapping to the X' scaffold. The average difference in % of unmapped reads in the two mapping exercises is ~1.9% for males and ~1.2% for androgenic females, showing low mismapping rate. 

Once reads have been mapped to the reference genome, the number of reads can be summarised to the gene level.
