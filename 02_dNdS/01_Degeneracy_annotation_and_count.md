# Annotating synonymous and nonsynonymous sites per gene #
This follows section 4.2. of [04_divergence.md](https://github.com/RossLab/PGE_slowerX/blob/main/04_divergence.md), to calculate the total number of generate vs degenerate sites in each gene. 
\
\
I am using `partitioncds.py` from [this github](https://github.com/A-J-F-Mackintosh/Mackintosh_et_al_2022_Ipod/tree/main). This script requires a reference genome, an annotation bed file, and a vcf. To generate a vcf, I align an Illumina library from a gynogenic female (with the inversion) to the reference genome (with the inversion) and call variants.



```
# bwa (v0.7.18)
# samtools (v1.21)
# picard (v3.3.0)
# bcftools v(1.21)

# Define file names for the Illumina librarys (KM1 and KM2), and the reference genome (GENOME). 
KM1='/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/KM_heterozygosity/inputs/KMG_1.trimmed.fq.gz'
KM2='/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/KM_heterozygosity/inputs/KMG_2.trimmed.fq.gz'
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
base='Bcop_v3_KM'

## Using bwa-mem to align Illumina library to reference genome
bwa index ${GENOME}
bwa mem -t 32 ${GENOME} ${KM1} ${KM2} | \
samtools view -b - > ${base}.bwa.bam && \
samtools sort -@ 8 -o ${base}.bwa.sorted.bam ${base}.bwa.bam
samtools flagstat ${base}.bwa.sorted.bam > ${base}.bwa.sorted.flagstat.txt
samtools index ${base}.bwa.sorted.bam

# mark and remove duplicates, and assign unique read groups
picard AddOrReplaceReadGroups INPUT=${base}.bwa.sorted.bam OUTPUT=${base}.bwa.rg.sorted.bam RGID=${base} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${base}
# Duplicates arise during sample preparation e.g. library construction using PCR
picard MarkDuplicates MAX_RECORDS_IN_RAM=400000 INPUT=${base}.bwa.rg.sorted.bam OUTPUT=${base}.bwa.rg.rmdup.sorted.bam M=${base}.met.txt REMOVE_DUPLICATES=true
samtools index ${base}.bwa.rg.rmdup.sorted.bam

## Use bcftools to call variants
# bcftools mpileup ouputs a raw pileup of sequencing depths and read base calls, to compute base call likelihoods
# bcftools call interprets the base call likelihoods to determine variants. -m allows calling of multiallelic sites, and -v outputs variant sites only. -Oz makes it output as a compressed vcf
bcftools mpileup -Ou -f ${GENOME} ${base}.bwa.rg.rmdup.sorted.bam | bcftools call -mv -Oz -o ${base}.vcf.gz
bcftools index ${base}.vcf.gz
```
\
I have a gff file for the genome but not a bed file, so I had to make that. The gff file was also missing info for some rows so had to fix it beforehand. 
```
# Define file names
GFF="/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3_fixed.gff3"
base="Bcop_v3_fixed"

cat ${GFF} | gt gff3 -sort -tidy -retainids -fixregionboundaries > ${base}.tidy.gff3
gff2bed < ${base}.tidy.gff3 > ${base}.tidy.gt.bed
cat ${base}.tidy.gt.bed | awk '$8=="CDS"' OFS="\t" | cut -f 1,2,3,5,6,10 | sort -k1,1V -k2,2n > ${base}.tidy.gt.sorted.bed
awk 'BEGIN {OFS="\t"} {
    split($6, arr, ";");
    split(arr[1], id_arr, "=");
    transcript_id = id_arr[2];
    print $1, $2, $3, transcript_id, $4, $5;
}' ${base}.tidy.gt.sorted.bed > ${base}.tidy.gt.sorted.input.bed
```
\
Once I have all the files I need, I can run `partitioncds.py`.
```
# Sync files in or define file names
GENOME="/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta"
BEDFILE="/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/output/Bcop_v3_fixed.tidy.gt.sorted.input.bed"
VCF="/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/output/Bcop_v3_KM.vcf.gz"
ANNO="/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3_fixed.gff3"

python /mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/scripts/partitioncds.py -f $GENOME -v $VCF -b $BEDFILE

### subset 0x, 2x, 3x and 4x degenerate sites:
grep "_0" gimble.cds.bed > gimble.cds.0xsites.bed
grep "_2" gimble.cds.bed > gimble.cds.2xsites.bed
grep "_3" gimble.cds.bed > gimble.cds.3xsites.bed
grep "_4" gimble.cds.bed > gimble.cds.4xsites.bed
gzip gimble.cds.bed

gff2bed < ${ANNO} > Bcop_v3_fixed.bed
ANNO_BED="Bcop_v3_fixed.bed"
### intersect degeneracy annotation with bedfile
bedtools intersect -a ${ANNO_BED} -b gimble.cds.0xsites.bed > genes.0x.sites.bed && rm gimble.cds.0xsites.bed
grep "mRNA" genes.0x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.0x.sites.fixed.bed && gzip genes.0x.sites.fixed.bed #&& rm genes.0x.sites.bed
bedtools intersect -a ${ANNO_BED} -b gimble.cds.2xsites.bed > genes.2x.sites.bed && rm gimble.cds.2xsites.bed
grep "mRNA" genes.2x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.2x.sites.fixed.bed && gzip genes.2x.sites.fixed.bed #&& rm genes.2x.sites.bed
bedtools intersect -a ${ANNO_BED} -b gimble.cds.3xsites.bed > genes.3x.sites.bed && rm gimble.cds.3xsites.bed
grep "mRNA" genes.3x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.3x.sites.fixed.bed && gzip genes.3x.sites.fixed.bed #&& rm genes.3x.sites.bed
bedtools intersect -a ${ANNO_BED} -b gimble.cds.4xsites.bed > genes.4x.sites.bed && rm gimble.cds.4xsites.bed
grep "mRNA" genes.4x.sites.bed | sed 's/ID=//g' | sed 's/;.*//g' | cut -f1,2,3,10 > genes.4x.sites.fixed.bed && gzip genes.4x.sites.fixed.bed #&& rm genes.4x.sites.bed
```
\
And summarise counts for each gene with this R script. 
```
# R

library(data.table)
library(dplyr)

degen_0_sites <- read.table('genes.0x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_0_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')
degen_2_sites <- read.table('genes.2x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_2_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')
degen_3_sites <- read.table('genes.3x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_3_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')
degen_4_sites <- read.table('genes.4x.sites.fixed.bed.gz', header=F, stringsAsFactors=F)
colnames(degen_4_sites) <- c('scaf', 'bed_start', 'bed_end', 'transcript')

degen_0_counts <- degen_0_sites %>% count(transcript)
colnames(degen_0_counts) <- c('transcript', '0x_sites')
degen_2_counts <- degen_2_sites %>% count(transcript)
colnames(degen_2_counts) <- c('transcript', '2x_sites')
degen_3_counts <- degen_3_sites %>% count(transcript)
colnames(degen_3_counts) <- c('transcript', '3x_sites')
degen_4_counts <- degen_4_sites %>% count(transcript)
colnames(degen_4_counts) <- c('transcript', '4x_sites')

degen_counts_02 <- merge(degen_0_counts, degen_2_counts, by=c('transcript'), all=T)
degen_counts_023 <- merge(degen_counts_02, degen_3_counts, by=c('transcript'), all=T)
all_degen_counts <- merge(degen_counts_023, degen_4_counts, by=c('transcript'), all=T)

# get N syn vs nonsyn sites
all_degen_counts$N_non <- (all_degen_counts$`0x_sites`)+(all_degen_counts$`2x_sites`*0.5)+(all_degen_counts$`3x_sites`*0.33)
all_degen_counts$N_syn <- (all_degen_counts$`4x_sites`)+(all_degen_counts$`2x_sites`*0.5)+(all_degen_counts$`3x_sites`*0.66)
nrow(all_degen_counts)
# write table
write.table(all_degen_counts, file='degen_site_counts.tsv', row.names=F, col.names=T, quote=F, sep='\t')
```
Now that I have the number of synonymous and nonsynonymous sites in each gene, I can move on to finding synonymous vs nonsynonymous mutations. 
