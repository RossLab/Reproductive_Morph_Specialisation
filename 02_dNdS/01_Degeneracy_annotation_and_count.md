# Annotating synonymous and nonsynonymous sites per gene #
This follows section 4.2. of [04_divergence.md](https://github.com/RossLab/PGE_slowerX/blob/main/04_divergence.md), to calculate the total number of generate vs degenerate sites in each gene. 
/
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
