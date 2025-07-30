# Creating an N-masked reference genome #
To be able to compare gynogenic and androgenic expression to males, I needed to map all reads to the same genome with the same X chromosome. However, there are sequence differences between the X and X' homologs in gynogenic females. Therefore, I first need to create an N-masked genome, which is a genome that is masked at the sites with variants between the X and X' homolog. This is adapted from aspects of Baird 2024 (thesis), Chapter 4. 
\
To keep things consistent, I want to use Bcop v3, which is the _Bradysia coprophila_ reference genome that I've been using in all my other analyses. However it has the X', so I will just subset for the autosomes and X for the purpose of making this N-masked genome. 
\
To figure out the variants between the X and the X', I need to map genomic reads with the X' haplotype (aka a gynogenic female or pooled gynogenic females) to the reference genome without the X'. In this instance I am using pooled data from 50 X'X females. I am using fastp (v 0.24.0) to trim the reads. 
```
# combine files
cat LR43_EDSW200011440-1a_HJ5JVDSXY_L2_1.fq.gz LR43_EDSW200011440-1a_HJ5FNDSXY_L2_1.fq.gz > H2_XpX_pool_reseq_1.fq.gz
cat LR43_EDSW200011440-1a_HJ5JVDSXY_L2_2.fq.gz LR43_EDSW200011440-1a_HJ5FNDSXY_L2_2.fq.gz > H2_XpX_pool_reseq_2.fq.gz

rm LR*

# trim
for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
  	fastp -i ${base}_1.fq.gz -I ${base}_2.fq.gz -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz
done
```
I then map it to the reference genome. As I said before, I need to subset for the autosomes and X chromosome only. I then use bwa to map the reads to the reference genome since it is meant to be better than minimap2 for short read (see https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa).
```
# Define file names
GENOME="/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta"
READ_1="/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence_FINAL/06_sexually_antagonistic_selection/N_masked_genome/outputs/H2_XpX_pool_1.trimmed.fq.gz"
READ_2="/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence_FINAL/06_sexually_antagonistic_selection/N_masked_genome/outputs/H2_XpX_pool_2.trimmed.fq.gz"

## Subset the genome for X, II, III, and IV chromosomes only
samtools faidx ${GENOME}
samtools faidx ${GENOME} II III IV X > Bcop-v3_subset_X.fasta

bwa index $GENOME

bwa mem -t 32 $GENOME ${READ_1} ${READ_2} | \
samtools view -b - > H2_XpX_pool.bam && \
samtools sort -@ 32 -o H2_XpX_pool.sorted.bam H2_XpX_pool.bam
```
I prefer to call variants with bcftools. I first mark duplicates with samtools to reduce the chance of false positives.
```
# Define file names
BAM="/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence_FINAL/06_sexually_antagonistic_selection/N_masked_genome/outputs/H2_XpX_pool.sorted.bam"
GENOME="/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta"

samtools index ${BAM}

# Use Picard to mark/remove duplicates
picard AddOrReplaceReadGroups \
    I=${BAM} \
    O=H2_XpX_pool.rg.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=sample1 \
    VALIDATION_STRINGENCY=SILENT

picard MarkDuplicates \
    I=H2_XpX_pool.rg.bam \
    O=H2_XpX_pool.rmdup.bam \
    M=H2_XpX_pool.dup_metrics.txt \
    REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=SILENT

samtools index H2_XpX_pool.rmdup.bam

# Call variants
bcftools mpileup --threads 16 -Ou -f $GENOME H2_XpX_pool.rmdup.bam | \
bcftools call --threads 16 -m -A -Oz -o H2_XpX_pool.vcf.gz

bcftools index H2_XpX_pool.vcf.gz

# Filter for variants with high quality
bcftools filter --threads 16 -i 'QUAL>=30 && DP>=10 && MQ>=40 && TYPE="snp"' H2_XpX_pool.vcf.gz -Oz -o H2_XpX_pool_filtered_snps.vcf.gz

bcftools index H2_XpX_pool_filtered_snps.vcf.gz
```
