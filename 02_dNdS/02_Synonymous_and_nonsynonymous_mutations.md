# Annotating synonymous and nonsynonymous variants per gene #
I need to find the number of synonymous and non-synonymous mutations that have accumulated in the X vs the Inversion. To do that I align the _Bradysia coprophila_ genome to the closest relative for which we have a genome, _B. odoriphaga_. _B. odoriphaga_ does not have a X' inversion. Importantly, I have to do this separately for two versions of the _B. coprophila_ genome: one with the autosomes (II, III, IV) plus the X chromosome, one with the autosomes plus the inversion (X'). The final dN/dS value for the autosomal genes should be the same for these two versions of the genome, but for the X and Inversion homologs they may be different. 
\
\
I first map PacBio long reads from _B. odoriphaga_ to the _B. coprophila_ genome with the autosome and X, then call variants. 
```
# samtools (v1.21)
# minimap2 (v2.28)
# picard (v3.3.0)
# bcftools (v1.21)

# Eefine file names, including the long reads for _B. odoriphaga_, and the reference genome for _B. coprophila_
Bodo='/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/input/SRR11366020_trimmed.fastq'
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
base='Bcop_v3_Bodo_X'

## Subset the genome for X, II, III, and IV chromosomes only
samtools faidx ${GENOME}
samtools faidx ${GENOME} II III IV X > Bcop-v3_subset_X.fasta

GENOME='Bcop-v3_subset_X.fasta'

## Use minimap2 to align pacbio long read to genome
minimap2 -t 30 -c --secondary=no -ax map-pb ${GENOME} ${Bodo} > ${base}.minimap2.sam
samtools view -b ${base}.minimap2.sam > ${base}.minimap2.bam && \
samtools sort -@ 8 -o ${base}.minimap2.sorted.bam ${base}.minimap2.bam
samtools flagstat ${base}.minimap2.sorted.bam > ${base}.minimap2.sorted.flagstat.txt
samtools index ${base}.minimap2.sorted.bam

# mark and remove duplicates
# assign unique read groups
picard AddOrReplaceReadGroups INPUT=${base}.minimap2.sorted.bam OUTPUT=${base}.minimap2.rg.sorted.bam RGID=${base} RGLB=lib1 RGPL=pacbio RGPU=unit1 RGSM=${base}
# Duplicates arise during sample preparation e.g. library construction using PCR
picard MarkDuplicates MAX_RECORDS_IN_RAM=400000 INPUT=${base}.minimap2.rg.sorted.bam OUTPUT=${base}.minimap2.rg.rmdup.sorted.bam M=${base}.met.txt REMOVE_DUPLICATES=true
samtools index ${base}.minimap2.rg.rmdup.sorted.bam

## Use bcftools to call variants
# genotype variants
# bcftools mpileup ouputs a raw pileup of sequencing depths and read base calls, to compute base call likelihoods
# bcftools call interprets the base call likelihoods to determine variants. -m allows calling of multiallelic sites, and -v outputs variant sites only. -Oz makes it output as a compressed vcf
bcftools mpileup -Ou -f ${GENOME} ${base}.minimap2.rg.rmdup.sorted.bam | bcftools call -mv -Oz -o ${base}.vcf.gz
bcftools index ${base}.vcf.gz
```
\
I do the same for the autosome and inversion version of the genome. 

```
# Define file names
Bodo='/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/input/SRR11366020_trimmed.fastq'
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
base='Bcop_v3_Bodo_Inv'

## Subset the genome for Inversion, II, III, and IV chromosomes only
samtools faidx ${GENOME}
samtools faidx ${GENOME} II III IV Inversion > Bcop-v3_subset_Inv.fasta

GENOME='Bcop-v3_subset_Inv.fasta'

## Using minimap2 to align pacbio long read to genome
minimap2 -t 30 -c --secondary=no -ax map-pb ${GENOME} ${Bodo} > ${base}.minimap2.sam
samtools view -b ${base}.minimap2.sam > ${base}.minimap2.bam && \
samtools sort -@ 8 -o ${base}.minimap2.sorted.bam ${base}.minimap2.bam
samtools flagstat ${base}.minimap2.sorted.bam > ${base}.minimap2.sorted.flagstat.txt
samtools index ${base}.minimap2.sorted.bam

# mark and remove duplicates
# assign unique read groups
picard AddOrReplaceReadGroups INPUT=${base}.minimap2.sorted.bam OUTPUT=${base}.minimap2.rg.sorted.bam RGID=${base} RGLB=lib1 RGPL=pacbio RGPU=unit1 RGSM=${base}
# Duplicates arise during sample preparation e.g. library construction using PCR
picard MarkDuplicates MAX_RECORDS_IN_RAM=400000 INPUT=${base}.minimap2.rg.sorted.bam OUTPUT=${base}.minimap2.rg.rmdup.sorted.bam M=${base}.met.txt REMOVE_DUPLICATES=true
samtools index ${base}.minimap2.rg.rmdup.sorted.bam

## Use bcftools to call variants
# genotype variants
# bcftools mpileup ouputs a raw pileup of sequencing depths and read base calls, to compute base call likelihoods
# bcftools call interprets the base call likelihoods to determine variants. -m allows calling of multiallelic sites, and -v outputs variant sites only. -Oz makes it output as a compressed vcf
bcftools mpileup -Ou -f ${GENOME} ${base}.minimap2.rg.rmdup.sorted.bam | bcftools call -mv -Oz -o ${base}.vcf.gz
bcftools index ${base}.vcf.gz
```
\
Then, I can use snpEff to annotate the vcf with the type of mutation the SNP produces (missense, nonsense, etc), and filter only for variants_effect_missense_variant (nonsynonymous), variants_effect_stop_retained_variant and variants_effect_synonymous_variant (synonymous). Again, I do this separately for the autosome + X and autosome + Inversion vcf. 
```
# gffread v0.12.7
# snpeff v5.2
# snpsift v5.2 

# Define file names
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
ANNO='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3.augustus.gtf'
VCF='/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/output/Bcop_v3_Bodo_X.vcf.gz'

# I have edited the config file of snpEff in my conda environment, which is ~/micromamba/envs/vcf/share/snpeff-5.2-1/snpEff.config to add Bradysia coprophila
# Bradysia coprophila genome (v3)
# Bradysia_coprophila.genome : Bradysia_coprophila

# I have then created a data/ folder in the same directory. Within that I have created a Bradysia_coprophila directory and copied in both the ref genome and the gtf annotation

gffread -x cds.fa -g ${GENOME} ${ANNO}
mv cds.fa ~/micromamba/envs/vcf/share/snpeff-5.2-1/data/Bradysia_coprophila/

gffread -y protein.fa -g ${GENOME} ${ANNO}
mv protein.fa ~/micromamba/envs/vcf/share/snpeff-5.2-1/data/Bradysia_coprophila/

# Build SnpEff database
snpEff build -gtf22 -v Bradysia_coprophila

# Annotate the vcf with snpEff
snpEff Bradysia_coprophila \
-v -no-downstream -no-intron -no-upstream -no-utr -no-intergenic \
-stats Bcop.filterstats.html \
-canon \
${VCF} > Bcop_vs_Bodo_anno.vcf

# Extract useful reasons with snpSift
SnpSift extractFields Bcop_vs_Bodo_anno.vcf \
CHROM ANN[0].GENEID POS REF ALT isHom"(GEN[0])" ANN[*].EFFECT \
> Bcop.just_genes_snps.txt

cat Bcop.filterstats.genes.txt | cut -f2,3,18,26,27 > Bcop_X.filterstats.genes.relevant.txt
```

```
# Define file names
GENOME='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta'
ANNO='/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3.augustus.gtf'
VCF='/mnt/loki/ross/flies/sciaridae/Bradysia_coprophila/B_coprophila_morph_gene_divergence/04_dnds/output/Bcop_v3_Bodo_Inv.vcf.gz'

# I have edited the config file of snpEff in my conda environment, which is ~/micromamba/envs/vcf/share/snpeff-5.2-1/snpEff.config to add Bradysia coprophila
# Bradysia coprophila genome (v3)
# Bradysia_coprophila.genome : Bradysia_coprophila

# I have then created a data/ folder in the same directory. Within that I have created a Bradysia_coprophila directory and copied in both the ref genome and the gtf annotation

gffread -x cds.fa -g ${GENOME} ${ANNO}
mv cds.fa ~/micromamba/envs/vcf/share/snpeff-5.2-1/data/Bradysia_coprophila/

gffread -y protein.fa -g ${GENOME} ${ANNO}
mv protein.fa ~/micromamba/envs/vcf/share/snpeff-5.2-1/data/Bradysia_coprophila/

# Build SnpEff database
snpEff build -gtf22 -v Bradysia_coprophila

# Annotate the vcf with snpEff
snpEff Bradysia_coprophila \
-v -no-downstream -no-intron -no-upstream -no-utr -no-intergenic \
-stats Bcop.filterstats.html \
-canon \
${VCF} > Bcop_vs_Bodo_anno.vcf

# Extract useful reasons with snpSift
SnpSift extractFields Bcop_vs_Bodo_anno.vcf \
CHROM ANN[0].GENEID POS REF ALT isHom"(GEN[0])" ANN[*].EFFECT \
> Bcop.just_genes_snps.txt

cat Bcop.filterstats.genes.txt | cut -f2,3,18,26,27 > Bcop_Inv.filterstats.genes.relevant.txt
```
\
Now I have the number of synonymous and nonsynonymous mutations for each gene on the autosomes, X chromosome and the inversion. To allow a comparison between the two, I run Orthofinder to find 1-1 homologs between the X chromosome and the inversion. 
```
# Define file names
GFF="/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3_fixed.gff3"
GENOME="/mnt/loki/ross/assemblies/flies/sciaridae/Bradysia_coprophila/Bcop_v3-chromosomes.fasta"

awk '$1 == "X" && $3 == "mRNA"' ${GFF} > Bcop_X.gff
awk '$1 == "Inversion" && $3 == "mRNA"' ${GFF} > Bcop_Inv.gff

gffread -w Bcop_X_genes.fasta -g ${GENOME} Bcop_X.gff
gffread -w Bcop_Inv_genes.fasta -g ${GENOME} Bcop_Inv.gff

mkdir genes/

OUT='X_vs_Inv'
mv Bcop_X_genes.fasta genes/
mv Bcop_Inv_genes.fasta genes/
orthofinder -o $OUT -n $OUT -t 12 -d -f genes/
mv $OUT/Results_$OUT/Orthogroups/Orthogroups.tsv ${OUT}_Orthogroups.tsv
mv $OUT/Results_$OUT/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ${OUT}_Orthogroups_SCOs.tsv
```
