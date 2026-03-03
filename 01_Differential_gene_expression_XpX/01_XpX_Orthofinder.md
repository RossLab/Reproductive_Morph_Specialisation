# Orthofinder #

I first use Orthofinder (v3.0.1b1) to find single copy orthologs between the X and X' inversion haplotypes. Beforehand, I have filtered my gff to include only the longest isoform of each gene. 
```
# Sync files in or define file names
GFF="Bcop_v3_fixed_longest_isoform.gff3"
GENOME="Bcop_v3-chromosomes.fasta"

awk '$1 == "X"' ${GFF} > Bcop_X.gff
awk '$1 == "Inversion"' ${GFF} > Bcop_Inv.gff

gffread -w Bcop_X_genes.fasta -g ${GENOME} Bcop_X.gff
gffread -w Bcop_Inv_genes.fasta -g ${GENOME} Bcop_Inv.gff

mkdir genes/

OUT='X_vs_Inv_longest'
mv Bcop_X_genes.fasta genes/
mv Bcop_Inv_genes.fasta genes/
orthofinder -o $OUT -n $OUT -t 12 -d -f genes/
mv $OUT/Results_$OUT/Orthogroups/Orthogroups.tsv ${OUT}_Orthogroups.tsv
mv $OUT/Results_$OUT/Orthogroups/Orthogroups_SingleCopyOrthologues.txt ${OUT}_Orthogroups_SCOs.tsv
```
Having this list of gametologs between the two haplotypes allows me to add up their read count in gynogenic females, so that I can compare X-linked gene expression between gynogenic and androgenic females.
