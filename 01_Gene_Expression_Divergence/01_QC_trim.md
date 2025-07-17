# QC and Trim RNA-seq reads #

The quality of the reads were assessed before trimming with FastQC (v0.12.1), trimmed with Fastp (v0.24.0), then QC'ed again. This is done with a bash script.

```
# Assess quality of reads with fastqc
fastqc -t 4 *.fq.gz
# -t 4 to use 4 threads

# Trim reads with fastp 
for file in $(ls /mnt/hel/ross/sequencing/raw/novogene_may_2023_gnats/X204SC23040157-Z01-F001_0*/01.RawData/B*/*_1.fq.gz)
do
        base=$(basename $file "_1.fq.gz")
        fastp -i /mnt/hel/ross/sequencing/raw/novogene_may_2023_gnats/X204SC23040157-Z01-F001_0*/01.RawData/B*/${base}_1.fq.gz \
        -I /mnt/hel/ross/sequencing/raw/novogene_may_2023_gnats/X204SC23040157-Z01-F001_0*/01.RawData/B*/${base}_2.fq.gz \
        -o ${base}_1.trimmed.fq.gz -O ${base}_2.trimmed.fq.gz
done

# Assess quality of tirmmed reads with fastqc
fastqc -t 4 *.trimmed.fq.gz
```

The FastQC reports looked normal, so I can proceed to the next stage, which is mapping the reads to the genome.
