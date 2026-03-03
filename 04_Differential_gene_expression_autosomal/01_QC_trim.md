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
