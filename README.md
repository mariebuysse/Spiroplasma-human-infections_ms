# Spiroplasma_ms

# Step 1: Retrieving novel *Spiroplasma* genomes from ticks
## 1.1. Adapters' trimming
```
gzip -d 24PF-R1.fastq.gz
gzip -d 24PF-R1.fastq.gz
atropos -T 4 -a file:24PF-adaptersF -A file:24PF-adaptersR -o 24PF-R1-trimmed.fastq -p 24PF-R2-trimmed.fastq -pe1 24PF-R1.fastq -pe2 24PF-R2.fastq
tar -czvf 24PF-R1-trimmed.fastq.gz 24PF-R1-trimmed.fastq
tar -czvf 24PF-R2-trimmed.fastq.gz 24PF-R2-trimmed.fastq
```
