# Spiroplasma_ms

# Step 1: Retrieving novel *Spiroplasma* genomes from ticks
## 1.1. Adapters' trimming
```
atropos -T 4 -a file:24PF-adaptersF -A file:24PF-adaptersR -o 24PF-R1-trimmed.fastq.gz -p 24PF-R2-trimmed.fastq.gz -pe1 24PF-R1.fastq.gz -pe2 24PF-R2.fastq.gz
```
