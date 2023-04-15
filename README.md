# Spiroplasma_ms

# Step 1: Retrieving novel *Spiroplasma* genomes from ticks
## 1.1. Adapters' trimming
```
atropos -T 4 -a file:$ech-adaptersF -A file:$ech-adaptersR -o $ech-R1-trimmed.fastq -p $ech-R2-trimmed.fastq -pe1 $ech-R1.fastq -pe2 $ech-R2.fastq
```
For `CAS1` and `CAS3`:
```
cutadapt -b file:CAS1-adaptersF -o CAS1-R1-trimmed.fastq.gz CAS1-R1.fastq.gz
```
