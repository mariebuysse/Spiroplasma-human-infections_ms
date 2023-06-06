### General command lines and analytic steps used in the production of the *Spiroplasma-human-cases* manuscript

# Step 1: Retrieving *Spiroplasma* genomes from biological samples
## 1.1. Adapters' trimming
```
cutadapt -b file:$ech-adaptersF -o $ech-R1-trimmed.fastq.gz $ech-R1.fastq.gz
```


