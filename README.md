# General command lines and analytic steps used in the production of the *Spiroplasma-human-cases* manuscript

# Step 1: Retrieving *Spiroplasma* genomes from biological samples
## 1.1. Adapters' trimming
For `CAS1` and `CAS3`:
```
cutadapt -b file:CAS1-adaptersF -o CAS1-R1-trimmed.fastq.gz CAS1-R1.fastq.gz
```


