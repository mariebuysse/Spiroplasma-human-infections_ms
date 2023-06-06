### Analytic steps and general command lines used in the production of the *Spiroplasma-human-cases* manuscript

# Step 1: Retrieving *Spiroplasma* genomes from biological samples
In this study, we aimed to retrieve and analyze the *Spiroplasma* MAGs (Metagenome-Assembled Genome) from biological samples collected from two symptomatic infants with ocular infections. These cases have been reported in ***Spiroplasma* species as a rare cause of congenital cataract and uveitis: a case series** from **Farassat N, Reich M, Serr A, Küchlin S, Erwemi M, Auw-Hädrich C, Krastel H, and Lagrèze WA (BMC Ophthalmol, 2021)**. 

For each sample, a set of paired-end reads is available. Details about the experimental and sequencing methods are available in our related manuscript. Following analyses are based on these reads' sets, and those are referred to as:
```
ech="CASE1 CASE3"
```

## 1.1. Adapters' trimming
Adapters' sequences were removed from raw reads using `Cutadapt` (https://cutadapt.readthedocs.io/en/stable/, Cutadapt removes adapter sequences from high-throughput sequencing reads, 
Martin M,  EMBnet.journal, 2011, doi: org/10.14806/ej.17.1.200).
```
cutadapt -b file:$ech-adaptersF -o $ech-R1-trimmed.fastq.gz $ech-R1.fastq.gz
```

## 1.2. Assembly
Reads were assembled using `MEGAHIT` (https://github.com/voutcn/megahit, MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph, Li D, Liu C-M, Luo R, Sadakane K, and Lam T-W, Bioinformatics, 2015, doi: 10.1093/bioinformatics/btv033). 
```
command line
```

## 1.3. Binning and retrieving *Spiroplasma* MAGs
*Spiroplasma* MAGs were retrieved from assemblies using `CONCOCT` (https://github.com/BinPro/CONCOCT, Binning metagenomic contigs by coverage and composition, Alneberg J, Smári Bjarnason B, de Bruijn I, Schirmer M, Quick J, Ijaz U, Lahti L, Loman N, Andersson A, and Quince C, Nature Methods, 2014, doi: 10.1038/nmeth.3103) and the `anvi'o` pipeline (https://anvio.org/, Community-led, integrated, reproducible multi-omics with anvi’o, Eren A, Kiefl E, Shaiber A et al., Nature Microbiology, 2021, doi: 10.1038/s41564-020-00834-3). 

First, the contigs were formated to match the requirements of `anvi'o`. 
```
command line
```
The contigs were binned using `CONCOCT`: 
```
command line
```
The `clustering_merged` CSV file produced by CONCOCT needed to be exported in a tabular-delimited TEXT file. 
The bins' names had to be renamed to match the requirements of `anvi'o` : 
```
command line
```
The `bins-renamed.txt` file had to be transformed again to correspond to a tabular-delimited TEXT file, called `bins.txt` hereafter.

The bins were taxonomically assigned using the `anvi'o` pipeline below: 
```
command line
```
Let's see the `$ech-SUMMARY` file (html format) to check the taxnonomy results and some stats about the bins.
