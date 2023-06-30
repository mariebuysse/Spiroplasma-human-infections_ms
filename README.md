### Analytic steps and general command lines used in the production of the *Spiroplasma-human-cases* manuscript

# Step 1: Retrieving *Spiroplasma* genomes from biological samples
In this study, we aimed to retrieve and analyze the *Spiroplasma* MAGs (Metagenome-Assembled Genome) from biological samples collected from two symptomatic infants with ocular infections (i.e. lensectomy and vitrectomy samples). These cases have been reported in ***Spiroplasma* species as a rare cause of congenital cataract and uveitis: a case series** from **Farassat N, Reich M, Serr A, Küchlin S, Erwemi M, Auw-Hädrich C, Krastel H, and Lagrèze WA (BMC Ophthalmol, 2021, doi: 10.1186/s12886-021-02201-0)**. 

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
module load bioinfo/MEGAHIT/1.2.9
megahit -1 $ech-R1-trimmed.fastq.gz -2 $ech-R2-trimmed.fastq.gz -o $ech-metaMEGAHIT 
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
The `clustering_merged` CSV file produced by `CONCOCT` needed to be exported in a tabular-delimited TEXT file. 
The bins' names had to be renamed to match the requirements of `anvi'o` : 
```
command line
```
The `bins-renamed.txt` file had to be transformed again to correspond to a tabular-delimited TEXT file, called `bins.txt` hereafter.

The bins were taxonomically assigned using the `anvi'o` pipeline below: 
```
command line
```
Let's see the `$ech-SUMMARY` file to check the taxnonomy results and some stats about each bin.

## 1.4. Quality check of the final *Spiroplasma* MAGs
Quality and multiple statistics were accessed using miComplete (https://pypi.org/project/micomplete/, Hugoson E., Lam W.T., Guy L. (2020) miComplete: weighted quality evaluation of assembled microbial genomes. Bioinformatics. doi: 10.1093/bioinformatics/btz664) and Quast(https://github.com/ablab/quast, Gurevich A., Saveliev V., Vyahhi N., Tesler G. (2013) QUAST: quality assessment tool for genome assemblies. Bioinformatics. doi: 10.1093/bioinformatics/btt086).
```
## With Quast
quast.py ./Spiro-CASES-genomes/*.fasta -o ./RESULTS-QUAST

## With miComplete
miComplete set.tab --hmms Bact105 #set.tab a tabular separated file containing per line both a path to a genome and the type (i.e. fna)
```
Results are:
```
## Quast
Assembly                    Spiro-CAS1-contigs  Spiro-CAS3-contigs
# contigs (>= 0 bp)         195                 229               
# contigs (>= 1000 bp)      195                 229               
# contigs (>= 5000 bp)      86                  95                
# contigs (>= 10000 bp)     42                  45                
# contigs (>= 25000 bp)     5                   2                 
# contigs (>= 50000 bp)     0                   0                 
Total length (>= 0 bp)      1292297             1359048           
Total length (>= 1000 bp)   1292297             1359048           
Total length (>= 5000 bp)   1039119             1037392           
Total length (>= 10000 bp)  721044              683917            
Total length (>= 25000 bp)  180050              55248             
Total length (>= 50000 bp)  0                   0                 
# contigs                   195                 229               
Largest contig              40341               28125             
Total length                1292297             1359048           
GC (%)                      23.01               24.58             
N50                         11017               10069             
N75                         5947                5176              
L50                         35                  45                
L75                         74                  92                
# N's per 100 kbp           0.00                0.00

## miComplete
```

# Step 2: Genomes' description and comparison with others *Spiroplasma* genomes 
## 2.1. Phylogeny-based taxonomic assignation 
### 2.1.1. Based on *Spiroplasma* MAGs and whole genomes
First, single-copy orthologs (SCO) were identified using `OrthoFinder` (https://github.com/davidemms/OrthoFinder, OrthoFinder: phylogenetic orthology inference for comparative genomics, Emms DM and Kelly S, Genome Biology, 2019, doi: 10.1186/s13059-019-1832-y) from a set of specimens' genomes chosen to study the phylogenetic relationships between the obtained *Spiroplasma* MAGs and other *Spiroplasma* representatives:
```
commmand line
```
For each SCO, sequences were individually aligned using `mafft` (https://github.com/GSLBiotech/mafft, MAFFT multiple sequence alignment software version 7: improvements in performance and usability, Katoh K and Standley DM, Molecular Biology and Evolution, 2013 doi: 10.1093/molbev/mst010):
```
command line
```
For each SCO, ambigious hypervariable regions were removed using `trimAl` (https://github.com/inab/trimal, trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses, Capella-Gutiérrez S, Silla-Martínez JM, and Gabaldón T, Bioinformatics, 2009, doi: 10.1093/bioinformatics/btp348):
```
command line
```
Then, all SCO sequences were concatenated using `Amas` (https://github.com/marekborowiec/AMAS, AMAS: a fast tool for alignment manipulation and computing of summary statistics, Borowiec ML, PeerJ, 2016, doi: 10.7717/peerj.1660) in a single file:
```
command line
```
Substitution models were evaluated using `modeltest-ng` (https://github.com/ddarriba/modeltest, ModelTest-NG: A new and scalable tool for the selection of DNA and protein evolutionary models, Darriba D, Posada D, Kozlov AM, Stamatakis A, Morel B, and Flouri T, Molecular Biology and Evolution, 2020, doi: 10.1093/molbev/msz189) in order to determinate the appropriate substitution model (according to AICc criterion) to use for the phylogenetic tree construction with `RAxML` (https://github.com/stamatak/standard-RAxML, RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies, Stamatakis A, Bioinformatics, 2014, doi: 10.1093/bioinformatics/btu033):
```
command line
```
The phylogenetic tree was visualized and modified using `figtree` (https://github.com/rambaut/figtree/).

### 2.1.2. Based on single gene (16S rDNA gene)
Single-gene phylogenies were produced using the 16S rDNA gene sequences. This gene was used to describe *Spiroplasma* strains associated with previously reported cases of ocular infection. As different fragment lengths of the 16S rDNA gene were sequenced in those studies, several phylogenies will be produced to include the diversity of *Spiroplasma* strains. Sequences from other *Spiroplasma* representatives were retrieved from GenBank (National Center for Biotechnology Information), while sequences from published genomes were obtained using a local `BLAST` (https://www.ncbi.nlm.nih.gov/books/NBK279690/, Basic local alignment search tool, Altschul SF, Gish W, Miller W, Myers EW, and Lipman DJ, Journal of Molecular Biology, 1990, doi: 10.1016/S0022-2836(05)80360-2) and the following command lines:
```
command line
```

The phylogenetic trees were produced in a similar way that previously described (with adaptations to nucleotidic sequences):
```
command line
```

### 2.1.3. Multi-locus sequencing typing (MLST)
This phylogeny is based on an alignement of concatenated sequences of 4 genes (16S rDNA, *dnaK*, *gyrA*, and *rpoB* genes) from *Spiroplasma* representatives and from our samples. For our samples, sequences where obtained from PCR assays, while other sequences were retrieved from fragment of single-gene sequences or from assemblies available on GenBank (National Center for Biotechnology Information) (see manuscript for detailed information). The fragments of genes were retrieved from published genomes using a local `BLAST` and the following command lines:
```
makeblastdb -in ./Genomes-ref/$ref-genome.fasta -dbtype nucl -out $ref-genome_db
blastn -query Multiquery_Spiro.fasta -out ./output/SpiroMLST-genes_vs_$ref-genome.out -num_threads 6 -db $ref-genome_db -outfmt "6 qseqid sseqid sseq qlen pident nident mismatch evalue sstart send gapopen" -perc_identity 40

## If needed, genes showing low similarity percentage were retrieved using a tblastn:
tblastn -query Multiquery_Spiro_prot.fasta -out ./Other-genes_prot/SpiroMLST-prot-tblastn_vs_$ref-genome.out -db $ref-genome -num_threads 6 -outfmt "6 qseqid sseqid sseq qlen pident nident mismatch evalue sstart send gapopen" -evalue 1e-10
python getFasta.py $ref-genome.fasta query-contigsfragments.bed > contigsfragments.fasta
```

The phylogenetic tree was produced in a similar way that previously described:
```
for file in /CONC/*
do mafft "$file" > "$file"
done

cp ./CONC/*_align.fasta ./CONC_trimal/
for file in /CONC_trimal/*
do trimal -in "$file" -out "$file" -fasta -gt 1 -cons 50
done

AMAS.py concat -f fasta -d dna --in-files ./CONC_trimal/*.fasta

modeltest-ng -i Spiro-MLST.fasta -p 12 -T raxml -d nt

raxmlHPC-PTHREADS -T 8 -f a -s Spiro-MLST.fasta -n Spiro-MLST.boot -m GTRGAMMAIX -x 1234 -# 1000 -p 123
```
