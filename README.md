# To do list
### Genome retrieving
GRM-P1 & GRM-P3 : binning concoct + anvio OK  
SWD-P1 & SWD-P2 : test d'extract reads via Kraken et bwa align avec CAS1 puis MEGAHIT puis binning  

### Genomic analyses
INCLUDING: all Sixodetis genomes (ticks & others), Smirum, Mycoplasma patho, OTHERS SPIRO CLADES?  
Virulence Mycoplasma:  
RIP:  
Ankyrin:   
MK genes:  
T4SS (Smirum already published and Mycoplasma):  
Genomic content organized by function for the 4 genomes:  

### Phylogenies
16S : v6 1-500; v6 1-700; v6 1-1400; v7 650-1400  
MLST : v6  


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
**NEED TO TEST: https://rhysnewell.github.io/rosella/usage/ & https://bitbucket.org/berkeleylab/metabat/src/master/**

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
miComplete miComplete_Spiro-human.tab --hmms Bact105 #miComplete_Spiro-human.tab a tabular separated file containing per line both a path to a genome and the type (i.e. fna)
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
Name	Length	GC-content	Present Markers	Completeness	Redundancy	ContigsN50	L50	N90	L90	CDs	
Spiro-CAS1-contigs	1292297	23.0	82	0.7810	1.1585	195	11017	35	3016	118	1971	
Spiro-CAS3-contigs	1359048	24.56	81	0.7714	1.1481	229	10069	45	2560	146	2025	
```

## 1.5. Annotation
```
prokka $ech-genome.fasta --locustag $ech --prefix $ech --outdir Prokka-$ech --rfam --compliant --cpus 6

bakta --db /share/banks/bakta_db_2401/db/ --verbose --compliant --output Bakta-SpiroGRMP1/ --prefix SpiroGRMP1 --locus-tag SPIROGRMP1 --genus Spiroplasma --species ixodetis --strain GRM-P1 --threads 4 Spiro-CAS1-contigs.fasta
```

## 1.6. Genomes' visualization
`GCview` 
```
perl ./cgview_xml_builder.pl -sequence ./genome-$ech.gbk -output Spiro-$ech.xml -gc_skew F -gc_content F -size large-v2 -gene_decoration arc -tick_density 0.02 -custom backboneThickness=20 featureThickness=200 featureSlotSpacing=60
java -jar ./cgview.jar -i Spiro-$ech.xml -o map_Spiro-$ech.png -f png
```

# Step 2: Genomes' description and comparison with others *Spiroplasma* genomes 
## 2.1. Phylogeny-based taxonomic assignation 
### 2.1.1. Based on *Spiroplasma* MAGs and whole genomes
First, single-copy orthologs (SCO) were identified using `OrthoFinder` (https://github.com/davidemms/OrthoFinder, OrthoFinder: phylogenetic orthology inference for comparative genomics, Emms DM and Kelly S, Genome Biology, 2019, doi: 10.1186/s13059-019-1832-y) from a set of specimens' genomes chosen to study the phylogenetic relationships between the obtained *Spiroplasma* MAGs and other *Spiroplasma* representatives:
```
orthofinder -f ./OrthoFinder_genomes/ -t 4 -S blast ## OrthoFinder_genomes being a directory including all .faa files of specimens of interest
```
For each SCO, sequences were individually aligned using `mafft` (https://github.com/GSLBiotech/mafft, MAFFT multiple sequence alignment software version 7: improvements in performance and usability, Katoh K and Standley DM, Molecular Biology and Evolution, 2013 doi: 10.1093/molbev/mst010):
```
for file in /Single_Copy_Orthologue_Sequences/*
do mafft "$file" > "$file"
done
```
For each SCO, ambigious hypervariable regions were removed using `trimAl` (https://github.com/inab/trimal, trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses, Capella-Gutiérrez S, Silla-Martínez JM, and Gabaldón T, Bioinformatics, 2009, doi: 10.1093/bioinformatics/btp348):
```
cp ./Single_Copy_Orthologue_Sequences/*_align.fasta ./Single_Copy_Orthologue_Sequences_trimal/
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do trimal -in "$file" -out "$file" -fasta -gt 1 -cons 50
done
```
Then, all SCO sequences were concatenated using `Amas` (https://github.com/marekborowiec/AMAS, AMAS: a fast tool for alignment manipulation and computing of summary statistics, Borowiec ML, PeerJ, 2016, doi: 10.7717/peerj.1660) in a single file:
```
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do awk '/^>/{print ">organism" ++i; next}{print}' < "$file" > "${file%_align.fasta}_rename.fasta"
done
cp ./Single_Copy_Orthologue_Sequences_trimal/*_rename.fasta ./Single_Copy_Orthologue_Sequences_AMAS/
AMAS.py concat -f fasta -d aa --in-files ./Single_Copy_Orthologue_Sequences_AMAS/*.fasta
```
Substitution models were evaluated using `modeltest-ng` (https://github.com/ddarriba/modeltest, ModelTest-NG: A new and scalable tool for the selection of DNA and protein evolutionary models, Darriba D, Posada D, Kozlov AM, Stamatakis A, Morel B, and Flouri T, Molecular Biology and Evolution, 2020, doi: 10.1093/molbev/msz189) in order to determinate the appropriate substitution model (according to AICc criterion) to use for the phylogenetic tree construction with `RAxML-NG` (https://github.com/stamatak/standard-RAxML, RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference, Kozlov A.M, Darriba D, Flouri T, Morel B, and Stamatakis A, Bioinformatics, 2019, doi: 10.1093/bioinformatics/btz305):
```
modeltest-ng -i Spiro-human_concatenated.faa -p 12 -T raxml -d aa
raxml-ng ...
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

raxml-ng --all --msa Spiro-MLST.fasta --model GTR+I+G4 --prefix Spiro-MLST --seed 5 --threads 8 --bs-trees 1000
raxml-ng --support --tree Spiro-MLST.raxml.bestTree --bs-trees 1000 --prefix SpiroMLST-boot --threads 2
```
## 2.2. fastANI and EZAAI
On the same dataset that used for phylogenomics:
```
fastANI --rl list_genomes_ANI_Spiro-humains.txt --ql list_genomes_ANI_Spiro-humains.txt -o fastani_Spiro_humains.txt --matrix

ezaai convert -i ./Genomes_Bakta_identity/SZISE.faa -s prot -o SZISE_db
ezaai calculate -i db_EZAAI_Spiro_humains/ -j db_EZAAI_Spiro_humains/ -o EZAAI_Spiro_humains_matrix.txt -t 6
```
