>This repository details command lines used in the production of the "**A unique evolutionary lineage of *Spiroplasma ixodetis* associated with infantile cataracts and adult febrile illness in Europe**" manuscript by Buysse et al. 

[![DOI](https://zenodo.org/badge/626857615.svg)](https://doi.org/10.5281/zenodo.15862248)

### Table of contents
- [1. Retrieving *Spiroplasma* genomes from biological samples](#1-retrieving-spiroplasma-genomes-from-biological-samples)
  - [1.1. Adapters' trimming](#11-adapters-trimming)
  - [1.2. Assembly](#12-assembly)
  - [1.3. Binning and retrieving *Spiroplasma* MAGs](#13-binning-and-retrieving-spiroplasma-mags)
  - [1.4. Quality check of the final *Spiroplasma* MAGs](#14-quality-check-of-the-final-spiroplasma-mags)
  - [1.5. Annotation](#15-annotation)
  - [1.6. MAGs' visualization](#16-mags-visualization)
- [2. Genomes' description and comparison with others *Spiroplasma* genomes](#2-genomes-description-and-comparison-with-others-spiroplasma-genomes)
  - [2.1. Identity matrix (ANI and AAI)](#21-identity-matrix-ani-and-aai)
  - [2.2. Genomic content comparison between MAGs, *S. ixodetis* and *S. platyhelix* genomes](#22-genomic-content-comparison-between-mags-s-ixodetis-and-s-platyhelix-genomes)
  - [2.3. COG categories prediction](#23-cog-categories-prediction)
  - [2.4. Phylogenomics](#24-phylogenomics)
  - [2.5. Detection of protein sequences](#25-detection-of-protein-sequences)
- [3. Gene-based phylogenetic analyses](#3-gene-based-phylogenetic-analyses)
  - [3.1. Based on single gene (16S rDNA gene)](#31-based-on-single-gene-16S-rdna-gene)
  - [3.2. Based on multi-locus sequencing typing (MLST)](#32-based-on-multi-locus-sequencing-typing-mlst)


# 1. Retrieving *Spiroplasma* genomes from biological samples
In this study, we aimed to retrieve and analyze the *Spiroplasma* MAGs (Metagenome-Assembled Genome) from biological samples collected from two symptomatic infants with ocular infections (i.e. lensectomy samples). These cases have been reported in [Farassat N et al. 2021](https://bmcophthalmol.biomedcentral.com/articles/10.1186/s12886-021-02201-0).

For each sample, a set of paired-end reads is available in NCBI (SRA: XXXX & XXXX). Details about the experimental and sequencing methods are available in the related manuscript.

Following analyses are based on these reads' sets, and those are referred to as:
```
ech="GRMP1 GRMP3"
```

## 1.1. Adapters' trimming
Adapters' sequences were removed from raw reads using `Cutadapt` (see details [here](https://cutadapt.readthedocs.io/en/stable/)):
```
cutadapt -b file:$ech-adaptersF -o $ech-R1-trimmed.fastq.gz $ech-R1.fastq.gz
```

## 1.2. Assembly
Reads were assembled using `MEGAHIT` (see [here](https://github.com/voutcn/megahit)):
```
megahit -1 $ech-R1-trimmed.fastq.gz -2 $ech-R2-trimmed.fastq.gz -o $ech-metaMEGAHIT 
```

## 1.3. Binning and retrieving *Spiroplasma* MAGs
*Spiroplasma* MAGs were retrieved from assemblies using `CONCOCT` (see details [here](https://github.com/BinPro/CONCOCT)) and the `anvi'o` pipeline (see details [here](https://anvio.org/)).

First, the contigs were formated to match the requirements of `anvi'o`:
```
sed '/^>/s/ .*//' $ech-contigs.fa > $ech-contigs-rename.fa
rm $ech-contigs.fa
mv $ech-contigs-rename.fa $ech-contigs.fa
```
The contigs were binned using `CONCOCT`: 
```
bwa index $ech-contigs.fa
bwa mem -t 4 $ech-contigs.fa $ech-R1.fastq.gz $ech-R2.fastq.gz | samtools sort -@ 4 -T mapped -O BAM -o $ech-reads-mapped.bam
samtools index $ech-reads-mapped.bam

cut_up_fasta.py $ech-contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed $ech-reads-mapped.bam > coverage_table.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b $ech-concoctMEGAHIT_output/ -t 4
merge_cutup_clustering.py $ech-concoctMEGAHIT_output/clustering_gt1000.csv > $ech-concoctMEGAHIT_output/clustering_merged.csv
mkdir $ech-concoctMEGAHIT_output/fasta_bins
extract_fasta_bins.py $ech-contigs.fa $ech-concoctMEGAHIT_output/clustering_merged.csv --output_path $ech-concoctMEGAHIT_output/fasta_bins
```
The `clustering_merged` CSV file produced by `CONCOCT` needed to be exported in a tabular-delimited TEXT file and the bins' names had to be renamed to match the requirements of `anvi'o` : 
```
awk '$2="bin"$2 {print}' bins-to-format.txt > bins-renamed.txt
```
The `bins-renamed.txt` file had to be transformed again to correspond to a tabular-delimited TEXT file, called `bins.txt` hereafter. The bins were taxonomically assigned using the `anvi'o` pipeline below: 
```
#To create the contigs database
anvi-gen-contigs-database -f $ech-contigs.fa -o $ech-metaMEGAHIT.db --ignore-internal-stop-codons -n Binning -T 4
anvi-run-hmms -c $ech-metaMEGAHIT.db

#To create the profile database
anvi-profile -i $ech-reads-mapped.bam -c $ech-metaMEGAHIT.db --min-contig-length 250 -T 4 -o $ech-PROFILE --cluster-contigs

#To import the bins as a collection
anvi-import-collection bins.txt -c $ech-metaMEGAHIT.db -p $ech-PROFILE/PROFILE.db -C bins --contigs-mode

#To assign the bins and visualize the results
anvi-run-scg-taxonomy -c $ech-metaMEGAHIT.db -T 2
anvi-estimate-scg-taxonomy -c $ech-metaMEGAHIT.db --output-file $ech-TAXONOMY.txt -p $ech-PROFILE/PROFILE.db -C bins --compute-scg-coverages -T 2
anvi-summarize -p $ech-PROFILE/PROFILE.db -c $ech-metaMEGAHIT.db -C bins -o $ech-SUMMARY 
```
Let's see the `$ech-SUMMARY` file to check the taxnonomy results and some stats about each bin.

## 1.4. Quality check of the final *Spiroplasma* MAGs
Quality and multiple statistics were accessed using `miComplete` (see details [here](https://pypi.org/project/micomplete/)) and `Quast` (see details [here](https://github.com/ablab/quast)):
```
#Quast
quast.py ./Spiro-GRM-genomes/*.fasta -o ./RESULTS-QUAST

#miComplete
miComplete miComplete_Spiro-human.tab --hmms Bact105 #miComplete_Spiro-human.tab a tabular separated file containing per line both a path to a genome and the type (i.e. fna)
```

Results are:
```
#Quast
Assembly                    Spiro-GRMP1         Spiro-GRMP3
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

#miComplete
Name	Length	GC-content	Present Markers	Completeness	Redundancy	Contigs	N50	L50	N90	L90	CDs	
Spiro-GRMP1	1292297	0.23	82	0.7810	1.1585	195	11017	35	3016	118	1971	
Spiro-GRMP3	1359048	0.25	81	0.7714	1.1481	229	10069	45	2560	146	2025
```

## 1.5. Annotation
Annotation was obtained using `Bakta` (see details [here](https://github.com/oschwengers/bakta)):
```
bakta --db /share/banks/bakta_db_2401/db/ --verbose --compliant --output Bakta-Spiro-$ech/ --prefix Spiro-$ech --locus-tag SPIRO-$ech --genus Spiroplasma --species ixodetis --strain -$ech --threads 4 Spiro--$ech.fasta --translation-table 4
```

## 1.6. MAGs' visualization
MAG representation was performed using `GCview` (see details [here](https://github.com/paulstothard/cgview)):
```
perl ./cgview_xml_builder.pl -sequence ./genome-$ech.gbk -output Spiro-$ech.xml -gc_skew F -gc_content F -size large-v2 -gene_decoration arc -tick_density 0.02 -custom backboneThickness=20 featureThickness=200 featureSlotSpacing=60
java -jar ./cgview.jar -i Spiro-$ech.xml -o map_Spiro-$ech.png -f png
```


# 2. Genomes' description and comparison with others *Spiroplasma* genomes 
## 2.1. Identity matrix (ANI and AAI)
Identity matrix were calcultated using `FastANI` (see details [here](https://github.com/ParBLiSS/FastANI)) and `EzAAI` (see details [here](https://github.com/endixk/ezaai)):
```
#fastANI
fastANI --rl list_genomes_ANI.txt --ql list_genomes_ANI.txt -o fastani_Spiro.txt --matrix

#EZAAI
ezaai convert -i $genome.faa -s prot -o $genome_db
ezaai calculate -i db_EZAAI_all_Spiro/ -j db_EZAAI_all_Spiro/ -o EZAAI_Spiro_matrix.txt -t 6
aai_data <- read.table("EZAAI_Spiro_matrix.txt", header=T,fill=T,sep="\t",dec=",")
colnames(aai_data) <- c("Label_1", "Label_2", "AAI")
aai_matrix <- aai_data %>% pivot_wider(names_from = Label_2, values_from = AAI)
aai_mat <- as.matrix(sapply(aai_matrix, as.numeric))
heatmaply(aai_mat)
```

## 2.2. Genomic content comparison between MAGs, *S. ixodetis* and *S. platyhelix* genomes
First, single-copy orthologs (SCO) were identified using `OrthoFinder` (see details [here](https://github.com/davidemms/OrthoFinder)), then differences were visualized using R environment and `ggVennDiagram` R package (see details [here](https://www.rdocumentation.org/packages/ggVennDiagram/versions/1.2.2)): 
```
#OrthoFinder
orthofinder -f ./OrthoFinder_genomes/ -t 4 -S blast ## OrthoFinder_genomes being a directory including all .faa files of specimens of interest

#R environment
ortho_tab <- read.table("GeneCount_GRMP_Sixo.txt", sep="\t", header=TRUE, fill=TRUE) 
SiRHIGRMP1 <- subset(ortho_tab, ortho_tab$GRMP1>0)
SiRHIGRMP1_list <- SiRHIGRMP1$Orthogroup
SiRHIGRMP3 <- subset(ortho_tab, ortho_tab$GRMP3>0)
SiRHIGRMP3_list <- SiRHIGRMP3$Orthogroup
Sixo <- subset(ortho_tab, ortho_tab$Sixo>0)
Sixo_list <- Sixo$Orthogroup
Spla <- subset(ortho_tab, ortho_tab$Spla>0)
Spla_list <- Spla$Orthogroup
Orthologs <- list(SiRHIGRMP1 = SiRHIGRMP1_list, SiRHIGRMP3 = SiRHIGRMP3_list, Sixo = Sixo_list, Spla = Spla_list)
p <- ggVennDiagram(Orthologs, label_percent_digit = 1, label_size = 4) 
plot_venn <- p + scale_fill_distiller(palette = "PuBu", direction = 1)
```

## 2.3. COG categories prediction
COG categories were predicted using `eggNOG-mapper` (see details [here](https://github.com/eggnogdb/eggnog-mapper)): 
```
emapper.py -i Spiro-$ech.faa --itype proteins --output $ech --trans_table 4 -m diamond
```

## 2.4. Phylogenomics
First, single-copy orthologs (SCO) were identified using `OrthoFinder` (see details [here](https://github.com/davidemms/OrthoFinder)) from a set of genomes chosen to study the phylogenetic relationships between the obtained *Spiroplasma* MAGs and other *Spiroplasma* representatives:
```
orthofinder -f ./OrthoFinder_genomes/ -t 4 -S blast ## OrthoFinder_genomes being a directory including all .faa files of specimens of interest
```
For each SCO, sequences were individually aligned using `mafft` (see details [here](https://github.com/GSLBiotech/mafft)):
```
for file in /Single_Copy_Orthologue_Sequences/*
do mafft "$file" > "$file"
done
```
For each SCO, ambigious hypervariable regions were removed using `trimAl` (see details [here](https://github.com/inab/trimal)):
```
cp ./Single_Copy_Orthologue_Sequences/*_align.fasta ./Single_Copy_Orthologue_Sequences_trimal/
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do trimal -in "$file" -out "$file" -fasta -gt 1 -cons 50
done
```
Then, all SCO sequences were concatenated using `Amas` (see details [here](https://github.com/marekborowiec/AMAS)):
```
for file in /Single_Copy_Orthologue_Sequences_trimal/*
do awk '/^>/{print ">organism" ++i; next}{print}' < "$file" > "${file%_align.fasta}_rename.fasta"
done
cp ./Single_Copy_Orthologue_Sequences_trimal/*_rename.fasta ./Single_Copy_Orthologue_Sequences_AMAS/
AMAS.py concat -f fasta -d aa --in-files ./Single_Copy_Orthologue_Sequences_AMAS/*.fasta
```
Substitution models were evaluated using `modeltest-ng` (see [here](https://github.com/ddarriba/modeltest)) in order to determinate the appropriate substitution model (according to AICc criterion) to use for the phylogenetic tree construction with `RAxML-NG` (see [here](https://github.com/amkozlov/raxml-ng)):
```
modeltest-ng -i Spiro-human_concatenated.faa -p 12 -T raxml -d aa
raxml-ng --all --msa Spiro-human_concatenated.faa --model LG+I+G4 --prefix Spirohuman-tree --seed 5 --threads 8 --bs-trees 1000
raxml-ng --support --tree Spirohuman-tree.raxml.bestTree --bs-trees 1000 --prefix Spirohuman-tree-boot --threads 2

```
The phylogenetic tree was visualized and modified using `figtree` (see details [here](https://github.com/rambaut/figtree)).


## 2.5. Detection of protein sequences
### 2.5.1. Detection of protein sequences with an OTU domain
To detect and summarize protein sequences containing an OTU domain across a set of genomes, we used HMM profile for PF02338 (id. for the OTU-like cysteine protease, <https://www.ebi.ac.uk/interpro/entry/pfam/PF02338/>). For this purpose, we ran a dedicated script that relies on `HMMER` (hmmsearch, see details [here](<https://github.com/EddyRivasLab/hmmer/tree/master>)) to search for the domain in a combined FASTA `.faa` file containing proteins from the target genomes.
To run the script, the files required are (i) a set of `Genomes.faa` (here SpiroGRMP1.faa, SpiroGRMP3.faa and other *Spiroplasma* genomes), and (ii) `PF02338.hmm` (to download [here](<https://www.ebi.ac.uk/interpro/entry/pfam/PF02338/>)).

```
# Concatenate all .faa files into one
cat *.faa > All_proteins.faa

# Run hmmsearch to get OTU domain table output
hmmsearch --tblout OTUmatches_All_proteins_resultats.tbl PF02338.hmm All_proteins.faa
hmmsearch --domtblout Work_OTUmatches_All_proteins_resultats.domtbl PF02338.hmm All_proteins.faa

# Extract protein IDs and their sequence lengths
awk '/^>/ {
  if (seq) print id "\t" length(seq);
  id = substr($1, 2); gsub(/ .*/, "", id); seq = ""; next
} { seq = seq $0 } END { if (seq) print id "\t" length(seq) }' All_proteins.faa > Work_protein_lengths.tsv

# Extract hits from the hmmsearch OTU domain table and calculate match lengths
awk '!/^#/ {
  prot = $1; evalue = $7;
  start = $20+0; end = $21+0;
  len = (end > start) ? end - start + 1 : start - end + 1;
  print prot "\t" evalue "\t" start "\t" end "\t" len
}' "Work_OTUmatches_All_proteins_resultats.domtbl" > Work_all_hits.tsv

# Sort protein lengths and hits by protein ID
sort -k1,1 Work_protein_lengths.tsv > Work_protein_lengths.sorted.tsv
sort -k1,1 Work_all_hits.tsv > Work_all_hits.sorted.tsv

# Join hits and protein lengths on protein ID
join -t $'\t' -1 1 -2 1 Work_all_hits.sorted.tsv Work_protein_lengths.sorted.tsv > final_results.tmp
# Add header and save final results
(echo -e "ProteinID\tEvalue\tStart\tEnd\tMatchLength\tProteinLength"; cat final_results.tmp) > "OTUmatches_final_results.tsv"
rm final_results.tmp
```

### 2.5.2. Phylogeny of the OTU domain
Homologs of the OTU-containing sequences from SiRHI were identified using `BLASTP` on the NCBI BLAST webserver (<https://blast.ncbi.nlm.nih.gov/Blast.cgi>, <https://doi.org/10.1016/S0022-2836(05)80360-2>), querying only the OTU domain sequence. For the OTU domain alignment, we also included sequences from Cif and Spaid. The sequences were aligned using `MAFFT` in `UGENE` (see details [here](<https://ugene.net/>)), and positions containing gaps (‘-‘) were subsequently removed. All the OTU sequences for the alignment are provided in the `OTU_ALIGNMENT.faa`.

Then, substitution models were evaluated using `modeltest` to determine the most appropriate ML substitution model (based on the AICc criterion), followed by phylogenetic tree construction with `raxml-ng`:

```
#Identification of the best-fitted ML substitution model
modeltest-ng -i OTU_ALIGNMENT.faa -p 12 -T raxml -d aa 
# The model: LG+G4

#Phylogeny and bootstrap estimation
raxml-ng --all --msa OTU_ALIGNMENT.faa --model LG+G4--prefix OTUtree-raxmlng --seed 5 --threads 4 --bs-trees 1000
raxml-ng --support --tree OTUtree -raxmlng.raxml.bestTree --bs-trees 1000 --prefix OTUtree -boot --threads 2
```

Finally, the phylogenetic tree was visualized and adapted using `figtree` and `MEGA11` (see details [here](<https://megasoftware.net/>))

### 2.5.3. Detection of protein sequences with other domains
Once the OTU-containing sequences were detected, we predicted additional protein domains in these sequences using the 
`HHpred` online webserver (see details [here](<https://toolkit.tuebingen.mpg.de/tools/hhpred>)) with default parameters, querying the `SCOPe70 (v2.08)`, `Pfam-A (v37.0)`, `SMART (v6.0)`, and `COG/KOG (v1.0)` databases.
For domain identification in Cif protein sequences, see details [here](<https://github.com/julien01A/cif_evolution/>).


# 3. Gene-based phylogenetic analyses
## 3.1. Based on single gene (16S rDNA gene)
Single-gene phylogenies were produced using the 16S rDNA gene sequences. This gene was used to describe *Spiroplasma* strains associated with previously reported cases. As different fragment lengths of the 16S rDNA gene were sequenced in those studies, several phylogenies will be produced to include the diversity of *Spiroplasma* strains. Sequences from other *Spiroplasma* representatives were retrieved from GenBank (National Center for Biotechnology Information), while sequences from published genomes were obtained using a local `BLAST` (see details [here](https://www.ncbi.nlm.nih.gov/books/NBK279690/)). 

Then, the phylogenetic trees were produced in a similar way that previously described (with adaptations to nucleotidic sequences):
```
modeltest-ng -i Spiro-16S-regions-regions.fasta -p 12 -T raxml -d nt
```

## 3.2. Based on multi-locus sequencing typing (MLST)
This phylogeny is based on an alignement of concatenated sequences of three genes (*dnaK*, *gyrA*, and *rpoB* genes) from *Spiroplasma* representatives and from our samples. For our samples, sequences where obtained from PCR assays, while other sequences were retrieved from GenBank and published genomes using a local `BLAST`. Then, the phylogenetic tree was produced in a similar way that previously described. 
