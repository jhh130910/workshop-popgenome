Resources and Preparation
============

Resources:

- Online: <http://tonig-evo.github.io/workshop-popgenome>
- Data files: [PopGenome_data.zip](https://github.com/tonig-evo/workshop-popgenome/blob/master/PopGenome_data.zip)
- Solutions: [solutions.md](https://github.com/tonig-evo/workshop-popgenome/blob/master/solutions.md)

PDFs:

- <https://github.com/tonig-evo/workshop-popgenome/blob/master/An_introduction_to_the_PopGenome_package.pdf>
- <https://github.com/tonig-evo/workshop-popgenome/blob/master/PopGenome.pdf>
- <https://github.com/tonig-evo/workshop-popgenome/blob/master/Whole_genome_analyses_using_VCF_files.pdf>

To work on iceberg, copy necessary files from fastdata
```
mkdir PopGenome
cd PopGenome
cp /usr/local/extras/Genomics/workshops/March2016/PopGenome/PopGenome_data.zip ./
unzip PopGenome_data.zip
rm PopGenome_data.zip
```

Introduction
============

PopGenome is an R package
(<https://cran.r-project.org/web/packages/PopGenome/index.html>) for
analyses of population genomic data. For this
tutorial, please make sure that your R working directory is set
correctly and you have all the packages installed, e.g.
install.packages("PopGenome"). The following files are necessary to
conduct this practical session:

-   fasta_file.txt, a fasta file for one locus from different
    *Arabidopsis thaliana* individuals (accessions) and the outgroup
    sequence from *Arabidopsis lyrata* in the folder **fasta**

-   LGE22.gff, LGE22.vcf, LGE22.fa in subfolders in the folder
    **great_tit**

-   variants.vcf, ind_species1.txt, ind_species2.txt, rad_assembly.fa
    in the folder **rad**

To get an overview about the file contents, inspect files with a text
Editor (e.g. Notepad+) or via the command line (e.g. more) before
starting to work on them.

Load library into R
===================

```R
# Loading module 
library(PopGenome)
```

Fasta Files
===========

Reading Fasta Files
-------------------

```R
# Reading data, read fasta from folder 
GENOME.class <- readData("fasta") 
get.sum.data(GENOME.class)
get.individuals(GENOME.class)
```

1.What statistics can one obtain from get.sum.data function?

2.Folders **fasta_a**, **fasta_b**, **fasta_c** contain modified alignments. Identify the differences between the datasets. Why PopGenome fails to load fasta files?

Obtaining summary statistics from alignments
--------------------------------------------

**Note:** Since calculation of certain population genetic parameter can
be computational intense, they have to be executed separately
beforehand. For this modules have to be run. Note that module **Fst**
has to be executed with **F_st**. The statistic Tajima’s D is part of
the module **neutrality** not **Fst**

```R
# Available statistics and examples 
show.slots(GENOME.class) 
# Run necessary module 
GENOME.class <- F_ST.stats(GENOME.class) 
GENOME.class <- neutrality.stats(GENOME.class) 
GENOME.class@n.sites 
GENOME.class@Pi
GENOME.class@Tajima.D
```

3.What different modules are available? (show.slots)

4.What module is necessary to be executed in order to obtain Wall.B?

5.How could one obtain a per site estimate of Pi?

Obtaining statistics for regions
--------------------------------

```R
# Available region data and statistics 
GENOME.class@region.data
GENOME.class@region.stats 
# Examples
GENOME.class@region.data@biallelic.sites[[1]][1:10]
GENOME.class@region.data@transitions[[1]][1:10]
```

6.How many sites have gaps?

7.How many singletons are in the dataset? (see also **An_introduction_to_the_PopGenome_package.pdf**, section
3.1)

8.What is the difference between *region.data* and *region.stats*? (see also **Whole_genome_analyses_using_VCF_files.pdf**, section 11 and 12)

Define outgroups and populations
--------------------------------

**Note:** If one ore more outgroup sequences are defined, PopGenome will
only consider SNPs where the outgroup is monomorphic; the monomorphic
nucleotide is then automatically defined as the major allele (encoded by
0).

```R
# Without defining populations 
get.individuals(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]] 
# Define populations with lists
GENOME.class <- set.populations(GENOME.class,list(
c("CON","KAS-1","RUB-1","PER-1","RI-0","MR-0","TUL-0"),
c("MH-0","YO-0","ITA-0","CVI-0","COL-2","LA-0","NC-1") )) 
# Check whether grouping is set correctly 
GENOME.class@region.data@populations
GENOME.class@region.data@populations2 
GENOME.class@region.data@outgroup
# Recalculate statistics for populations 
GENOME.class <-neutrality.stats(GENOME.class,detail=TRUE) 
GENOME.class@Tajima.D 
# Each population 
get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]] 
# Set an outgroup 
GENOME.class <-set.outgroup(GENOME.class,c("Alyrata"))
GENOME.class@region.data@outgroup 
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]] 
get.neutrality(GENOME.class)[[2]]
```

9.Name implemented statistics that require an outgroup, e.g. that are
calculated after defining the outgroup.

10.What do you have to pay attention to when applying the McDonald-Kreitman
test?
(see **Whole_genome_analyses_using_VCF_files.pdf**)

Analysing VCF files for whole genome data
=========================================

Loading VCF files
-----------------

There are two ways to read in VCF files, either a folder of VCF files
with **readData** or a single VCF file with **readVCF**. To read a VCF
file using *readVCF* it needs to be compressed with *bgzip* and indexed
with *tabix*. The tabix files need to be placed in the same folder as
the vcf file.
```R
# What parameters need to be defined 
GENOME.class <-readVCF("great_tit/vcf/LGE22.vcf.gz", 6000,"chrLGE22_Parus_Major_build_1.0.2",1,773534)
GENOME.class@region.names 
GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE) 
get.sum.data(GENOME.class) 
GENOME.class@region.data
```
11.What parameters need to be defined to *readVCF*? (see **PopGenome.pdf**)

12.What is the overall diversity (theta and pi) of chromosome LGE22?

Loading VCF files with annotation
---------------------------------
```R
GENOME2.class <- readData("great_tit/vcf2",format="VCF", gffpath="great_tit/gff") 
get.sum.data(GENOME2.class)
GENOME2.class@region.data 
GENOME2.class <- set.synnonsyn(GENOME2.class, ref.chr="great_tit/fasta/LGE22.fasta")
GENOME2.class@region.data@synonymous
GENOME2.class@region.data@CodingSNPS 
GENOME2.class.syn <- neutrality.stats(GENOME2.class,subsites="syn")
GENOME2.class.syn@Tajima.D 
GENOME2.class.syn@theta_Watterson
```
13.What is theta Watterson and Tajima’s D of chromosome LGE22 for synonymous and nonsynonymous sites?

Analysing RADseq data using VCF 
================================

```R
# SPLIT VCF FILE
VCF_split_into_scaffolds("rad/variants.vcf","rad_split_vcf") 
# READ IN DATA, smaller subset
GENOME.class <- readData("rad_split_vcf_small",format="VCF") 
pop1<-as.character(read.table("rad/ind_species1.txt")[[1]]) 
pop2<-as.character(read.table("rad/ind_species2.txt")[[1]]) 
GENOME.class<- set.populations(GENOME.class,list(pop1,pop2),diploid=TRUE) 
# CHECK
GENOME.class@populations
```
Obtaining statistics from multiple VCFs derived from RADseq
-----------------------------------------------------------
```R
# NEUTRALITY STATISTICS 
GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE) 
get.neutrality(GENOME.class)[[1]] GENOME.class@Tajima.D 
# FST 
GENOME.class <- F_ST.stats(GENOME.class)
get.F_ST(GENOME.class)[[1]] 
GENOME.class@nucleotide.F_ST 
# DIVERSITY
GENOME.class <- diversity.stats(GENOME.class)
get.diversity(GENOME.class) 
GENOME.class@nuc.diversity.within 
# SFS
GENOME.class <-detail.stats(GENOME.class,site.spectrum=TRUE,site.FST=TRUE) 
results <- get.detail(GENOME.class) 
GENOME.class@region.stats@minor.allele.freqs
```

14.Plot a site-frequency-spectrum for each population.

```R
# Concatenate loci
CON <- concatenate.regions(GENOME.class) 
CON <- detail.stats(CON,site.spectrum=TRUE,site.FST=TRUE) 
results <-get.detail(CON) 
allele_Freqs <- CON@region.stats@minor.allele.freqs[[1]] 
freq.table <- list()
freq.table[[1]] <- table(allele_Freqs) 
sfs <- data.frame(freq.table)

library(ggplot2) 
ggplot(sfs, aes(x=allele_Freqs, y=Freq)) + geom_bar(stat = ’identity’)
```
Additional aspects
==================

More information are available in three pdfs accompanied by the program
(see folder **pdf**): **An introduction to the PopGenome package**:
Sliding window analysis, reading SNP data files, coalescent simulations;
**Whole genome analyses using PopGenome and VCF files**: Details about
reading tabixed VCF files, examples, graphical output, parallel read-in,
pre-filtering VCF files; **Package PopGenome**: Documentation about all
implemented functions with examples

Including features from GFF files to Fasta files
------------------------------------------------

If no gff-file was specified when the data was read in, it is assumed
that the alignment is in the correct reading frame (starting at a first
codon position). The GFF folder contains GFF-files for each alignment
stored in the FASTA folder. The GFF files should have the same names
(without any extensions like .fas or .gff) as the corresponding FASTA
files to ensure that sequence and annotation are matched correctly.

Handling missing data and differences between readData and readVCF
------------------------------------------------------------------

PopGenome can use missing data, e.g. positions with gaps. Also note the
differences between *readData* and *readVCF* for your own analysis (for
further details see **Whole_genome_analyses_using_VCF_files.pdf**).
