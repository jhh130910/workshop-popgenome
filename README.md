<http://tonig-evo.github.io/workshop-popgenome>

Introduction
============

PopGenome is an R package
(<https://cran.r-project.org/web/packages/PopGenome/index.html>) for
analyses of population genomic data \citep{Pfeifer2014}. For this
tutorial, please make sure that your R working directory is set
correctly and you have all the packages installed, e.g.
install.packages("PopGenome"). The following files are necessary to
conduct this practical session:

-   fasta\_file.txt, a fasta file for one locus from different
    *Arabidopsis thaliana* individuals (accessions) and the outgroup
    sequence from *Arabidopsis lyrata* in the folder **fasta**

-   LGE22.gff, LGE22.vcf, LGE22.fa in subfolders in the folder
    **great\_tit**

-   variants.vcf, ind\_species1.txt, ind\_species2.txt, rad\_assembly.fa
    in the folder **rad**

To get an overview about the file contents, inspect files with a text
Editor (e.g. Notepad+) or via the command line (e.g. more) before
starting to work on them.

Load library into R
===================

```R
\# Loading module 
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

1. What statistics can one obtain from get.sum.data function?

Folders **fasta\_a**, **fasta\_b**, **fasta\_c** contain modified
alignments. Identify the differences between the datasets. Why PopGenome
fails to load fasta files?

Obtaining summary statistics from alignments
--------------------------------------------

**Note:** Since calculation of certain population genetic parameter can
be computational intense, they have to be executed separately
beforehand. For this modules have to be run. Note that module **Fst**
has to be executed with **F\_st**. The statistic Tajima’s D is part of
the module **neutrality** not **Fst**

\# Available statistics and examples show.slots(GENOME.class) \# Run
necessary module GENOME.class \<- F~S~T.stats(GENOME.class) GENOME.class
\<- neutrality.stats(GENOME.class) GENOME.class@n.sites GENOME.class@Pi
GENOME.class@Tajima.D

What different modules are available? (show.slots)

What module is necessary to be executed in order to obtain Wall.B
\citep{Wall1999}?

How could one obtain a per site estimate of Pi?

Obtaining statistics for regions
--------------------------------

\# Available region data and statistics GENOME.class@region.data
GENOME.class@region.stats \# Examples
GENOME.class@region.data@biallelic.sites[[1]][1:10]
GENOME.class@region.data@transitions[[1]][1:10]

How many sites have gaps?

How many singletons are in the dataset?\
(see also **An\_introduction\_to\_the\_PopGenome\_package.pdf**, section
3.1)

What is the difference between *region.data* and *region.stats*?\
(see also **Whole\_genome\_analyses\_using\_VCF\_files.pdf**, section 11
and 12)

Define outgroups and populations
--------------------------------

**Note:** If one ore more outgroup sequences are defined, PopGenome will
only consider SNPs where the outgroup is monomorphic; the monomorphic
nucleotide is then automatically defined as the major allele (encoded by
0).

\# Without defining populations get.individuals(GENOME.class)
GENOME.class \<- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]] \# Define populations with lists
GENOME.class \<- set.populations(GENOME.class,list(
c("CON","KAS-1","RUB-1","PER-1","RI-0","MR-0","TUL-0"),
c("MH-0","YO-0","ITA-0","CVI-0","COL-2","LA-0","NC-1") )) \# Check
whether grouping is set correctly GENOME.class@region.data@populations
GENOME.class@region.data@populations2 GENOME.class@region.data@outgroup
\# Recalculate statistics for populations GENOME.class \<-
neutrality.stats(GENOME.class,detail=TRUE) GENOME.class@Tajima.D \# Each
population get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]] \# Set an outgroup GENOME.class \<-
set.outgroup(GENOME.class,c("Alyrata"))
GENOME.class@region.data@outgroup GENOME.class \<-
neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]] get.neutrality(GENOME.class)[[2]

Name implemented statistics that require an outgroup, e.g. that are
calculated after defining the outgroup.

What do you have to pay attention to when applying the McDonald-Kreitman
test?\
(see **Whole\_genome\_analyses\_using\_VCF\_files.pdf**)

Analysing VCF files for whole genome data
=========================================

Loading VCF files
-----------------

There are two ways to read in VCF files, either a folder of VCF files
with **readData** or a single VCF file with **readVCF**. To read a VCF
file using *readVCF* it needs to be compressed with *bgzip* and indexed
with *tabix*. The tabix files need to be placed in the same folder as
the vcf file.

\# What parameters need to be defined GENOME.class\<-
readVCF("great~t~it/vcf/LGE22.vcf.gz",
6000,"chrLGE22~P~arus~M~ajor~b~uild~1~.0.2",1,773534)
GENOME.class@region.names GENOME.class\<- neutrality.stats(GENOME.class,
FAST=TRUE) get.sum.data(GENOME.class) GENOME.class@region.data

What parameters need to be defined to *readVCF*? (see **PopGenome.pdf**)

What is the overall diversity (theta and pi) of chromosome LGE22?

Loading VCF files with annotation
---------------------------------

GENOME2.class \<- readData("great~t~it/vcf2",format="VCF",
gffpath="great~t~it/gff") get.sum.data(GENOME2.class)
GENOME2.class@region.data \#get.codons(GENOME.class,1) \#split \<-
splitting.data(GENOME2.class, subsites="CDS") GENOME2.class \<-
set.synnonsyn(GENOME2.class, ref.chr="great~t~it/fasta/LGE22.fasta")
GENOME2.class@region.data@synonymous
GENOME2.class@region.data@CodingSNPS GENOME2.class.syn \<-
neutrality.stats(GENOME2.class,subsites="syn")
GENOME2.class.syn@Tajima.D GENOME2.class.syn@theta~W~atterson

What is theta Watterson and Tajima’s D of chromosome LGE22 for
synonymous and nonsynonymous sites?

Analysing RADseq data using VCF 
================================

\# SPLIT VCF FILE
VCF~s~plit~i~nto~s~caffolds("rad/variants.vcf","rad~s~plit~v~cf") \#
READ IN DATA GENOME.class \<-
readData("rad~s~plit~v~cf~s~mall",format="VCF") pop1
\<-as.character(read.table("rad/ind~s~pecies1.txt")[[1]]) pop2
\<-as.character(read.table("rad/ind~s~pecies2.txt")[[1]]) GENOME.class
\<- set.populations(GENOME.class,list(pop1,pop2),diploid=TRUE) \# CHECK
GENOME.class@populations

Obtaining statistics from multiple VCFs derived from RADseq
-----------------------------------------------------------

\# NEUTRALITY STATISTICS GENOME.class \<- neutrality.stats(GENOME.class,
FAST=TRUE) get.neutrality(GENOME.class)[[1]] GENOME.class@Tajima.D \#
FST GENOME.class \<- F~S~T.stats(GENOME.class)
get.F~S~T(GENOME.class)[[1]] GENOME.class@nucleotide.F~S~T \# DIVERSITY
GENOME.class \<- diversity.stats(GENOME.class)
get.diversity(GENOME.class) GENOME.class@nuc.diversity.within \# SFS
GENOME.class \<-
detail.stats(GENOME.class,site.spectrum=TRUE,site.FST=TRUE) results \<-
get.detail(GENOME.class) GENOME.class@region.stats@minor.allele.freqs

Plot a site-frequency-spectrum for each population.

CON \<- concatenate.regions(GENOME.class) CON \<-
detail.stats(CON,site.spectrum=TRUE,site.FST=TRUE) results \<-
get.detail(CON) allele~F~reqs \<-
CON@region.stats@minor.allele.freqs[[1]] freq.table \<- list()
freq.table[[1]] \<- table(allele~F~reqs) sfs \<- data.frame(freq.table)

library(ggplot2) ggplot(sfs, aes(x=allele~F~reqs, y=Freq)) +
geom~b~ar(stat = ’identity’)

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
further details see **Whole\_genome\_analyses\_using\_VCF\_files.pdf**).
