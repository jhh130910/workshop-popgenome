```R
# Loading module
library(PopGenome)
# Reading data, read fasta from folder
GENOME.class <- readData("fasta")
get.sum.data(GENOME.class)
get.individuals(GENOME.class)

# Question 1
# n.sites n.biallelic.sites n.gaps n.unknowns n.valid.sites n.polyallelic.sites trans.transv.ratio

# Question 2 a
GENOME.class <- readData("fasta_a")
get.sum.data(GENOME.class) # one sites less, important when data supposed to be in frame

# Question 2 b
GENOME.class <- readData("fasta_b")
get.sum.data(GENOME.class) # different number of individuals per sample
get.individuals(GENOME.class)

# Question 2 c
#GENOME.class <- readData("fasta_c")
#get.sum.data(GENOME.class) # incorrect alignment

# Available statistics and examples
GENOME.class <- readData("fasta")
show.slots(GENOME.class)
# Run necessary module
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class)
GENOME.class@n.sites
GENOME.class@Pi
GENOME.class@Tajima.D

# Question 3
# Data, FST, Neutrality, linkage, Recomb, Sweeps, Detail, Gff/Gtf

# Question 4
GENOME.class <- linkage.stats(GENOME.class)
GENOME.class@Wall.B

# Question 5
GENOME.class@Pi / GENOME.class@n.valid.sites

# Available region data and statistics
GENOME.class@region.data
GENOME.class@region.stats
# Examples
GENOME.class@region.data@biallelic.sites[[1]][1:10]
GENOME.class@region.data@transitions[[1]][1:10]

# Question 6, number of gapped sites
GENOME.class@region.data@sites.with.gaps
length(GENOME.class@region.data@sites.with.gaps[[1]]) # 1454

# Question 7, singletons
GENOME.class@region.data@n.singletons
sum(GENOME.class@region.data@n.singletons[[1]]) # 67

# Question 8

#11 The slot region.data
#During the reading process PopGenome will store some SNP specific information in the
#slot GENOME.class@region.data. This slot will for example store the genomic posi-
#tion of each SNP GENOME.class@region.data@biallelic.sites. In general, all in-
#formations here are stored as numeric vectors of length = n.biallelic.sites. Just typing
#GENOME.class@region.data will print a summary of the available slots. When multi-
#ple files have been read in the slots of the object of class region.data are organized
#as lists. Each element of the list is accessible via [[region.id]], where region.id
#is the identifier of the file of interest. The corresponding information is stored in the
#slot GENOME.class@region.names. In case of transformed GENOME objects e.g per-
#formed by sliding.window.transform [[region.id]] will be the identifier for the window
#of interest.
#12 The slot region.stats
#In some cases a multi-locus-scale representation of the statistic values is not possible and
#we were forced to organize those values as a list. In the slot GEOME.class@region.stats
#for example we can find the slot haplotype.counts which contains the haplotype distri-
#bution of each population. Here, the haplotypes regarding the whole population (whole
#data set) was specified (n.haplotypes=n.columns). Each row corresponds to one pop-
#ulation and the sum of each line is the sample size of each population. Obviously,
#the dimension of this matrix can differ between regions/windows. As described in the
#previous section specific files or regions/windows are accessible via [[region.id]].

# Without defining populations
get.individuals(GENOME.class)
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]]
# Define populations with lists
GENOME.class <- set.populations(GENOME.class,list(
c("CON","KAS-1","RUB-1","PER-1","RI-0","MR-0","TUL-0"),
c("MH-0","YO-0","ITA-0","CVI-0","COL-2","LA-0","NC-1")
))
# Check whether grouping is set correctly
GENOME.class@region.data@populations
GENOME.class@region.data@populations2
GENOME.class@region.data@outgroup
# Recalculate statistics for populations
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
GENOME.class@Tajima.D
# Each population
get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]]
# Set an outgroup
GENOME.class <- set.outgroup(GENOME.class,c("Alyrata"))
GENOME.class@region.data@outgroup
GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
get.neutrality(GENOME.class)[[1]]
get.neutrality(GENOME.class)[[2]

# Question 9
# Zeng.E, Fay.Wu.H

# Question 10
# The outgroup needs to be defined as a population, not as an outgroup



# What parameters need to be defined
library(PopGenome)
GENOME.class<- readVCF("great_tit/vcf/LGE22.vcf.gz", 
6000,"chrLGE22_Parus_Major_build_1.0.2",1,773534)
GENOME.class@region.names
GENOME.class<- neutrality.stats(GENOME.class, FAST=TRUE)
get.sum.data(GENOME.class)
GENOME.class@region.data

# Question 11
# filename, no SNPs per chunk, chr, start and stop

# Question 12
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class@theta_Watterson # 1161.065
GENOME.class@Pi # 1153.292
GENOME.class@Tajima.D # -2.529237

GENOME2.class <- readData("great_tit/vcf2",format="VCF",
gffpath="great_tit/gff")
get.sum.data(GENOME2.class)
GENOME2.class@region.data
#get.codons(GENOME.class,1)
#split <- splitting.data(GENOME2.class, subsites="CDS")
GENOME2.class <- set.synnonsyn(GENOME2.class, 
ref.chr="great_tit/fasta/LGE22.fasta")
GENOME2.class@region.data@synonymous
GENOME2.class@region.data@CodingSNPS

#Question 13

GENOME2.class.syn <- neutrality.stats(GENOME2.class,subsites="syn")
GENOME2.class.syn@Tajima.D
GENOME2.class.syn@theta_Watterson
#syn  sites
#D= -0.1095918
#theta= 44.2859

GENOME2.class.non <- neutrality.stats(GENOME2.class,subsites="nonsyn")
GENOME2.class.non@Tajima.D
GENOME2.class.non@theta_Watterson
#nonsyn  sites
#D=-0.4723517
#theta=13.82584




library(PopGenome)

# SPLIT VCF FILE
VCF_split_into_scaffolds("rad/variants.vcf","rad_split_vcf")
# READ IN DATA 
GENOME.class <- readData("rad_split_vcf_small",format="VCF")
pop1 <-as.character(read.table("rad/ind_species1.txt")[[1]])
pop2 <-as.character(read.table("rad/ind_species2.txt")[[1]])
GENOME.class <- set.populations(GENOME.class,list(pop1,pop2),diploid=TRUE)
# CHECK
GENOME.class@populations


# NEUTRALITY STATISTICS
GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE)
get.neutrality(GENOME.class)[[1]]
GENOME.class@Tajima.D
# FST
GENOME.class <- F_ST.stats(GENOME.class)
get.F_ST(GENOME.class)[[1]]
GENOME.class@nucleotide.F_ST
# DIVERSITY
GENOME.class <- diversity.stats(GENOME.class)
get.diversity(GENOME.class)
GENOME.class@nuc.diversity.within
# SFS
GENOME.class <- detail.stats(GENOME.class,site.spectrum=TRUE,site.FST=TRUE)
results <- get.detail(GENOME.class)
GENOME.class@region.stats@minor.allele.freqs

# PLOT
CON <- concatenate.regions(GENOME.class)
CON <- detail.stats(CON,site.spectrum=TRUE,site.FST=TRUE)
results <- get.detail(CON)
allele_Freqs <- CON@region.stats@minor.allele.freqs[[1]]
freq.table <- list()
freq.table[[1]] <- table(allele_Freqs)
sfs <- data.frame(freq.table)

library(ggplot2)
ggplot(sfs, aes(x=allele_Freqs, y=Freq)) + geom_bar(stat = 'identity')
```
