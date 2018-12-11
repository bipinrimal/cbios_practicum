# qiime tools import \
# --input-path emp_cr_silva_16S_123.subset_10k.biom \
# --type 'FeatureTable[Frequency]' \
# --input-format BIOMV210Format \
# --output-path 10_biom.qza


# qiime feature-table filter-samples \
# --i-table 10_biom.qza --m-metadata-file \
# emp_qiime_mapping_subset_10k.tsv --p-where "host_taxid='9606'" \
# --o-filtered-table subset_10k.qza

# qiime tools extract \
# --input-path subset_10k.qza \
# --output-path feature-table.biom


########################################
###### Script to be run in bash ########
########################################

# Set the working directory
setwd("~/Classes/fall_2018/cbios_practicum")

########################################
###### Load the Libraries ##############
########################################

#source("https://bioconductor.org/biocLite.R")
library(phyloseq)
library(biomformat)
library(ape)
library(phytools)
library(ggplot2)
library(rhdf5)
source('read_hdf5.r')

biomfile <- 'feature-table.biom' #subsetted biom file
treefile <- '97_otus.tre'
mapfile <- 'metadata.tsv'
taxfile<-'taxonomy.txt' #obtained from biom file


# Loading mapping file and biom object.
map<-read.csv(mapfile,sep='\t',row.names = 1)
# Subset it to human samples
map_human<-subset(map,map$host_taxid=='9606')
map_human<-sample_data(map_human)
rm(map)

# Reading the subsetted biomfile
biom_object <- read_biom(biomfile)
otumat<-as.matrix(biom_data(biom_object))
OTU=otu_table(otumat,taxa_are_rows = TRUE)
rm(biom_object)

# the taxonomy file 
taxmat <- read.table("taxonomy.txt", sep=";",fill=TRUE)
colnames(taxmat)<-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
first_split<-strsplit(as.character(taxmat$Domain),"\t")
taxmat$Domain <-sapply(first_split,'[',2)
rownames(taxmat)<-sapply(first_split,'[',1)
taxmat<-as.matrix(taxmat)
TAX<-tax_table(taxmat)

physeq1<-merge_phyloseq(OTU,TAX,map_human)
# check file formats
class(otumat)
class(taxmat)

phytree<-read.tree(treefile)
physeq2 <-merge_phyloseq(physeq1,map_human,phytree)
