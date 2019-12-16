---
title: "CBIOS_Practicum"
author: "Sofia Roitman"
date: "9/24/2018"
output: html_document
---

## Set working directory
```{r message=FALSE, error=FALSE}
setwd("/Users/sofiaroitman/Google_Drive/GRAD_SCHOOL/Fall_2018/CBIOS_practicum/")
```

#R housekeeping
Load libraries and functions
```{r message=FALSE, error=FALSE}
# libraries and functions
library(tidyverse)
library(vegan)
library(xlsx)
library(BiocInstaller)
library(phyloseq)
library(biomformat)
library(ape)
library(Biostrings)
library(genefilter)
library(magrittr)
library(DESeq2)
library(plyr)
library(data.table)
library(phytools)
library(stringi)
library(ggplot2)
library(rhdf5)
library(microbiome)
library("RVAideMemoire")
source('read_hdf5.r')
```

## Import and load data
```{r}
# Download subsetted phyloseq file
rt_raw<-readRDS('subsetted_physeq2.rds')
rt_highcover_ns_abun <- rt_raw


```

## Subset samples for later use
```{r}
# Subset: just village
rt_highcover_ns_abun_vil = subset_samples(rt_highcover_ns_abun, env_biome == "village biome")

# Subset: just urban
rt_highcover_ns_abun_urb = subset_samples(rt_highcover_ns_abun, env_biome == "urban biome")

# Subset: just sebum
rt_highcover_ns_abun_sebum = subset_samples(rt_highcover_ns_abun, env_material == "sebum")

# Subset: just mucus
rt_highcover_ns_abun_mucus = subset_samples(rt_highcover_ns_abun, env_material == "mucus")

# Subset: just saliva
rt_highcover_ns_abun_saliva = subset_samples(rt_highcover_ns_abun, env_material == "saliva")

# Subset: just feces
rt_highcover_ns_abun_feces = subset_samples(rt_highcover_ns_abun, env_material == "feces")

```


## Create distance matrices
```{r message=FALSE, error=FALSE}
# All samples
bc_rt_highcover_ns_abun <- phyloseq::distance(otu_table(rt_highcover_ns_abun), "bray")

# Just village samples
bc_rt_highcover_ns_abun_vil <- phyloseq::distance(otu_table(rt_highcover_ns_abun_vil), "bray")

# Just urban samples
bc_rt_highcover_ns_abun_urb <- phyloseq::distance(otu_table(rt_highcover_ns_abun_urb), "bray")

# Just sebum samples
bc_rt_highcover_ns_abun_seb <- phyloseq::distance(otu_table(rt_highcover_ns_abun_sebum), "bray")

# Just mucus samples
bc_rt_highcover_ns_abun_muc <- phyloseq::distance(otu_table(rt_highcover_ns_abun_mucus), "bray")

# Just saliva samples
bc_rt_highcover_ns_abun_sal <- phyloseq::distance(otu_table(rt_highcover_ns_abun_saliva), "bray")

# Just feces samples
bc_rt_highcover_ns_abun_fec <- phyloseq::distance(otu_table(rt_highcover_ns_abun_feces), "bray")

```

## Use ggplot theme to make everything more legible
```{r}
ggtheme <- theme(axis.title = element_text(colour="black",family = "Helvetica",
                                           size = rel(1.5)), 
                 axis.text = element_text(family = "Helvetica",colour = "black",
                                          size = rel(1)), 
                 axis.line = element_line(size = 0.5,colour = "black"), 
                 axis.ticks = element_blank(),
                 panel.grid.major = element_line(colour="grey",size = rel(0.25)), 
                 panel.grid.minor = element_blank(), 
                 panel.background = element_blank(),
                 plot.title = element_text(colour = "black", face = "bold",
                                           size = rel(2),family = "Helvetica",hjust = 0.5,
                 geom_point(size = 15)))

```


# Group significance testing

## All samples
### Test for group differencesin all data: (~ env_biome+country+env_material)
Nothing treated as crossed due to complex design, could explore at a later date.
Here are some thoughts: season and field_host_name are crossed. season and Depth_m_ are crossed. field_host_name and Depth_m_ are NOT crossed, b/c cross design missing in pre-transplantion (i.e., not all combinations of host and depth)
```{r}
print(date())

# Make a data frame from the otu_data
df <- data.frame(sample_data(rt_highcover_ns_abun))
distance_methods <-c("bc_rt_highcover_ns_abun")

# Run for loop in distance matrices         
for (i in distance_methods){ 
  form <- as.formula(paste(i, "env_biome+country+env_material", sep="~"))
  print(form)
 adonis(form, data=df)->result
 print(result)
 capture.output(result, file = paste0("/Users/sofiaroitman/Google_Drive/GRAD_SCHOOL/Fall_2018/CBIOS_practicum/output/adonis/",i,"env_biome+country+env_material",date(),".txt"))
} 

```

**RESULT: biome, country, and env_material are both significant (in all distance metrics) when all data are considered together.**


### Posthoc: group differences in all data: (~ sampling_reef_name+parent_colony+transplantation_status)
```{r}
##### All samples, pairwise for biome and country

# Make a data frame from the otu_data
df <- data.frame(sample_data(rt_highcover_ns_abun))

# Pariwise Permutation MANOVAs: env_biome
pairwise.perm.manova(bc_rt_highcover_ns_abun,df$env_biome,nperm=1000)

# Pariwise Permutation MANOVAs: country
pairwise.perm.manova(bc_rt_highcover_ns_abun,df$country,nperm=1000)


###### Subsetted samples by sample type
## Make data frames from subsetted data

# Sebum
df1 <- data.frame(sample_data(rt_highcover_ns_abun_sebum))

# Mucus
df2 <- data.frame(sample_data(rt_highcover_ns_abun_mucus))

# Feces
df3 <- data.frame(sample_data(rt_highcover_ns_abun_feces))

# Saliva
df4 <- data.frame(sample_data(rt_highcover_ns_abun_saliva))


# Pariwise Permutation MANOVAs: Sebum subset, country
pairwise.perm.manova(bc_rt_highcover_ns_abun_seb,df1$country,nperm=1000)
# Pariwise Permutation MANOVAs: Sebum subset, biome
pairwise.perm.manova(bc_rt_highcover_ns_abun_seb,df1$env_biome,nperm=1000)

# Pariwise Permutation MANOVAs: Mucus subset, country
pairwise.perm.manova(bc_rt_highcover_ns_abun_muc,df2$country,nperm=1000)
# Pariwise Permutation MANOVAs: Mucus subset, biome
pairwise.perm.manova(bc_rt_highcover_ns_abun_muc,df2$env_biome,nperm=1000)

# Pariwise Permutation MANOVAs: Feces subset, country
pairwise.perm.manova(bc_rt_highcover_ns_abun_fec,df3$country,nperm=1000)
# Pariwise Permutation MANOVAs: Feces subset, biome
pairwise.perm.manova(bc_rt_highcover_ns_abun_fec,df3$env_biome,nperm=1000)

# Pariwise Permutation MANOVAs: Saliva subset, country
pairwise.perm.manova(bc_rt_highcover_ns_abun_sal,df4$country,nperm=1000)
# Pariwise Permutation MANOVAs: Saliva subset, biome
pairwise.perm.manova(bc_rt_highcover_ns_abun_sal,df4$env_biome,nperm=1000)



```
**RESULT: All countries and biomes are significantly different from each other across all sample types.**

## All samples: group signifiance testing (above) revealed  to all be significant
### Unconstrained ordination (PCoA)
[link](https://bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html)
```{r}
#### All samples: country and biome
#Calculate the PCoA on this distance matrix
rt.pcoa = ordinate(rt_highcover_ns_abun, method="PCoA", distance=bc_rt_highcover_ns_abun)
plot_scree(rt.pcoa, "Screen plot")
# PCoA ordination
pcoa<-0
pcoa <- plot_ordination(rt_highcover_ns_abun, rt.pcoa, "samples", color="country", shape = "env_biome") + ggtheme
pcoa

# Create output folder
dir.create("output/plots/")

pdf("output/plots/PCoA_rt_highcover_ns_abun_country+biome.pdf",width=9,height=5)
pcoa
dev.off()


#################
##### All samples: sample type and boime
#Calculate the PCoA on this distance matrix
rt.pcoa = ordinate(rt_highcover_ns_abun, method="PCoA", distance=bc_rt_highcover_ns_abun)
plot_scree(rt.pcoa, "Screen plot")
# PCoA ordination
pcoa<-0
pcoa <- plot_ordination(rt_highcover_ns_abun, rt.pcoa, "samples", color="env_material", shape = "env_biome") + stat_ellipse(type = "norm", lnietype = 2) + ggtheme
pcoa

# Create output folder
dir.create("output/plots/")

pdf("output/plots/PCoA_rt_highcover_ns_abun_env_material+biome.pdf",width=9,height=5)
pcoa
dev.off()

```


# PCoA based on village and country
```{r}
#### Village: country
#Calculate the PCoA on this distance matrix
rt.pcoa = ordinate(rt_highcover_ns_abun_vil, method="PCoA", distance=bc_rt_highcover_ns_abun_vil)
plot_scree(rt.pcoa, "Screen plot")
# PCoA ordination
pcoa<-0
pcoa <- plot_ordination(rt_highcover_ns_abun_vil, rt.pcoa, "samples", color="country") + ggtheme
pcoa

# Create output folder
dir.create("output/plots/")

pdf("output/plots/PCoA_rt_highcover_ns_abun_village+country.pdf",width=9,height=5)
pcoa
dev.off()


#################

#### Urban: country
#Calculate the PCoA on this distance matrix
rt.pcoa = ordinate(rt_highcover_ns_abun_urb, method="PCoA", distance=bc_rt_highcover_ns_abun_urb)
plot_scree(rt.pcoa, "Screen plot")
# PCoA ordination
pcoa<-0
pcoa <- plot_ordination(rt_highcover_ns_abun_urb, rt.pcoa, "samples", color="country") + ggtheme
pcoa

# Create output folder
dir.create("output/plots/")

pdf("output/plots/PCoA_rt_highcover_ns_abun_urban+country.pdf",width=9,height=5)
pcoa
dev.off()


```



# Violin Alpha div plots
```{r}
##### ALL SAMPLES
# Return table with selected diversity indicators
tab <- diversities(rt_highcover_ns_abun, index = "all")
kable(head(tab))

# Get metadata from pyloseq object
rt_highcover_ns_abun.meta <- meta(rt_highcover_ns_abun)
kable(head(rt_highcover_ns_abun.meta))

# Add diversity table to metadata
rt_highcover_ns_abun.meta$Shannon <- tab$shannon
rt_highcover_ns_abun.meta$InverseSimpson <- tab$inverse_simpson

write.csv(tab, file = "output/Adiv_tab.csv")

## Want to compare differences in shannon index between biome types
# Create list of pairwise comparisons
biome <- levels(rt_highcover_ns_abun.meta$env_biome) # get the variables
biome.pairs <- combn(seq_along(biome), 2, simplify = FALSE, FUN = function(i)biome[i])
print(biome.pairs)

p1 <- ggviolin(rt_highcover_ns_abun.meta, x = "env_biome", y = "Shannon", add = "boxplot", fill = "env_biome", palette = c("#fc9272", "#fa9fb5")) + ggtheme

print(p1)

p2 <- p1 + stat_compare_means(comparisons = biome.pairs, method = "wilcox.test")

print(p2)

pdf("output/plots/All_biome.pdf")
p2
dev.off()



####### SEBUM
# Return table with selected diversity indicators
tab <- diversities(rt_highcover_ns_abun_sebum, index = "all")
kable(head(tab))

# Get metadata from pyloseq object
rt_highcover_ns_abun_sebum.meta <- meta(rt_highcover_ns_abun_sebum)
kable(head(rt_highcover_ns_abun_sebum))

# Add diversity table to metadata
rt_highcover_ns_abun_sebum.meta$Shannon <- tab$shannon
rt_highcover_ns_abun_sebum.meta$InverseSimpson <- tab$inverse_simpson

## Want to compare differences in shannon index between biome types
# Create list of pairwise comparisons
sebum_biome <- levels(rt_highcover_ns_abun_sebum.meta$env_biome) # get the variables
sebum_biome.pairs <- combn(seq_along(sebum_biome), 2, simplify = FALSE, FUN = function(i)sebum_biome[i])
print(sebum_biome.pairs)

p1 <- ggviolin(rt_highcover_ns_abun_sebum.meta, x = "env_biome", y = "Shannon", add = "boxplot", fill = "env_biome", palette = c("orchid1", "mediumpurple")) + ggtheme

print(p1)

p2 <- p1 + stat_compare_means(comparisons = sebum_biome.pairs, method = "wilcox.test")

print(p2)

pdf("output/plots/Adiv_sebum_biome_wilcoxon.pdf")
p2
dev.off()



####### MUCUS
# Return table with selected diversity indicators
tab <- diversities(rt_highcover_ns_abun_mucus, index = "all")
kable(head(tab))

# Get metadata from pyloseq object
rt_highcover_ns_abun_mucus.meta <- meta(rt_highcover_ns_abun_mucus)
kable(head(rt_highcover_ns_abun_mucus))

# Add diversity table to metadata
rt_highcover_ns_abun_mucus.meta$Shannon <- tab$shannon
rt_highcover_ns_abun_mucus.meta$InverseSimpson <- tab$inverse_simpson

## Want to compare differences in shannon index between health states
# Create list of pairwise comparisons
mucus_biome <- levels(rt_highcover_ns_abun_mucus.meta$env_biome) # get the variables
mucus_biome.pairs <- combn(seq_along(mucus_biome), 2, simplify = FALSE, FUN = function(i)mucus_biome[i])
print(mucus_biome.pairs)

p3 <- ggviolin(rt_highcover_ns_abun_mucus.meta, x = "env_biome", y = "Shannon", add = "boxplot", fill = "env_biome", palette = c("mediumseagreen", "palegreen")) + ggtheme

print(p3)

p4 <- p3 + stat_compare_means(comparisons = mucus_biome.pairs, method = "wilcox.test")

print(p4)

pdf("output/plots/Adiv_mucus_biome_wilcoxon.pdf")
p4
dev.off()



####### SALIVA
# Return table with selected diversity indicators
tab <- diversities(rt_highcover_ns_abun_saliva, index = "all")
kable(head(tab))

# Get metadata from pyloseq object
rt_highcover_ns_abun_saliva.meta <- meta(rt_highcover_ns_abun_saliva)
kable(head(rt_highcover_ns_abun_mucus))

# Add diversity table to metadata
rt_highcover_ns_abun_saliva.meta$Shannon <- tab$shannon
rt_highcover_ns_abun_saliva.meta$InverseSimpson <- tab$inverse_simpson

## Want to compare differences in shannon index between health states
# Create list of pairwise comparisons
saliva_biome <- levels(rt_highcover_ns_abun_saliva.meta$env_biome) # get the variables
saliva_biome.pairs <- combn(seq_along(saliva_biome), 2, simplify = FALSE, FUN = function(i)saliva_biome[i])
print(saliva_biome.pairs)

p5 <- ggviolin(rt_highcover_ns_abun_saliva.meta, x = "env_biome", y = "Shannon", add = "boxplot", fill = "env_biome", palette = c("dodgerblue", "lightskyblue")) + ggtheme

print(p5)

p6 <- p5 + stat_compare_means(comparisons = mucus_biome.pairs, method = "t.test")

print(p4)

pdf("output/plots/Adiv_mucus_biome.pdf")
p4
dev.off()



####### FECES
# Return table with selected diversity indicators
tab <- diversities(rt_highcover_ns_abun_feces, index = "all")
kable(head(tab))

# Get metadata from pyloseq object
rt_highcover_ns_abun_feces.meta <- meta(rt_highcover_ns_abun_feces)
kable(head(rt_highcover_ns_abun_feces))

# Add diversity table to metadata
rt_highcover_ns_abun_feces.meta$Shannon <- tab$shannon
rt_highcover_ns_abun_feces.meta$InverseSimpson <- tab$inverse_simpson

## Want to compare differences in shannon index between health states
# Create list of pairwise comparisons
feces_biome <- levels(rt_highcover_ns_abun_feces.meta$env_biome) # get the variables
feces_biome.pairs <- combn(seq_along(feces_biome), 2, simplify = FALSE, FUN = function(i)feces_biome[i])
print(feces_biome.pairs)

p7 <- ggviolin(rt_highcover_ns_abun_feces.meta, x = "env_biome", y = "Shannon", add = "boxplot", fill = "env_biome", palette = c("yellow4", "khaki3")) + ggtheme

print(p7)

p8 <- p7 + stat_compare_means(comparisons = mucus_biome.pairs, method = "wilcox.test")

print(p8)

pdf("output/plots/Adiv_feces_biome_wilcoxon.pdf")
p8
dev.off()
```









