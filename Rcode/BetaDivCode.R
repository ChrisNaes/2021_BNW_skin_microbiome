rm(list=ls())#clears the crap

library("dplyr")
library("ggpubr")
library("vegan")
library("tidyr")
library("knitr")
library("tidyverse")
library("svglite")
library("ggplot2")
library("qiime2R")

########################################################################
## Beta diversity stat test ##


##16S##
#######

metadata <- read_tsv("data/16S/BNW-skin-metadata-16S_Final.tsv")
dat <- read.csv("data2.csv") # bacterial data
head(metadata)
head(dat)

dat=dat[order(dat$SampleID),]
metadata=metadata[order(metadata$SampleID),]

dat$SampleID
metadata$SampleID

metadat=data.frame(metadata[,c(1:2,8)])
metadat

dat1=cbind(metadat,dat)
dim(dat1)
head(dat1)

## Permanova with wombat id nested within severity
permanova=adonis(dat1[,5:351]~dat1$Severity/dat1$wombatid)
permanova

##Subset the data to evaluate pairwise coparisions
HvM=dat1[-which(dat1$Severity=="Severe mange"),]
HvSM=dat1[-which(dat1$Severity=="Mangy"),]
MvSM=dat1[-which(dat1$Severity=="Confidently healthy"),]

##Pairwise PERMANOVAs
adonis(HvM[,5:351]~HvM$Severity/HvM$wombatid)
adonis(HvSM[,5:351]~HvSM$Severity/HvSM$wombatid)
adonis(MvSM[,5:351]~MvSM$Severity/MvSM$wombatid)

###########################
##Graphing Bacterial data##
###########################

metadata <- read_tsv("data/16S/BNW-skin-metadata-16S_Final.tsv")
dat <- read.csv("data2.csv")
head(dat)

dim(dat)
dim(metadata)

dat=dat[order(dat$SampleID),]
metadata=metadata[order(metadata$SampleID),]

dat$SampleID
metadata$SampleID

dat=dat[,2:dim(dat)[2]]
metadat=data.frame(metadata[,c(1:2,8)])
metadat

#initial MDS plot to check it works
nmds=metaMDS(dat,k=2)

#better MDS plot

cols=ifelse(metadat$Severity=="Severe mange","#8DA0CB",
            ifelse(metadat$Severity=="Mangy","#FC8D62",
                   ifelse(metadat$Severity=="Confidently healthy","#66C2A5",0)))
ordiplot(nmds,type="n")
orditorp(nmds,display="sites",labels=metadat$wombatid,col=cols, cex = 1)
p1 <- ordispider(nmds,display="sites",groups=metadat$Severity,col="grey")
par(mfrow=c(1,2))
nmds=metaMDS(dat,k=2)
cols=ifelse(metadat$Severity=="Severe mange","#8DA0CB",
            ifelse(metadat$Severity=="Mangy","#FC8D62",
                   ifelse(metadat$Severity=="Confidently healthy","#66C2A5",0)))
ordiplot(nmds,type="n")
orditorp(nmds,display="sites",labels=metadat$wombatid,col=cols, cex = 1)
p1 <- ordispider(nmds,display="sites",groups=metadat$Severity,col="grey")

########################################################################
## Beta diversity stat test ##

##ITS##
#######

metadata <- read_tsv("BNW-skin-metadata-ITS2_Final.tsv")
dat <- read.csv("data_fungi.csv")
head(metadata)
head(dat)

dat=dat[order(dat$SampleID),]
metadata=metadata[order(metadata$SampleID),]

dat$SampleID
metadata$SampleID

metadat=data.frame(metadata[,c(1:2,8)])
metadat

dat1=cbind(metadat,dat)
dim(dat1)
head(dat1)

## Permanova with wombat id nested within severity
permanova=adonis(dat1[,5:184]~dat1$Severity/dat1$wombatid)
permanova

##Subset the data to evaluate pairwise coparisions
HvM=dat1[-which(dat1$Severity=="Severe mange"),]
HvSM=dat1[-which(dat1$Severity=="Mangy"),]
MvSM=dat1[-which(dat1$Severity=="Confidently healthy"),]

##Pairwise PERMANOVAs
adonis(HvM[,5:184]~HvM$Severity/HvM$wombatid)
adonis(HvSM[,5:184]~HvSM$Severity/HvSM$wombatid)
adonis(MvSM[,5:184]~MvSM$Severity/MvSM$wombatid)

###########################
##Graphing Fungal data##
###########################

metadata1 <- read_tsv("data/ITS2/BNW-skin-metadata-ITS2_Final.tsv")
dat1 <- read.csv("data_fungi.csv")
head(dat1)

dim(dat1)
dim(metadata1)

dat1=dat1[order(dat1$SampleID),]
metadata1=metadata1[order(metadata1$SampleID),]

dat1$SampleID
metadata1$SampleID

dat1=dat1[,2:dim(dat1)[2]]
metadat1=data.frame(metadata1[,c(1:2,8)])
metadat1

#initial MDS plot to check it works
nmds1=metaMDS(dat1,k=2)

#better MDS plot

cols=ifelse(metadat1$Severity=="Severe mange","#8DA0CB",
            ifelse(metadat1$Severity=="Mangy","#FC8D62",
                   ifelse(metadat1$Severity=="Confidently healthy","#66C2A5",0)))
ordiplot(nmds1,type="n")
orditorp(nmds1,display="sites",labels=metadat1$wombatid,col=cols, cex = 1)
p2 <- ordispider(nmds1,display="site",groups=metadat1$Severity,col="grey")
