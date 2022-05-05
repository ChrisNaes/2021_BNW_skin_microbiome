
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
library("car")

###################################################################
## Alpha diversity stat test ##
##Data entry, sorting and processing for graphing and analyses##

##16S##
#######

metadata <- read_tsv("data/16S/BNW-skin-metadata-16S_Final.tsv")
dat <- read.csv("data/16S/OTUs-16S.csv")
head(dat)
head(metadata)

dim(dat)
dim(metadata)

##Simple data merging and reducing the dataset to focal variables
dat1=merge(dat,metadata,by=c("SampleID"),all=T)
head(dat1)
colnames(dat1)
dat2=dat1[,c(1,3:5,10)]
colnames(dat2)=c("SampleID","OTU",colnames(dat2[3:5]))
head(dat2)

##Box and whisker plot
plot(factor(dat2$Severity),dat2$OTU) #simple box plot

##Repeated measures ANOVA with wombatid as the repeated measure
mod=aov(OTU~Severity+Error(wombatid),data=dat2)
summary(mod)

##Tukey's posthoc test for Severity levels - the emmeans package allows you to do this when you have repeated measures and can't use TukeyHSD
library(emmeans)
options(contrasts = c("contr.sum", "contr.poly"))# set orthogonal contrasts
emm <- emmeans(mod, ~ Severity)
pairs(emm)  # adjust argument not specified -> default p-value adjustment in this case is "tukey"

##Levene Test - Homogeneity of variance##

leveneTest(OTU~Severity, data = dat2) #overall
leveneTest(OTU~Severity, data = dat2[-which(dat2$Severity=="Severe mange"),]) #healthy v mangy
leveneTest(OTU~Severity, data = dat2[-which(dat2$Severity=="Mangy"),]) #healthy v severe
leveneTest(OTU~Severity, data = dat2[-which(dat2$Severity=="Confidently healthy"),]) #mangy v severe.

##ITS##
#######

metadata <- read_tsv("data/ITS2/BNW-skin-metadata-ITS2_Final.tsv")
dat <- read.csv("data/ITS2/OTUs-ITS2.csv")
head(dat)
head(metadata)

dim(dat)
dim(metadata)

##Simple data merging and reducing the dataset to focal variables
dat1=merge(dat,metadata,by=c("SampleID"),all=T)
head(dat1)
colnames(dat1)
dat2=dat1[,c(1,3:5,10)]
colnames(dat2)=c("SampleID","OTU",colnames(dat2[3:5]))
head(dat2)

##Box and whisker plot
plot(factor(dat2$Severity),dat2$OTU) #simple box plot

##Repeated measures ANOVA with wombatid as the repeated measure
mod=aov(OTU~Severity+Error(wombatid),data=dat2)
summary(mod)

##Tukey's posthoc test for Severity levels - the emmeans package allows you to do this when you have repeated measures and can't use TukeyHSD
library(emmeans)
options(contrasts = c("contr.sum", "contr.poly"))# set orthogonal contrasts
emm <- emmeans(mod, ~ Severity)
pairs(emm)  # adjust argument not specified -> default p-value adjustment in this case is "tukey"


##Levene Test - Homogeneity of variance##

leveneTest(OTU~Severity, data = dat2) #overall
leveneTest(OTU~Severity, data = dat2[-which(dat2$Severity=="Severe mange"),]) #healthy v mangy
leveneTest(OTU~Severity, data = dat2[-which(dat2$Severity=="Mangy"),]) #healthy v severe
leveneTest(OTU~Severity, data = dat2[-which(dat2$Severity=="Confidently healthy"),]) #mangy v severe


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



##ITS##
#######

metadata <- read_tsv("data/ITS2/BNW-skin-metadata-ITS2_Final.tsv")
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



