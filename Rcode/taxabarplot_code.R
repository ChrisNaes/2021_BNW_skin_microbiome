library("phyloseq")
library("qiime2R")
library("viridis")
library("RColorBrewer")
library("ggplot2")

########## Bacterial Stacked Taxa Barplot ############
Physeq_16S <- qza_to_phyloseq(
  features = "data/16S/BNW-skin-table-NoCont-wombat-COLLASPED.qza",
  tree = "data/16S/BNW-skin-sepp-tree.qza",
  taxonomy = "data/16S/BNW-skin-silva-138-taxonomy.qza",
  metadata = "data/16S/wombatid-group-sev-16S.tsv"
)

merge_low_abundance16S <- function(Physeq_16S, threshold=0.001){
  transformed <- transform_sample_counts(Physeq_16S, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.list <- row.names(otu.table[rowMeans(otu.table) < threshold,])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Less abundant"
      tax_table(merged)[i,1:7] <- "Less abundant"}
  }
  return(merged)
}

merge_less_than_top16S <- function(Physeq_16S, top=10){
  transformed <- transform_sample_counts(Physeq_16S, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Less abundant"
      tax_table(merged)[i,1:7] <- "Less abundant"}
  }
  return(merged)
}

ps.fam16S <- tax_glom(Physeq_16S, "Family")
ps.fam.top1016S <- merge_less_than_top16S(ps.fam16S, top=10)

glom_16S <- tax_glom(ps.fam.top1016S, taxrank = 'Family')
glom_16S # should list # taxa as # phyla
data_glom_16S <- psmelt(glom_16S) # create dataframe from phyloseq object
data_glom_16S$Family <- as.character(data_glom_16S$Family) #convert to character

#Count # phyla to set color palette
Count_16S = length(unique(data_glom_16S$Family))
Count_16S

data_glom_16S$Family <- factor(data_glom_16S$Family,
                               levels = c("Microbacteriaceae", "Sphingomonadaceae", "Pseudomonadaceae",
                                          "Sphingobacteriaceae", "Staphylococcaceae", "Alcaligenaceae",
                                          "Brevibacteriaceae", "Corynebacteriaceae", "Micrococcaceae",
                                          "Caulobacteraceae", "Less abundant"))

plot <- ggplot(data=data_glom_16S, aes(x=Sample, y=Abundance, fill=Family)) +
  facet_grid(~Severity, scales = "free") + xlab("Wombat ID") + ylab("Relative Abundance")

plot1 <- plot +
  geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=5)) +
  theme(axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 40, face = "bold"),
        legend.text = element_text(size = 30),
        legend.title.align = 0.5,
        legend.title = element_text(size = 30, face = "bold", angle = 90),
        strip.text = element_text(size = 30, face = "bold"))

plot1 +  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                                      "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                                      "#FFFF99", "#B15928"))

ggsave(filename = "output/Fig_taxa_barplot-Family_Top10_16S.svg", height=11, width=20, dpi = 300)

########## Fungal STacked Taxa Barplot ##############
PhyseqITS <- qza_to_phyloseq(
  features = "data/ITS2/BNW-skin-table-NoCont-wombat-COLLASPED.qza",
  tree = "data/ITS2/unite-rooted-tree.qza",
  taxonomy = "data/ITS2/BNW-skin-unite-taxonomy-paired.qza",
  metadata = "data/ITS2/wombatid-group-sev-ITS2.tsv"
)

merge_low_abundance <- function(PhyseqITS, threshold=0.001){
  transformed <- transform_sample_counts(PhyseqITS, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.list <- row.names(otu.table[rowMeans(otu.table) < threshold,])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Less abundant"
      tax_table(merged)[i,1:7] <- "Less abundant"}
  }
  return(merged)
}

merge_less_than_top <- function(PhyseqITS, top=10){
  transformed <- transform_sample_counts(PhyseqITS, function(x) x/sum(x))
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Less abundant"
      tax_table(merged)[i,1:7] <- "Less abundant"}
  }
  return(merged)
}

brewer.pal(n = 12, name = "Paired")

ps.fam <- tax_glom(PhyseqITS, "Family")
ps.fam.top10 <- merge_less_than_top(ps.fam, top=10)

glom <- tax_glom(ps.fam.top10, taxrank = 'Family')
glom # should list # taxa as # phyla
data_glom<- psmelt(glom) # create dataframe from phyloseq object
data_glom$Family <- as.character(data_glom$Family) #convert to character

#Count # phyla to set color palette
Count = length(unique(data_glom$Family))
Count

data_glom$Family <- factor(data_glom$Family, levels = c("Cladosporiaceae", "Pseudeurotiaceae", "Debaryomycetaceae",
                                                      "Umbelopsidaceae", "Didymellaceae", "Rhynchogastremataceae",
                                                      "Aspergillaceae", "Bulleribasidiaceae",
                                                      "Phaeosphaeriaceae", "Myxotrichaceae", "Less abundant"))

plot <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill=Family)) +
  facet_grid(~Severity, scales = "free") + xlab("Wombat ID") + ylab("Relative Abundance")

plot1 <- plot +
  geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=5)) +
  theme(axis.text.x = element_text(size = 30, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 30),
        axis.title = element_text(size = 40, face = "bold"),
        legend.text = element_text(size = 30),
        legend.title.align = 0.5,
        legend.title = element_text(size = 30, face = "bold", angle = 90),
        strip.text = element_text(size = 30, face = "bold"))

plot1 +  scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                                      "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                                      "#FFFF99", "#B15928"))

ggsave(filename = "output/Fig_taxa_barplot-Family_Top10_ITS2.svg", height=11, width=20, dpi = 300)
