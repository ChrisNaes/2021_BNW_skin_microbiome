library(tidyverse)
library(qiime2R)
library(dplyr)
library(reshape2)
library(phyloseq)
############# Trimmed heatmap 16S ############
Physeq <- qza_to_phyloseq(
  features = "data/16S/BNW-skin-table-NoCont-wombat-COLLASPED.qza",
  tree = "data/16S/BNW-skin-sepp-tree.qza",
  taxonomy = "data/16S/BNW-skin-silva-138-taxonomy.qza",
  metadata = "data/16S/wombatid-group-sev-16S.tsv"
)

subset = subset_taxa(Physeq, Family=="Staphylococcaceae" | Family=="Streptococcaceae"
                     | Family=="Dermabacteraceae" | Family=="Corynebacteriaceae"
                     | Family=="Pseudomonadaceae" | Family=="Caulobacteraceae"
                     | Family=="Burkholderiaceae" | Family=="Brevibacteriaceae" | Family=="Sphingomonadaceae")

subset2 <- filter_taxa(subset, function(x) mean(x) > 0.01, TRUE)
subset2
subset3 <- transform_sample_counts(subset2, function(x) x/sum(x))
subset3
#create dataframe
final1 <- psmelt(subset3)
#turn all OTUs into family counts
Genus <- tax_glom(subset3, taxrank = "Genus")

subset_Genus_total = subset_taxa(Genus, Genus=="Staphylococcus" | Genus=="Brevibacterium" | Genus=="Sphingomonas"
                           | Genus=="Corynebacterium" | Genus=="Pseudomonas" | Genus=="Brachybacterium"
                           | Genus=="Breviundimonas")

subset_Genus = subset_taxa(Genus, Genus=="Staphylococcus" | Genus=="Brevibacterium"
                           | Genus=="Corynebacterium" | Genus=="Brachybacterium"
                           | Genus=="Breviundimonas")

Genus_df <- psmelt(subset_Genus) # create dataframe from phyloseq object

#### Graphing #####

plot_gen <- ggplot(Genus_df, aes(x=Sample, y=Genus, fill = Abundance)) +
  geom_tile()
plot_gen
plot_gen2 <- plot_gen+xlab("Wombat ID")+
  ylab("Genus")+
  labs(fill = "Relative abundance")+
  facet_grid(~Severity, scales = "free", space = "free_x")+
  theme(axis.text.x = element_text(size = 30, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.title = element_text(size = 40, face = "bold"),
        legend.text = element_text(size = 30),
        legend.title.align = 0.5,
        legend.title = element_text(size = 30, face = "bold"),
        strip.text = element_text(size = 30, face = "bold"))

heat_16S <- plot_gen2 + scale_fill_distiller(palette = "Blues", direction = 1) +
  guides(fill = guide_colorbar(barheight = 11))
heat_16S

ggsave(filename = "output/Fig3_heat_16S.svg", width = 20, height = 5, dpi = 300)

# Trimmed heatmap ITS2
PhyseqITS <- qza_to_phyloseq(
  features = "data/ITS2/BNW-skin-table-NoCont-wombat-COLLASPED.qza",
  tree = "data/ITS2/unite-rooted-tree.qza",
  taxonomy = "data/ITS2/BNW-skin-unite-taxonomy-paired.qza",
  metadata = "data/ITS2/wombatid-group-sev-ITS2.tsv"
)
subset = subset_taxa(PhyseqITS, Family=="Pseudeurotiaceae" | Family=="Didymellaceae"
                     | Family=="Debaryomycetaceae" | Family=="Cladosporiaceae"
                     | Family=="Bulleribasidiaceae" | Family=="Aspergillaceae" )
subset2 <- filter_taxa(subset, function(x) mean(x) > 0.01, TRUE)
subset2
subset3 <- transform_sample_counts(subset2, function(x) x/sum(x))
subset3
#create dataframe
final <- psmelt(subset3)
#turn all OTUs into family counts
Genus <- tax_glom(subset3, taxrank = "Genus")
Genus
subset_Genus_total = subset_taxa(Genus, Genus=="Vishniacozyma" | Genus=="Pseudogymnoascus" | Genus=="Penicillium"
                           | Genus=="Neoascochyta" | Genus=="Debaryomyces" | Genus=="Cladosporium")

subset_Genus = subset_taxa(Genus, Genus=="Pseudogymnoascus" | Genus=="Debaryomyces" | Genus=="Cladosporium")

Genus_df <- psmelt(subset_Genus) # create dataframe from phyloseq object

#Make plots
plot1 <- ggplot(Genus_df, aes(x=Sample, y=Genus, fill = Abundance)) +
  geom_tile()
plot1
plot2 <- plot1+xlab("Wombat ID") +
  ylab("Genus") +
  labs(fill = "Relative abundance") +
  facet_grid(~Severity, scales = "free", space = "free_x") +
  theme(axis.text.x = element_text(size = 30, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 30, face = "bold"),
        axis.title = element_text(size = 40, face = "bold"),
        legend.text = element_text(size = 30),
        legend.title.align = 0.5,
        legend.title = element_text(size = 30, face = "bold"),
        strip.text = element_text(size = 30, face = "bold")
  )

plot2 + scale_fill_distiller(palette = "Blues", direction = 1) +
  guides(fill = guide_colorbar(barheight = 11))


ggsave(filename = "output/Fig_heatmap_ITS2_Subset_Genus_relative_abun_.svg", width = 20, height = 5, dpi = 300)
