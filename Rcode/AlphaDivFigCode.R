library("phyloseq")
library("qiime2R")
library("tidyverse")

# Alpha diversity 16S
meta_16S <- read_tsv("data/16S/BNW-skin-metadata-16S_Final.tsv")
OTUs <- read_qza("data/16S/observed_features_vector.qza")
OTUs <- OTUs$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged
gplots::venn(list(metadata=meta_16S$SampleID, observed_features=OTUs$SampleID))

meta_16S<-
  meta_16S %>%
  left_join(OTUs)
head(meta_16S)

level_order <- c("Confidently healthy", "Mangy", "Severe mange")

# overview of number of wombats per category
sample_size <- meta_16S %>%
  group_by(Severity, wombatid) %>%
  summarize(num=n())

df <- meta_16S %>%
  group_by(Severity) %>%
  summarise(num=n())

df1 <- meta_16S %>%
  group_by(Severity, observed_features) %>%
  summarise(num=n())

fig_OTUs_16S <- ggplot(meta_16S, aes(x = factor(Severity, level = level_order), y = observed_features))

cbp2_adiv2 <- palette("Set2")

otu <- fig_OTUs_16S +
  #Boxplot
  geom_boxplot(size=2, outlier.shape=5, outlier.size=3, outlier.stroke=3, aes(colour=Severity)) +
  #Jitter, size, colour
  geom_jitter(position=position_dodge2(0.3), size=5, aes(colour=Severity)) +
  #geom_text(data = df, aes(label = c("a", "b", "b"),
  #                               x = Severity, y = c(300,150,175)), size = 15) +
  #geom_signif(data = res.tukey, comparisons = list(c("Confidently healthy", "Mangy", "Severe mange")),
  #            map_signif_level = TRUE) +
  #Custom manual colours
  scale_colour_manual(values=cbp2_adiv2) +
  #Tick labels
  theme(axis.text.x = element_text(face="bold", size=30, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", size=30),
        axis.title.x = element_text(size=30, face="bold"),
        axis.title.y = element_text(size=30, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        #Legend
        legend.position = "none") +
  labs(x = "\nSeverity") +
  labs(y = "Bacterial 16S ASV richness\n")
otu1 <- otu + theme(plot.title = element_text(size = 40, face = "bold"),
                    axis.text = element_text(size = 40),
                    axis.title = element_text(size = 40, face = "bold"),
                    legend.text = element_text(size = 40),
                    legend.title = element_text(size = 20))
otu1 <- otu1 + scale_x_discrete(labels=c("Conf. healthy, n=6", "Mangy, n=3", "Severe mange, n=4"))

g1 <- data.frame(a = c(1, 1:3,3), b = c(350,351,351,351,350))
g2 <- data.frame(a = c(1, 1,2, 2), b = c(300,301,301,300))
g3 <- data.frame(a = c(2, 2,3, 3), b = c(200,201,201,200))

otu1 <- otu1 +
  geom_path(data=g1, aes(x=a, y=b)) + annotate("text", x=1.5, y=368, label = "**", size = 10)+
  geom_path(data=g2, aes(x=a, y=b)) + annotate("text", x=1.5, y=315, label = "*", size = 10)+
  geom_path(data=g3, aes(x=a, y=b)) + annotate("text", x=2.5, y=215, label = "ns", size = 10)
##Save image as .svg
ggsave(filename = "output/Fig_AlphaDiv-OTUs-16S-wombat-Final-samplesize2.svg", width = 20, height = 11, dpi = 300)

# Alpha diversity ITS2
meta_ITS2 <- read_tsv("data/ITS2/BNW-skin-metadata-ITS2_Final.tsv")
OTUs_ITS <- read_qza("data/ITS2/observed_features_vector.qza")
OTUs_ITS <- OTUs_ITS$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged
gplots::venn(list(metadata=meta_ITS2$SampleID, observed_features=OTUs_ITS$SampleID))

meta_ITS2<-
  meta_ITS2 %>%
  left_join(OTUs_ITS)
head(meta_ITS2)

level_order <- c("Confidently healthy", "Mangy", "Severe mange")
fig_OTUs_ITS <- ggplot(meta_ITS2, aes(x = factor(Severity, level = level_order), y =observed_features))

sample_size <- meta_ITS2 %>%
  group_by(Severity, wombatid) %>%
  summarize(num=n())

otu_ITS <- fig_OTUs_ITS +
  #Boxplot
  geom_boxplot(size=2, outlier.shape=5, outlier.size=3, outlier.stroke=3, aes(colour=Severity)) +
  #Jitter, size, colour
  geom_jitter(position=position_dodge2(0.3), size=5, aes(colour=Severity)) +
  #Custom manual colours
  scale_colour_manual(values=cbp2_adiv2) +
  #Tick labels
  theme(axis.text.x = element_text(face="bold", size=30, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", size=30),
        axis.title.x = element_text(size=30, face="bold"),
        axis.title.y = element_text(size=30, face="bold"),
        axis.line = element_line(colour = "black"),
        #Background panel
        panel.background = element_rect(fill = "White"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        #Legend
        legend.position = "none") +
  #legend.title = element_text(size=0),
  #legend.text = element_text(size=0),
  #legend.key = element_rect(fill = "white", color = NA),
  #legend.key.size = unit(0, "line")) +
  #Axis labels
  #scale_y_continuous(breaks=seq(0,400,25)) +
  labs(x = "\nSeverity") +
  labs(y = "Fungal ITS2 ASV richness\n")
otu2 <- otu_ITS + theme(plot.title = element_text(size = 20, face = "bold"),
                        axis.text = element_text(size = 15),
                        axis.title = element_text(size = 18, face = "bold"),
                        legend.text = element_text(size = 20),
                        legend.title = element_text(size = 20))
otu2 <- otu2 + scale_x_discrete(labels=c("Conf. healthy, n=6", "Mangy, n=3", "Severe mange, n=4"))

g4 <- data.frame(a = c(1, 1:3,3), b = c(170,171,171,171,170))
g5 <- data.frame(a = c(1, 1,2, 2), b = c(150,151,151,150))
g6 <- data.frame(a = c(2, 2,3, 3), b = c(100,101,101,100))

otu2 <- otu2 +
  geom_path(data=g4, aes(x=a, y=b)) + annotate("text", x=1.5, y=179, label = "ns", size = 10)+
  geom_path(data=g5, aes(x=a, y=b)) + annotate("text", x=1.5, y=158, label = "ns", size = 10)+
  geom_path(data=g6, aes(x=a, y=b)) + annotate("text", x=2.5, y=107, label = "ns", size = 10)

##Save image as .svg
ggsave(filename = "output/Fig_AlphaDiv-OTUs-ITS-wombat-Final-samplesize2.svg", width = 20, height = 11, dpi = 300)

library("ggpubr")
OTU_com <- ggarrange(otu1 + rremove("xlab") +
                       theme(axis.text = element_text(size = 20)),
                     otu2 + rremove("xlab") +
                       theme(plot.margin = margin(r = 1)),
                     labels = c("A", "B"),
                     font.label = list(size = 40, face = "bold"),
                     align = "v"
)
OTU_com
ggsave(filename = "output/Fig2_AlphaDiv_otu_combined-samplesize.svg", width = 20, height = 11, dpi = 300)
