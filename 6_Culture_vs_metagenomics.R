---
# title: "PIPOVAP_AK_10/03/2022 Taxonomic composition"
# date: "10/03/2022"
# output: html_document
---
  
# Libraries

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
library(BiocManager)
install.packages("phyloseq")
library("phyloseq")
# BiocManager::install("phyloseq")
BiocManager::install("GenomeInfoDb",force=TRUE)
library(GenomeInfoDb)
install.packages('rlang')
library(rlang)
library(tidyverse)
library(readr)
library(ggplot2)
library(hablar)
library(lubridate)
library(magrittr)
library(RColorBrewer)
library(backports)
library(ggpubr)
library(scales)
library(openxlsx)
library(remotes)
library(rmarkdown)
library(VennDiagram)
library(svglite)
install.packages("GUniFrac")
library(GUniFrac)
install.packages("devtools")
library(devtools) # Load the devtools package
install_github("microbiome/microbiome")
library(microbiome)
devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
remotes::install_github("jfq3/ggordiplots")
library(ggordiplots)
remotes::install_github("MadsAlbertsen/ampvis2")
library(ampvis2)
BiocManager::install("DESeq2")
library(DESeq2)
install.packages("readxl")
library("readxl")
BiocManager::install("miaViz")
library(mia)
library(miaViz)
library(tidyselect)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
remotes::install_github("gmteunisse/Fantaxtic")
library(fantaxtic)
remotes::install_github("david-barnett/microViz")
library(microViz)


# Loading data ----

# Load Lung Data with No merge, Rarefaction = 10000, No collapse
physeqLung1<-readRDS("//file3.intranet.chuv/data3/SHARE/PIPOVAP/220203_pipovap_lung/220203_lung_10k.rds")
OTU_lung1 <- otu_table(physeqLung1)
sam_lung1 <- sample_data(physeqLung1)
tax_lung1 <- tax_table(physeqLung1)
tree_lung1 <- phy_tree(physeqLung1)
refseq_lung1 <- refseq(physeqLung1)

head(rownames(OTU_lung1))
head(colnames(OTU_lung1))

head(rownames(sam_lung1))
head(colnames(sam_lung1))

head(rownames(tax_lung1))
head(colnames(tax_lung1))

# Load filtered Lung metadata 
metadata_lung <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/220203_pipovap_lung/metadata_lung.xlsx")
View(metadata_lung)

metadata_lung<-metadata_lung%>%dplyr::rename(number_of_analyses=...1) #rename column "...1" in the metadata file
metadata_lung<-metadata_lung%>%dplyr::rename(smpl=Sample) 


names(metadata_lung)

spec(metadata_lung)
class(metadata_lung)
head(metadata_lung)
colnames(metadata_lung)

View(metadata_lung)

# Vector of Columns to be converted into factors 
ColAsFactors <- colnames(metadata_lung)[c(1:12)]
metadata_lung <- metadata_lung %>% hablar::convert(fct(all_of(ColAsFactors)))


## Production of Heatmaps to compare abundance on inclusion vs D1 of infection of bacteria isolated in conventional cultures on D1 of infection

# Load metadata containing microbiological information
micro<-read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/6.Comparison of metagenomics with conventional cultures/microbiology.xlsx")
View(micro)

#subset "micro" dataframe to exclude IPP and Group columns (information on Group is already present in metadata_lung dataframe)
micro<-(select(micro,!c(ipp,Group)))

#Rename record_id column to match the metadata_lung record_ID column
micro<-dplyr::rename(micro,"record_ID"="record_id")
View(micro)

#Join the tow dataframes (metadata_lung and micro dataframes based on record_ID column)
class(micro$record_ID)
class(metadata_lung$record_ID)
micro$record_ID<-as.factor(micro$record_ID)

metadata_lung<-left_join(metadata_lung,micro,by="record_ID")
View(metadata_lung)

# Use phyloseq::sample_data function in view of preparing a phyloseq object
sam_lung2 <- sample_data(metadata_lung)
View(sam_lung2)


# For consistency with the other elements of phyloseq object, use "sample_label"
# as sample_names

sample_names(sam_lung2) <- as.factor(sam_lung2$sample_label)
head(sample_names(sam_lung2))

# Create new phyloseq object with sam_lung2
physeqLung2 <- phyloseq(OTU_lung1, tax_lung1, sam_lung2, refseq_lung1, tree_lung1)
head(sample_names(physeqLung2))
sample_variables(physeqLung2)

# Correct the order of levels for different variables if necessary
levels(get_variable(physeqLung2, "Filtered_Data_AK_Dec2020"))

levels(get_variable(physeqLung2, "Subgroup"))

Correct_order_filtered_data <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation")
Correct_order_subgroup <- c("Inclusion","Infection_D5", "Discharge")

# Edit names of the levels
Correct_labels <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation")

sample_data(physeqLung2)$Filtered_Data_AK_Dec2020 <- factor(sample_data(physeqLung2)$Filtered_Data_AK_Dec2020,
                                                            levels = Correct_order_filtered_data,
                                                            labels = Correct_labels)
levels(get_variable(physeqLung2, "Filtered_Data_AK_Dec2020"))

sample_data(physeqLung2)$Subgroup <- factor(sample_data(physeqLung2)$Subgroup,
                                            levels = Correct_order_subgroup)

levels(get_variable(physeqLung2,"Subgroup"))


# Set manual scales of colors for the following groups: 1)Filtered_data_Ak_Dec2020, 2)Subgroup, 3)Group
InclInfectColors4 <- c("Inclusion"="gold2", "Infection_D1"="coral2",  #color palette for Filtered_data_AK_Dec2020
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3")


InclInfectColors3 <- c("Inclusion"="gold2", "Infection_D5"="chartreuse4",   # Color palette for Subgroup
                       "Discharge"="cadetblue3")

brewer.pal(n = 12, name = "Set3")                     

Infectcolorgroup<-c("Control"="#8DD3C7","Other_infection"="#FDB462","Pneumonia"="#80B1D3")  # color palette for Groups



## Now we have our phyloseq object "physeqLung2" containing our data of interest

# We are now going to produce phylogenetic trees showing the above mentioned species to see if there are any other closely phylogenetically related species that we have to take into account.


# Transform counts to relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqLung2Rel <- transform_sample_counts(physeqLung2, count_to_rel_abund)
View(psmelt(physeqLung2))
View(psmelt(physeqLung2Rel))



### Subset to Genus or Family level the bacteria found in cultures to see how they are classified in phylogenetic trees in order
### to select ASVs of interest to search for presence/absence


#### ENTEROBACTERIACEAE ####

#subset Proteobacteria phylum to see how the Enterobacteriaceae family is classified on the phylogenetic tree
proteobacteria_ps <- subset_taxa(physeqLung2, Phylum == "Proteobacteria")
View(tax_table(proteobacteria_ps))
taxa_sums(proteobacteria_ps)
otu_table(proteobacteria_ps)

ps_enterobac <- subset_taxa(physeqLung2Rel, Family == "Enterobacteriaceae")
View(tax_table(ps_enterobac))

ps_enterobac <- subset_taxa(physeqLung2Rel, Class == "Gammaproteobacteria")
View(tax_table(ps_enterobac))


# Plot tree of Enterobacteriaceae family 

#First color based on species
tree_enterobac_1<-plot_tree(ps_enterobac,nodelabf=nodeplotboot(),label.tips="Species",
          color="Species",
          size="Abundance",
          sizebase=10,
          plot.margin=0.1,
          ladderize="right")

tree_enterobac_1

ggsave("tree_enterobacteriaceae.pdf",
       width = 60, height = 45, dpi = 300, units = "cm")


#Second color based on confidence of identification
tree_enterobac_2<-plot_tree(ps_enterobac,nodelabf=nodeplotboot(),label.tips="Species",
          color="Confidence",
          size="Abundance",
          sizebase=10,
          plot.margin=0.1,
          ladderize="right")

tree_enterobac_2

ggsave("tree_enterobacteriaceae_conf.tiff",
       width = 60, height = 45, dpi = 300, units = "cm")


ggarrange(tree_enterobac_1,tree_enterobac_2, labels = c("A", "B"),common.legend = FALSE, legend = "right" )

ggsave("tree_enterobac.tiff",
       width = 60, height = 45, dpi = 300, units = "cm")



#### Morganellaceae ####

## Subset phyloseq object for Morganellaceae family

morgan<-subset_taxa(physeqLung2Rel, Family == "Morganellaceae")
View(tax_table(morgan))


# Plot tree of Aeromonadaceae family 

#First color based on species
plot_tree(morgan,nodelabf=nodeplotboot(),label.tips="Species",
          color="Species",
          size="Abundance",
          sizebase=10,
          plot.margin=0.1,
          ladderize="right")



#### AEROMONADACEAE ####

## Subset phyloseq object for Aeromonadaceae family

aeromon<-subset_taxa(physeqLung2Rel, Family == "Aeromonadaceae")
View(tax_table(aeromon))


# Plot tree of Aeromonadaceae family 

#First color based on species
plot_tree(aeromon,nodelabf=nodeplotboot(),label.tips="Species",
                            color="Species",
                            size="Abundance",
                            sizebase=10,
                            plot.margin=0.1,
                            ladderize="right")




#### STREPTOCOCCACEAE ####

## Subset phyloseq object for Streptococcaceae family

strepto<-subset_taxa(physeqLung2Rel, Species == "Streptococcus pneumoniae")
View(tax_table(strepto))

#First color based on species
plot_tree(strepto,nodelabf=nodeplotboot(),label.tips="Species",
          color="Species",
          size="Abundance",
          sizebase=10,
          plot.margin=0.1,
          ladderize="right")




#### STAPHYLOCCOCCACEAE ####


#subset Firmicutes phylum to see how the Staphyloccoccaceae family is classified on the phylogenetic tree

staphylo_ps <- subset_taxa(physeqLung2, Phylum == "Firmicutes")
View(tax_table(staphylo_ps))
taxa_sums(staphylo_ps)
otu_table(staphylo_ps)

ps_staph <- subset_taxa(physeqLung2Rel, Family == "Staphylococcaceae")
View(tax_table(ps_staph))


# Plot tree of Staphylococcaceae family 

#First color based on species
tree_staph_1<-plot_tree(ps_staph,nodelabf=nodeplotboot(),label.tips="Species",
                        color="Species",
                        size="Abundance",
                        sizebase=10,
                        plot.margin=0.1,
                        ladderize="right")

tree_staph_1

ggsave("tree_staph.tiff",
       width = 60, height = 45, dpi = 300, units = "cm")



#### PASTEURELACEAE ####

#subset Proteobacteria phylum to see how the Pasteurelaceae family is classified on the phylogenetic tree


proteobacteria_ps <- subset_taxa(physeqLung2, Phylum == "Proteobacteria")
View(tax_table(proteobacteria_ps))
taxa_sums(proteobacteria_ps)
otu_table(proteobacteria_ps)

ps_pasteurel <- subset_taxa(physeqLung2Rel, Family == "Pasteurellaceae")
View(tax_table(ps_pasteurel))


# Plot tree of Enterobacteriaceae family 

#First color based on species
tree_pasteur_1<-plot_tree(ps_pasteurel,nodelabf=nodeplotboot(),label.tips="Species",
                          color="Species",
                          size="Abundance",
                          sizebase=10,
                          plot.margin=0.1,
                          ladderize="right")

tree_pasteur_1

ggsave("tree_pasteurellaceae.tiff",
       width = 60, height = 45, dpi = 300, units = "cm")





## Now we can subset our data based on record_ID and the identified germ and produce heatmaps on inclusion and D1 of infection

## Patient 1 : Infection with Citrobacter freundii

## Now we can subset our data based on record_ID and the identified germ and produce heatmaps on inclusion and D1 of infection

p1<-subset_samples(physeqLung2,record_ID=="1")
p1<-microbiome::transform(p1, "compositional")
View(psmelt(p1))

p1df<-psmelt(p1)
p1df<-p1df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p1df)


#Plot heatmap

p1_taxbin<-psmelt(tax_transform(p1, trans = "binary", rank = "Species")%>%otu_get())

View(p1_taxbin)
p1_taxbin<-p1_taxbin%>%filter(OTU=="Pseudocitrobacter_sp")

class(p1_taxbin)
class(p1_taxbin$Sample)
class(p1_taxbin$Abundance)
class(p1_taxbin$OTU)
p1_taxbin$OTU<-as.factor(p1_taxbin$OTU)
p1_taxbin$Abundance<-as.factor(p1_taxbin$Abundance)
p1_taxbin$Sample<-as.factor(p1_taxbin$Sample)
levels(p1_taxbin$Sample)
levels(p1_taxbin$Abundance)


p1_taxbin$Sample<-plyr::revalue(p1_taxbin$Sample,c("S2"="Inclusion","S6"="Day 1 of infection"))
p1_taxbin$Abundance<-plyr::revalue(p1_taxbin$Abundance,c("0"="Absent"))
p1_taxbin$OTU<-plyr::revalue(p1_taxbin$OTU,c("Citrobacater spp"="Citrobacter spp"))
View(p1_taxbin)

Patient1<-ggplot(p1_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 1",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))


Patient1 




######################
######################
#####################
## Patient 2 : Infection with Citrobacter freundii

p2<-subset_samples(physeqLung2,record_ID=="2")
p2<-microbiome::transform(p2, "compositional")
View(psmelt(p2))

p2df<-psmelt(p2)
p2df<-p2df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p2df)


#Plot heatmap

p2_taxbin<-psmelt(tax_transform(p2, trans = "binary", rank = "Species")%>%otu_get())

View(p2_taxbin)
p2_taxbin<-p2_taxbin%>%filter(OTU=="Pseudocitrobacter_sp")

class(p2_taxbin)
class(p2_taxbin$Sample)
class(p2_taxbin$Abundance)
class(p2_taxbin$OTU)
p2_taxbin$OTU<-as.factor(p2_taxbin$OTU)
p2_taxbin$Abundance<-as.factor(p2_taxbin$Abundance)
p2_taxbin$Sample<-as.factor(p2_taxbin$Sample)
levels(p2_taxbin$Sample)
levels(p2_taxbin$Abundance)


p2_taxbin$Sample<-plyr::revalue(p2_taxbin$Sample,c("S49"="Day 1 of infection","S53"="Day 5 of infection", "S57"="Extubation"))
p2_taxbin$Abundance<-plyr::revalue(p2_taxbin$Abundance,c("0"="Absent"))
p2_taxbin$OTU<-plyr::revalue(p2_taxbin$OTU,c("Pseudocitrobacter_sp"="Citrobacter spp"))
View(p2_taxbin)

Patient2<-ggplot(p2_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 2",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))


Patient2 



ggpubr::ggarrange(Patient1, Patient2, heights = c(3, 3), ncol = 1, align = "hv",common.legend = TRUE,legend="right",vjust = "2.5")

######################
######################
#####################


## Patient 4 : Infection with Haemophilus influenzae

p4<-subset_samples(physeqLung2,record_ID=="4")
p4<-microbiome::transform(p4, "compositional")
View(psmelt(p4))

p4df<-psmelt(p4)
p4df<-p4df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p4df)


#Plot heatmap

p4_taxbin<-psmelt(tax_transform(p4, trans = "binary", rank = "Species")%>%otu_get())

View(p4_taxbin)
p4_taxbin<-p4_taxbin%>%filter(OTU=="Haemophilus influenzae"|OTU=="Haemophilus aegyptius")

class(p4_taxbin)
class(p4_taxbin$Sample)
class(p4_taxbin$Abundance)
class(p4_taxbin$OTU)
p4_taxbin$OTU<-as.factor(p4_taxbin$OTU)
p4_taxbin$Abundance<-as.factor(p4_taxbin$Abundance)
p4_taxbin$Sample<-as.factor(p4_taxbin$Sample)
levels(p4_taxbin$Sample)
levels(p4_taxbin$Abundance)


p4_taxbin$Sample<-plyr::revalue(p4_taxbin$Sample,c("S134"="Day 1 of infection","S508"="Extubation"))
p4_taxbin$Abundance<-plyr::revalue(p4_taxbin$Abundance,c("0"="Absent","1"="Present"))
p4_taxbin$OTU<-plyr::revalue(p4_taxbin$OTU,c("Haemophilus aegyptius"="Haemophilus aegyptius"))
View(p4_taxbin)

Patient4<-ggplot(p4_taxbin%>%filter(OTU=="Haemophilus influenzae"), aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 4",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient4



######################
######################
#####################

## Patient 8 : Infection with Aeromonas and E.coli

p8<-subset_samples(physeqLung2,record_ID=="8")
p8<-microbiome::transform(p8, "compositional")
View(psmelt(p8))

p8df<-psmelt(p8)
p8df<-p8df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p8df)


#Plot heatmap

p8_taxbin<-psmelt(tax_transform(p8, trans = "binary", rank = "Species")%>%otu_get())

View(p8_taxbin)
p8_taxbin<-p8_taxbin%>%filter(OTU=="Escherichia_sp"|OTU=="Aeromonas_sp")

class(p8_taxbin)
class(p8_taxbin$Sample)
class(p8_taxbin$Abundance)
class(p8_taxbin$OTU)
p8_taxbin$OTU<-as.factor(p8_taxbin$OTU)
p8_taxbin$Abundance<-as.factor(p8_taxbin$Abundance)
p8_taxbin$Sample<-as.factor(p8_taxbin$Sample)
levels(p8_taxbin$Sample)
levels(p8_taxbin$Abundance)
levels(p8_taxbin$OTU)

p8_taxbin$Sample<-plyr::revalue(p8_taxbin$Sample,c("S801"="Day 1 of infection","S807"="Day 5 of infection","S812"="Extubation"))
p8_taxbin$Abundance<-plyr::revalue(p8_taxbin$Abundance,c("0"="Absent","1"="Present"))
p8_taxbin$OTU<-plyr::revalue(p8_taxbin$OTU,c("Escherichia_sp"="E.coli", "Aeromonas_sp"="Aeromonas spp."))
View(p8_taxbin)

Patient8<-ggplot(p8_taxbin%>%filter(OTU=="Aeromonas spp."), aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 8",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))


Patient8


######################
######################
#####################

## Patient 9 : Infection with Staphylococcus aureus

p9<-subset_samples(physeqLung2,record_ID=="9")
p9<-microbiome::transform(p9, "compositional")
View(psmelt(p9))

p9df<-psmelt(p9)
p9df<-p9df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p9df)


#Plot heatmap

p9_taxbin<-psmelt(tax_transform(p9, trans = "binary", rank = "Species")%>%otu_get())

View(p9_taxbin)
p9_taxbin<-p9_taxbin%>%filter(OTU=="Staphylococcus aureus/schweitzeri/argenteus(3)")

class(p9_taxbin)
class(p9_taxbin$Sample)
class(p9_taxbin$Abundance)
class(p9_taxbin$OTU)
p9_taxbin$OTU<-as.factor(p9_taxbin$OTU)
p9_taxbin$Abundance<-as.factor(p9_taxbin$Abundance)
p9_taxbin$Sample<-as.factor(p9_taxbin$Sample)
levels(p9_taxbin$Sample)
levels(p9_taxbin$Abundance)


p9_taxbin$Sample<-plyr::revalue(p9_taxbin$Sample,c("S905"="Inclusion","S907"="Day 1 of infection"))
p9_taxbin$Abundance<-plyr::revalue(p9_taxbin$Abundance,c("1"="Present"))
p9_taxbin$OTU<-plyr::revalue(p9_taxbin$OTU,c("Staphylococcus aureus/schweitzeri/argenteus(3)"="Staphylococcus aureus"))
View(p9_taxbin)

Patient9<-ggplot(p9_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 9",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient9  

ggsave("patient1_culture_vs_meta.tiff",
       width = 25, height = 15, dpi = 300, units = "cm")


ggpubr::ggarrange(Patient1, Patient9, heights = c(3, 3), ncol = 1, align = "v",common.legend = TRUE,legend="right")


ggsave("culture_vs_meta.tiff",
       width = 25, height = 15, dpi = 300, units = "cm")


######################
######################
#####################

## Patient 11 : Infection with Klebsiella pneumoniae and MOrganella morganii

p11<-subset_samples(physeqLung2,record_ID=="11")
p11<-microbiome::transform(p11, "compositional")
View(psmelt(p11))

p11df<-psmelt(p11)
p11df<-p11df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p11df)


#Plot heatmap

p11_taxbin<-psmelt(tax_transform(p11, trans = "binary", rank = "Species")%>%otu_get())

View(p11_taxbin)
p11_taxbin<-p11_taxbin%>%filter(OTU=="Klebsiella pneumoniae"|OTU=="Morganella morganii")

class(p11_taxbin)
class(p11_taxbin$Sample)
class(p11_taxbin$Abundance)
class(p11_taxbin$OTU)
p11_taxbin$OTU<-as.factor(p11_taxbin$OTU)
p11_taxbin$Abundance<-as.factor(p11_taxbin$Abundance)
p11_taxbin$Sample<-as.factor(p11_taxbin$Sample)
levels(p11_taxbin$Sample)
levels(p11_taxbin$Abundance)
levels(p11_taxbin$OTU)

p11_taxbin$Sample<-plyr::revalue(p11_taxbin$Sample,c("S1107"="Day 1 of infection","S1112"="Extubation"))
p11_taxbin$Abundance<-plyr::revalue(p11_taxbin$Abundance,c("0"="Absent","1"="Present"))
View(p11_taxbin)

Patient11<-ggplot(p11_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 11",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))


Patient11

######################
######################
#####################

## Patient 12 : Infection with Staphylococcus aureus

p12<-subset_samples(physeqLung2,record_ID=="12")
p12<-microbiome::transform(p12, "compositional")
View(psmelt(p12))

p12df<-psmelt(p12)
p12df<-p12df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p12df)


#Plot heatmap

p12_taxbin<-psmelt(tax_transform(p12, trans = "binary", rank = "Species")%>%otu_get())

View(p12_taxbin)
p12_taxbin<-p12_taxbin%>%filter(OTU=="Staphylococcus aureus/schweitzeri/argenteus(3)")

class(p12_taxbin)
class(p12_taxbin$Sample)
class(p12_taxbin$Abundance)
class(p12_taxbin$OTU)
p12_taxbin$OTU<-as.factor(p12_taxbin$OTU)
p12_taxbin$Abundance<-as.factor(p12_taxbin$Abundance)
p12_taxbin$Sample<-as.factor(p12_taxbin$Sample)
levels(p12_taxbin$Sample)
levels(p12_taxbin$Abundance)


p12_taxbin$Sample<-plyr::revalue(p12_taxbin$Sample,c("S1204"="Day 1 of infection"))
p12_taxbin$Abundance<-plyr::revalue(p12_taxbin$Abundance,c("1"="Present"))
p12_taxbin$OTU<-plyr::revalue(p12_taxbin$OTU,c("Staphylococcus aureus/schweitzeri/argenteus(3)"="Staphylococcus aureus"))
View(p12_taxbin)

Patient12<-ggplot(p12_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 12",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient12  


######################
######################
#####################

## Patient 26 : Infection with Serratia marcescens

p26<-subset_samples(physeqLung2,record_ID=="26")
p26<-microbiome::transform(p26, "compositional")
View(psmelt(p26))

p26df<-psmelt(p26)
p26df<-p26df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p26df)


#Plot heatmap

p26_taxbin<-psmelt(tax_transform(p26, trans = "binary", rank = "Species")%>%otu_get())

View(p26_taxbin)
p26_taxbin<-p26_taxbin%>%filter(OTU=="Serratia marcescens(1)")

class(p26_taxbin)
class(p26_taxbin$Sample)
class(p26_taxbin$Abundance)
class(p26_taxbin$OTU)
p26_taxbin$OTU<-as.factor(p26_taxbin$OTU)
p26_taxbin$Abundance<-as.factor(p26_taxbin$Abundance)
p26_taxbin$Sample<-as.factor(p26_taxbin$Sample)
levels(p26_taxbin$Sample)
levels(p26_taxbin$Abundance)


p26_taxbin$Sample<-plyr::revalue(p26_taxbin$Sample,c("S2607"="Day 1 of infection", "S2608"="Extubation"))
p26_taxbin$Abundance<-plyr::revalue(p26_taxbin$Abundance,c("1"="Present"))
p26_taxbin$OTU<-plyr::revalue(p26_taxbin$OTU,c("Serratia marcescens(1)"="Serratia marcescens"))
View(p26_taxbin)

Patient26<-ggplot(p26_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Present"="darkred","Absent"="grey"))+
  labs(title="Patient 26",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient26  

######################
######################
#####################

## Patient 28 : Infection with Haemophilus influenzae

p28<-subset_samples(physeqLung2,record_ID=="28")
p28<-microbiome::transform(p28, "compositional")
View(psmelt(p28))

p28df<-psmelt(p28)
p28df<-p28df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p28df)


#Plot heatmap

p28_taxbin<-psmelt(tax_transform(p28, trans = "binary", rank = "Species")%>%otu_get())

View(p28_taxbin)
p28_taxbin<-p28_taxbin%>%filter(OTU=="Haemophilus influenzae"|OTU=="Haemophilus aegyptius")

class(p28_taxbin)
class(p28_taxbin$Sample)
class(p28_taxbin$Abundance)
class(p28_taxbin$OTU)
p28_taxbin$OTU<-as.factor(p28_taxbin$OTU)
p28_taxbin$Abundance<-as.factor(p28_taxbin$Abundance)
p28_taxbin$Sample<-as.factor(p28_taxbin$Sample)
levels(p28_taxbin$Sample)
levels(p28_taxbin$Abundance)


p28_taxbin$Sample<-plyr::revalue(p28_taxbin$Sample,c("S2802"="Inclusion","S2806"="Day 1 of infection","S2812"="Day 5 of infection","S2815"="Extubation"))
p28_taxbin$Abundance<-plyr::revalue(p28_taxbin$Abundance,c("0"="Absent","1"="Present"))

View(p28_taxbin)
# A that point we are going to drop "Haemophilus aegyptius" level since haemophilus is present to both samples

p28_taxbin<-p28_taxbin%>%filter(OTU=="Haemophilus influenzae")
View(p28_taxbin)

Patient28<-ggplot(p28_taxbin%>%filter(OTU=="Haemophilus influenzae"), aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 28",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient28  



######################
######################
#####################
## Patient 30 : Infection with Haemophilus influenzae and E.coli

p30<-subset_samples(physeqLung2,record_ID=="30")
p30<-microbiome::transform(p30, "compositional")
View(psmelt(p30))

p30df<-psmelt(p30)
p30df<-p30df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p30df)


#Plot heatmap

p30_taxbin<-psmelt(tax_transform(p30, trans = "binary", rank = "Species")%>%otu_get())

View(p30_taxbin)
p30_taxbin<-p30_taxbin%>%filter(OTU=="Haemophilus influenzae"|OTU=="Escherichia_sp")

class(p30_taxbin)
class(p30_taxbin$Sample)
class(p30_taxbin$Abundance)
class(p30_taxbin$OTU)
p30_taxbin$OTU<-as.factor(p30_taxbin$OTU)
p30_taxbin$Abundance<-as.factor(p30_taxbin$Abundance)
p30_taxbin$Sample<-as.factor(p30_taxbin$Sample)
levels(p30_taxbin$Sample)
levels(p30_taxbin$Abundance)


p30_taxbin$Sample<-plyr::revalue(p30_taxbin$Sample,c("S3005"="Day 1 of infection","S3011"="Day 5 of infection","S3013"="Extubation"))
p30_taxbin$Abundance<-plyr::revalue(p30_taxbin$Abundance,c("0"="Absent","1"="Present"))

View(p30_taxbin)

Patient30<-ggplot(p30_taxbin%>%filter(OTU=="Haemophilus influenzae"), aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 30",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient30 


######################
######################
#####################

## Patient 33: Infection with Haemophilus influenzae

p33<-subset_samples(physeqLung2,record_ID=="33")
p33<-microbiome::transform(p33, "compositional")
View(psmelt(p33))

p33df<-psmelt(p33)
p33df<-p33df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p33df)


#Plot heatmap

p33_taxbin<-psmelt(tax_transform(p33, trans = "binary", rank = "Species")%>%otu_get())

View(p33_taxbin)
p33_taxbin<-p33_taxbin%>%filter(OTU=="Haemophilus influenzae"|OTU=="Haemophilus aegyptius")

class(p33_taxbin)
class(p33_taxbin$Sample)
class(p33_taxbin$Abundance)
class(p33_taxbin$OTU)
p33_taxbin$OTU<-as.factor(p33_taxbin$OTU)
p33_taxbin$Abundance<-as.factor(p33_taxbin$Abundance)
p33_taxbin$Sample<-as.factor(p33_taxbin$Sample)
levels(p33_taxbin$Sample)
levels(p33_taxbin$Abundance)


p33_taxbin$Sample<-plyr::revalue(p33_taxbin$Sample,c("S3305"="Inclusion","S3307"="Day 1 of infection"))
p33_taxbin$Abundance<-plyr::revalue(p33_taxbin$Abundance,c("0"="Absent","1"="Present"))

View(p33_taxbin)


Patient33<-ggplot(p33_taxbin%>%filter(OTU=="Haemophilus influenzae"), aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 33",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient33  



######################
######################
#####################


## Patient 36: Infection with Haemophilus influenzae and Streètococcus pneumoniae

p36<-subset_samples(physeqLung2,record_ID=="36")
p36<-microbiome::transform(p36, "compositional")
View(psmelt(p36))

p36df<-psmelt(p36)
p36df<-p36df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p36df)


#Plot heatmap

p36_taxbin<-psmelt(tax_transform(p36, trans = "binary", rank = "Species")%>%otu_get())

View(p36_taxbin)
p36_taxbin<-p36_taxbin%>%filter(OTU=="Haemophilus influenzae"|OTU=="Streptococcus pneumoniae")

class(p36_taxbin)
class(p36_taxbin$Sample)
class(p36_taxbin$Abundance)
class(p36_taxbin$OTU)
p36_taxbin$OTU<-as.factor(p36_taxbin$OTU)
p36_taxbin$Abundance<-as.factor(p36_taxbin$Abundance)
p36_taxbin$Sample<-as.factor(p36_taxbin$Sample)
levels(p36_taxbin$Sample)
levels(p36_taxbin$Abundance)


p36_taxbin$Sample<-plyr::revalue(p36_taxbin$Sample,c("S3605"="Day 1 of infection","S3609"="Day 5 of infection"))
p36_taxbin$Abundance<-plyr::revalue(p36_taxbin$Abundance,c("0"="Absent","1"="Present"))

View(p36_taxbin)


Patient36<-ggplot(p36_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 36",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient36  



######################
######################
#####################

# Patient 40 who had an infection due to Klebsiella pneumoniae

p40<-subset_samples(physeqLung2,record_ID=="41")
p40<-microbiome::transform(p40, "compositional")
View(psmelt(p40))


p40df<-psmelt(p40)
p40df<-p40df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
View(p40df)


#Plot heatmap

p40_taxbin<-psmelt(tax_transform(p40, trans = "binary", rank = "Species")%>%otu_get())
View(p40_taxbin)
p40_taxbin<-p40_taxbin%>%filter(OTU=="Klebsiella pneumoniae")

class(p40_taxbin)
class(p40_taxbin$Sample)
class(p40_taxbin$Abundance)
class(p40_taxbin$OTU)
p40_taxbin$OTU<-as.factor(p40_taxbin$OTU)
p40_taxbin$Abundance<-as.factor(p40_taxbin$Abundance)
p40_taxbin$Sample<-as.factor(p40_taxbin$Sample)
levels(p40_taxbin$Sample)
levels(p40_taxbin$Abundance)


p40_taxbin$Sample<-plyr::revalue(p40_taxbin$Sample,c("S4005"="Inclusion","S4010"="Day 1 of infection"))
p40_taxbin$Abundance<-plyr::revalue(p40_taxbin$Abundance,c("0"="Absent"))

View(p40_taxbin)

Patient40<-ggplot(p40_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 40",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient40  


######################
######################
#####################

# Patient 41 who had an infection due to E.coli and Staphylococcus aureus

p41<-subset_samples(physeqLung2,record_ID=="41")
p41<-microbiome::transform(p41, "compositional")
View(psmelt(p41))

#Plot heatmap

p41_taxbin<-psmelt(tax_transform(p41, trans = "binary", rank = "Species")%>%otu_get())
View(p41_taxbin)
p41_taxbin<-p41_taxbin%>%filter(OTU=="Staphylococcus aureus/schweitzeri/argenteus(3)"|OTU=="Escherichia_sp")

class(p41_taxbin)
class(p41_taxbin$Sample)
class(p41_taxbin$Abundance)
class(p41_taxbin$OTU)
p41_taxbin$OTU<-as.factor(p41_taxbin$OTU)
p41_taxbin$Abundance<-as.factor(p41_taxbin$Abundance)
p41_taxbin$Sample<-as.factor(p41_taxbin$Sample)
levels(p41_taxbin$Sample)
levels(p41_taxbin$Abundance)


p41_taxbin$Sample<-plyr::revalue(p41_taxbin$Sample,c("S4104"="Inclusion","S4110"="Day 1 of infection","S4112"="Day 5 of infection"))
p41_taxbin$Abundance<-plyr::revalue(p41_taxbin$Abundance,c("0"="Absent","1"="Present"))

View(p41_taxbin)

Patient41<-ggplot(p41_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 41",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient41  

######################
######################
#####################

# Patient 42 who had an infection due to Haemophilus influenzae and Klebsiella pneumoniae

p42<-subset_samples(physeqLung2,record_ID=="42")
p42<-microbiome::transform(p42, "compositional")
View(psmelt(p42))

#Plot heatmap

p42_taxbin<-psmelt(tax_transform(p42, trans = "binary", rank = "Species")%>%otu_get())
View(p42_taxbin)
p42_taxbin<-p42_taxbin%>%filter(OTU=="Haemophilus influenzae"|OTU=="Klebsiella pneumoniae")

class(p42_taxbin)
class(p42_taxbin$Sample)
class(p42_taxbin$Abundance)
class(p42_taxbin$OTU)
p42_taxbin$OTU<-as.factor(p42_taxbin$OTU)
p42_taxbin$Abundance<-as.factor(p42_taxbin$Abundance)
p42_taxbin$Sample<-as.factor(p42_taxbin$Sample)
levels(p42_taxbin$Sample)
levels(p42_taxbin$Abundance)


p42_taxbin$Sample<-plyr::revalue(p42_taxbin$Sample,c("S4204"="Inclusion"))
p42_taxbin$Abundance<-plyr::revalue(p42_taxbin$Abundance,c("0"="Absent","1"="Present"))

View(p42_taxbin)

Patient42<-ggplot(p42_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 42",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient42  


######################
######################
#####################


# Patient 43 who had an infection due to Klebsiella pneumoniae and Morganella morganii

p43<-subset_samples(physeqLung2,record_ID=="43")
p43<-microbiome::transform(p43, "compositional")
View(psmelt(p43))

p43df<-psmelt(p43)
p43df<-p43df%>%select(OTU,Abundance,Filtered_Data_AK_Dec2020,Species)
p43df

#Plot heatmap

p43_taxbin<-(psmelt(tax_transform(p43, trans = "binary", rank = "Species")%>%otu_get()))
View(p43_taxbin)
p43_taxbin<-p43_taxbin%>%filter(OTU=="Klebsiella pneumoniae"|OTU=="Morganella morganii")

class(p43_taxbin)
class(p43_taxbin$Sample)
class(p43_taxbin$Abundance)
class(p43_taxbin$OTU)
p43_taxbin$OTU<-as.factor(p43_taxbin$OTU)
p43_taxbin$Abundance<-as.factor(p43_taxbin$Abundance)
p43_taxbin$Sample<-as.factor(p43_taxbin$Sample)
levels(p43_taxbin$Sample)
levels(p43_taxbin$Abundance)


p43_taxbin$Sample<-plyr::revalue(p43_taxbin$Sample,c("S4305"="Day 1 of infection","S4306"="Day 5 of infection"))
p43_taxbin$Abundance<-plyr::revalue(p43_taxbin$Abundance,c("1"="Present", "0"="Absent"))
View(p43_taxbin)

Patient43<-ggplot(p43_taxbin, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="Patient 43",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))

Patient43  



######################
######################
#####################


## Create a new dataframe containing all previously create dataframes with presence or absence of bacteria of interest to
## create the presence/absence heatmap

View(p1_taxbin)
View(p2_taxbin)

p1_taxbin$id<-sprintf("1")
p2_taxbin$id<-sprintf("2")
p4_taxbin$id<-sprintf("4")
p8_taxbin$id<-sprintf("8")
p9_taxbin$id<-sprintf("9")
p11_taxbin$id<-sprintf("11")
p12_taxbin$id<-sprintf("12")
p26_taxbin$id<-sprintf("26")
p28_taxbin$id<-sprintf("28")
p30_taxbin$id<-sprintf("30")
p33_taxbin$id<-sprintf("33")
p36_taxbin$id<-sprintf("36")
p40_taxbin$id<-sprintf("40")
p41_taxbin$id<-sprintf("41")
p42_taxbin$id<-sprintf("42")
p43_taxbin$id<-sprintf("43")

## Filter to exclude OTUs of 2nd microbiological sample that were not taken on the same time as metagenomic data
p4_taxbin<-p4_taxbin%>%filter(!OTU=="Haemophilus aegyptius")
p8_taxbin<-p8_taxbin%>%filter(!OTU=="E.coli")
p28_taxbin<-p28_taxbin%>%filter(!OTU=="Haemophilus aegyptius")
p30_taxbin<-p30_taxbin%>%filter(!OTU=="Escherichia_sp")
p33_taxbin<-p33_taxbin%>%filter(!OTU=="Haemophilus aegyptius")


p33_taxbin

all<-rbind(p1_taxbin,p2_taxbin,p4_taxbin,p8_taxbin,p9_taxbin,p11_taxbin,p12_taxbin,p26_taxbin,
           p28_taxbin,p30_taxbin,p33_taxbin,p36_taxbin,p40_taxbin,p41_taxbin,p42_taxbin,p43_taxbin)

View(all)

all$OTU<-plyr::revalue(all$OTU,c("Staphylococcus aureus/schweitzeri/argenteus(3)"="Staphylococcus aureus","Escherichia_sp"="E.coli"))
all$OTU<-plyr::revalue(all$OTU,c("E.coli"="Escherichia coli"))

all_plot<-ggplot(all, aes(Sample, OTU,fill=Abundance)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1)+
  coord_fixed()+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("Absent"="grey","Present"="darkred"))+
  labs(title="",
       subtitle="",
       fill="Result")+
  theme(axis.title.y= element_text(size = rel(1), angle = 0,face="bold"),axis.text = element_text(size=10),legend.text = element_text(size=10),
        legend.title = element_text(size=rel(1)),
        axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90),
        plot.title = element_text(face = "bold",size=12))+
  facet_wrap(~factor(id,levels=c(1,2,4,8,9,11,12,26,28,30,33,36,40,41,42,43)),ncol=8)

all_plot
all_plot<-all_plot+theme_bw()

all_plot+
  theme(axis.text.y =element_text(face="bold.italic"),
        axis.text.x =element_text(size=10,angle=90,face="bold"),
        legend.title = element_text(size=rel(1.8)))

ggsave("Presence/absence_heatmap.tiff",
       width = 6, height = 5, dpi = 300,units = "cm")


######################
######################
#####################

## Now we are going to filter the phyloseq object for the patients that interests us and the above mentioned germs to produce the same type of 
# graph but using a colored heatmap with annotation of the values of relative abundance in the squares.

pneum<-subset_samples(physeqLung2,record_ID=="1"|record_ID=="2"|record_ID=="4"|record_ID=="8"|record_ID=="9"
                      |record_ID=="11"|record_ID=="12"|record_ID=="26"|record_ID=="28"|record_ID=="30"
                      |record_ID=="33"|record_ID=="36"|record_ID=="40"|record_ID=="41"|record_ID=="42"|record_ID=="43")

pneum<-microbiome::transform(pneum, "compositional")
View(psmelt(pneum))

pneumdf<-psmelt(pneum)
levels(pneumdf$record_ID)

View(pneumdf)

pneumdf<-pneumdf%>%filter(record_ID=="1"&Species=="Pseudocitrobacter_sp"|record_ID=="2"&Species=="Pseudocitrobacter_sp"
                      |record_ID=="4"&Species=="Haemophilus influenzae"
                      |record_ID=="8"&Species=="Aeromonas_sp"
                      |record_ID=="9"&Species=="Staphylococcus aureus/schweitzeri/argenteus(3)"
                      |record_ID=="11"&Species=="Klebsiella pneumoniae"|record_ID=="11"&Species=="Morganella morganii"
                      |record_ID=="12"&Species=="Staphylococcus aureus/schweitzeri/argenteus(3)"
                      |record_ID=="26"&Species=="Serratia marcescens(1)"
                      |record_ID=="28"&Species=="Haemophilus influenzae"
                      |record_ID=="30"&Species=="Haemophilus influenzae"
                      |record_ID=="33"&Species=="Haemophilus influenzae"
                      |record_ID=="36"&Species=="Haemophilus influenzae"|record_ID=="36"&Species=="Streptococcus pneumoniae"
                      |record_ID=="40"&Species=="Klebsiella pneumoniae"
                      |record_ID=="41"&Species=="Staphylococcus aureus/schweitzeri/argenteus(3)"|record_ID=="41"&Species=="Escherichia_sp"
                      |record_ID=="42"&Species=="Klebsiella pneumoniae"|record_ID=="42"&Species=="Haemophilus influenzae"
                      |record_ID=="43"&Species=="Klebsiella pneumoniae"|record_ID=="43"&Species=="Morganella morganii")



pneumdf$Species<-plyr::revalue(pneumdf$Species,c("Staphylococcus aureus/schweitzeri/argenteus(3)"="Staphylococcus aureus","Escherichia_sp"="E.coli",
                                                 "Serratia marcescens(1)"="Serratia marcescens",
                                                 "Pseudocitrobacter_sp"="Citrobacter spp.",
                                                 "Aeromonas_sp"="Aeromonas spp."))


heat<-ggplot(pneumdf, aes(x = Filtered_Data_AK_Dec2020, y = Species, fill = Abundance)) +
  geom_tile(color="white",lwd=1.5,linetype=1)+
  coord_fixed()+
  facet_wrap(~factor(record_ID,levels=c(1,2,4,8,9,11,12,26,28,30,33,36,40,41,42,43)),ncol=8)+
  geom_text(aes(label = Abundance), color = "white", size = 4)+
  scale_fill_gradient()

heat

pneumdf$Abundance<-round(pneumdf$Abundance,digits=0)
pneumdf$Abundance<-pneumdf$Abundance*100

View(pneumdf)

heat+theme(axis.text.y =element_text(face="bold.italic"),
            axis.text.x =element_text(size=10,angle=90,face="bold"),
            legend.title = element_text(size=rel(1.8)))

  

######################
######################
#####################

######################
######################
#####################






############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 
############################## SEEN HERE LAST TIME  :-) 






















# If present remove any samples with a total number of reads below 1000
sort(sample_sums(physeqLung2), decreasing = T)
sample_sums(physeqLung2) #In our case since alla data have been previously filtered there are no samples with less than 10.000 reads 


View(otu_table(physeqLung2))
tax_table(physeqLung2)

phy_tree(physeqLung2)

taxa_names(physeqLung2)
rank_names(physeqLung2)


#PLot phylogenetic tree of all samples
plot_tree(physeqLung2, "treeonly", nodeplotblank, ladderize="left",label.tips="taxa_names")
plot_tree(physeqLung2, "anythingelse")
plot_tree(physeqLung2, nodelabf=nodeplotboot(), ladderize="left", color="Sample")
plot_tree(physeqLung2, nodelabf=nodeplotboot(), ladderize="left", color="Class")

#Tree to Phylum level
plot_tree(physeqLung2, nodelabf=nodeplotboot(), ladderize="left", color="Phylum",base.spacing=0.01,plot.margin=0.01)+
  labs(title="Phylogenetic tree (all samples)-Phylum level",subtitle="Rooted tree, 2188 tips and 2027 internal nodes",caption="*(Presence of a long branched root belonging to either Synergistetes or Tenericutes phylum)")

ggsave("Phylum_tree_all_samples.pdf",
       width = 40, height = 30, dpi = 300, units = "cm")

#Tree to class level
plot_tree(physeqLung2, nodelabf=nodeplotboot(), ladderize="left", color="Class",base.spacing=0.01,plot.margin=0.01)+
  labs(title="Phylogenetic tree (all samples)-Class level",subtitle="Rooted tree, 2188 tips and 2027 internal nodes")

ggsave("Class_tree_all_samples.pdf",
       width = 40, height = 30, dpi = 300, units = "cm")


#Plot tree of the 20 most abundant taxa
myTaxa = names(sort(taxa_sums(physeqLung2), decreasing = TRUE)[1:20])

ex1 = prune_taxa(myTaxa, physeqLung2)

plot(phy_tree(ex1), show.node.label = TRUE)

phy_tree(ex1)

#Stratification based on timepoint
plot_tree(ex1, color = "Subgroup",label.tips="Family",ladderize="left",base.spacing=0.03,plot.margin= 0.1,size="Abundance",sizebase=10,nodelabf=nodeplotboot(),text.size=4,
          title="Phylogenetic tree of the 20 most abundant taxa (stratified by timepoint)")+
  scale_color_manual(values = InclInfectColors3)

ggsave("Tree_20mostabund_based_on_timepoint.pdf",
       width = 35, height = 20, dpi = 300, units = "cm")


#Stratiication based on infection group
plot_tree(ex1,size="abundance",color="Group",base.spacing=0.02,shape="Subgroup",label.tips="Species",plot.margin=0.1,text.size=4,sizebase=10,
          title="Phylogenetic tree of the 20 most abundant taxa (stratified by infection group)")+
  scale_color_manual(values = Infectcolorgroup)+
  coord_polar(theta="y")

ggsave("Tree_20mostabund_based_on_infection_group.pdf",
       width = 35, height = 20, dpi = 300, units = "cm")


## Keep only taxa with total abundance at least 2 in at least 1 samples
# keep_function <- function(x) {x >= 2}
# TaxaToKeep <- genefilter_sample(physeqLung2NoExclNoLo, keep_function, A = 1)

# TaxaToKeep
# TaxaToKeep[1:20]

# physeqLung4 <- prune_taxa(TaxaToKeep, physeqLung2NoExclNoLo)

# sort(taxa_sums(physeqLung4), decreasing = F)[1:50]


################################################
###############################################
################################################

### LUNG DATA ANALYSIS   ####

## 1st step: describe the taxonomic composition at inclusion of patients to see if there are any differences in the 3 groups of patients

## Plot taxa prevalence according to abundance of features. This will be helpful to apply filtering based on downstream analyses.

## Create Taxa Prevalence plots at inclusion

#Keep only inclusion samples
ps1<- subset_samples(physeqLung2,Filtered_Data_AK_Dec2020=="Inclusion")

View(psmelt(ps1))
ps2<-psmelt(ps1) # transform phyloseq object to dataframe
class(ps2)


# compute prevalence of each feature,store as data.frame
prevdf=apply(X=otu_table(ps1),
             MARGIN = ifelse(taxa_are_rows(ps1),yes=1,no=2),
             FUN=function(x){sum(x>0)})

class(prevdf)
View(prevdf)

# Add taxonomy and total read counts to this data.frame
prevdf=data.frame(Prevalence=prevdf,TotalAbundance=taxa_sums(ps1),tax_table(ps1),otu_table(ps1))
View(prevdf)


# Examine if any phyla only have low abundance features
plyr::ddply(prevdf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})


# Now we have our 2 dataframes: the first coming from a phyloseq object containing metadata that was trnasformed to df
# and the second the dataframe of prevalence

View(ps2)
View(prevdf)

# Ad rownames as column to prevdf dataframe
prevdf$OTU<-rownames(prevdf)
View(prevdf)
colnames(prevdf)
prevdf<-prevdf%>%relocate(OTU,.before=Prevalence) # relocate to be the first column
View(prevdf)

#Select only columns of intester. Taxonomix features are eliminated because they a re already contained in the ps2 dataframe
prevdf2<-prevdf%>%select(OTU,Prevalence,TotalAbundance,S2,S88,S702,S905,S1605,S1805,S1902,S2005,S2105,S2402,S2802,S3205,S3305,S3905,S4005,S4104,S4204)

View(prevdf2)

#Transform form wide to long format
prevdf2<-melt(prevdf2, id.vars=c("OTU", "Prevalence","TotalAbundance"))
View(prevdf2)

prevdf2<-prevdf2%>%relocate(variable,.after=OTU) #relocate "variable" column

prevdf2<-dplyr::rename(prevdf2,Sample=variable)

#create a new column to stock a new OTU_Sample code that will serve as to merge the 2 dataframes later

prevdf2$OTU_Sample<-sprintf("")
prevdf2<-prevdf2%>%relocate(OTU_Sample,.after=Sample) #relocate OTU_Sample column
View(prevdf2)

prevdf2$OTU_Sample<-str_c(prevdf2$OTU,"_",prevdf2$Sample) #fill-in the OTU_Sample column by concatenating the two other columns (OTU and Sample)


# Create the same column OTU_Sample in the ps2 dataframe to be able to merge the data later on

ps2$OTU_Sample<-sprintf("")
ps2<-ps2%>%relocate(OTU_Sample,.after=Sample)
View(ps2)
ps2$OTU_Sample<-str_c(ps2$OTU,"_",ps2$Sample)


# Now it is  time to merge the two dataframes by the "OTU__Sample" column
# But first subset once again the prevdf2 dataframe to keep only prevalence and abundance columns
prevdf2<-prevdf2%>%select(OTU_Sample,Prevalence,TotalAbundance)
View(prevdf2)


# Join the two dataframes using the left_join function
prev<-left_join(ps2,prevdf2,by="OTU_Sample")
prev<-prev%>%relocate(Prevalence,.after=OTU_Sample)
prev<-prev%>%relocate(TotalAbundance,.after=Abundance)
View(prev)


#Plot taxa prevalence based on abundance using ggplot

ggplot(prev,aes(TotalAbundance,Prevalence/nsamples(ps1)*100,color=Phylum))+
  geom_point(size=2,alpha=0.7)+
  scale_x_log10()+
  xlab("Total Abundance (log10)")+
  ylab("Prevalence (%)")+
  facet_wrap(~Phylum)+
  labs(title="Taxa prevalence plot at inclusion")+
  theme(legend.position="none")

ggsave("Taxa_prevalence_plot_phylum_level_inclusion.tiff",
       width = 25, height = 20, dpi = 300, units = "cm")


ggplot(prev,aes(TotalAbundance,Prevalence/nsamples(ps1)*100,color=Group))+
  geom_point(size=2,alpha=0.7)+
  geom_hline(yintercept=5,alpha=0.5,linetype=2)+
  scale_x_log10()+
  xlab("Total Abundance (log10)")+
  ylab("Prevalence (%)")+
  facet_wrap(~Phylum)+
  labs(title="Taxa prevalence plot at inclusion")

ggsave("Taxa_prevalence_plot_phylum_level_inclusion_group.tiff",
       width = 25, height = 20, dpi = 300, units = "cm")



## Create boxplots of taxa stratified by group pf patients using the "ggstripchart" function

# Phylm level abundance per group : all timepoints/all samples

ps_df<-microbiomeutilities::phy_to_ldf(physeqLung2, 
                                transform.counts = "compositional")

View(ps_df)

ggstripchart(ps_df,"Group","Abundance",facet.by="Phylum",color="Group",palette="jco")+
  stat_compare_means(label.y=0.8,label.x=1.5,hide.ns=TRUE)+
  ylab("Relative abundance")+
  xlab("")+
  theme(legend.position="bottom")+
  labs(title="Phyla relative abundance per group (all samples)")

ggsave("Phyla_abundance_per_group.tiff",
       width = 30, height = 25, dpi = 300, units = "cm")




# Phylum level abundance per group: Only inclusion timepoint

ps1<- subset_samples(physeqLung2,Filtered_Data_AK_Dec2020=="Inclusion")

ps_df2<-microbiomeutilities::phy_to_ldf(ps1, 
                                       transform.counts = "compositional")

View(ps_df2)

ggstripchart(ps_df2,"Group","Abundance",facet.by="Phylum",color="Group",palette="jco")+
  stat_compare_means(label.y=0.6,label.x=1.7,hide.ns=TRUE)+
  ylab("Relative abundance")+
  xlab("")+
  theme(legend.position="bottom")+
  labs(title="Phyla relative abundance per group (Inclusion)")

ggsave("Phyla_abundance_per_group_inclusion.tiff",
       width = 30, height = 25, dpi = 300, units = "cm")



###############################
###############################



# Keep only "Inclusion" samples to see taxonomical composition on admission 

ps1<- subset_samples(physeqLung2,Filtered_Data_AK_Dec2020=="Inclusion")


# Create barplots of taxa representation on inclusion across the 3 different groups

#Get the 30 most abundant OTUs / ASVs
ps_tmp<-get_top_taxa(physeq_obj=ps1, n = 30, relative = TRUE,
                     discard_other = FALSE, other_label = "Other")


#Create labels for missing taxonomic ranks
ps_tmp <- name_taxa(ps_tmp, label = "Unkown", species= T, other_label = "Other")

#Generate a barplot that is colored by Phylum and labeled by Species, coloring
#collapsed taxa grey.
fantaxtic_bar(ps_tmp, color_by ="Phylum",label_by="Species", other_label = "Other")

#Generate a barplot that is colored by Phylum and labeled by Species. 
#As multiple ASVs have the same family annotation, generate unique labels.
# then facet by Group of patients
fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Species", other_label = "Other",
              gen_uniq_lbls = TRUE,grid_by = "Group",facet_cols = 3,order_alg="other.abnd")+
  labs(title="Abundance barplots of the 30 most abundant taxa on inclusion")

ggsave("Barplots_30mostabundanttaxa.tiff",
       width = 30, height = 25, dpi = 300, units = "cm")



## Create barplots of the 20 most abundant families on inclusion

ps_tmp_fam<-microbiomeutilities::aggregate_top_taxa2(ps1, top=20, "Family")

ps_tmp_fam

ps_tmp_fam <- name_taxa(ps_tmp_fam, label = "Unkown", species= T, other_label = "Other")


#Generate a barplot that is colored by Phylum and labeled by Species, coloring
#collapsed taxa grey.
fantaxtic_bar(ps_tmp_fam, color_by ="Family",label_by="Family", other_label = "Other")


#Generate a barplot that is colored by Phylum and labeled by Species. As multiple ASVs have the same family annotation, generate unique labels.
fantaxtic_bar(ps_tmp_fam, color_by = "Family", label_by = "Family", other_label = "Other",
              gen_uniq_lbls = TRUE,grid_by = "Group",facet_cols = 3,order_alg="other.abnd")+
  labs(title="Abundance barplots of the 20 most abundant Families on inclusion")

ggsave("Barplots_20mostabundantFam.tiff",
       width = 30, height = 25, dpi = 300, units = "cm")





### Barplots generated by another script
ps3=tax_glom(physeqLung2,"Genus",NArm=TRUE)
ps3

ps3ra=transform_sample_counts(ps3,function(x){x/sum(x)})

ps4<-psmelt(ps3ra)

ggplot(ps4,aes(x=Sample,y=Abundance,fill=Class,order=as.factor(Phylum)))+
  geom_bar(stat="identity",color="black")+
  ggtitle("All Samples-Relative abundance, sorted by Class abundance")+
  theme(axis.text.x=element_blank())

ggsave("Barplots_class_re_abnd.tiff",
       width = 30, height = 25, dpi = 300, units = "cm")



################################################
################################################


### Define the Core Microbiome  ###


microbiome::transform(ps1, "compositional")

head(prevalence(ps1, detection = 0, sort = TRUE))

View(prevalence(ps1, detection = 0, sort = TRUE))

head(prevalence(ps1, detection = 0, sort = TRUE, count = TRUE))


core.taxa.standard <- core_members(ps1, detection = 0, prevalence = 3/100)

core.taxa.standard

pseq.core <- core(physeqLung2, detection = 0, prevalence = .2)

pseq.core2 <- aggregate_rare(physeqLung2, "Genus", detection = 0, prevalence = .2)


core.taxa <- taxa(pseq.core)
core.taxa


core.abundance <- sample_sums(core(physeqLung2, detection = .01, prevalence = .2))

core.abundance

det <- c(0, 0.1, 0.5, 2, 5, 20,100)/100
prevalences <- seq(.1, 1, .1)
plot_core(core.taxa.standard, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")




plot_core(transform(physeqLung2, "compositional"),
               prevalences=seq(0.1, 1, .1), detections=det)

library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
p <- plot_core(pseq.rel,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(pseq.rel, sort = TRUE)[100]) +
  labs(x = "Detection Threshold/n(Relative Abundance (%))") +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)



################################################
################################################

################################################
###############################################
################################################




## Gut data preprocessing ----
OTU_Gut1 <- otu_table(physeqGut1)
sam_Gut1 <- sample_data(physeqGut1)
tax_Gut1 <- tax_table(physeqGut1)
tree_Gut1 <- phy_tree(physeqGut1)
refseq_Gut1 <- refseq(physeqGut1)

# Use phyloseq::sample_data function in view of preparing a phyloseq object
# with the new file of secondary data
sam_Gut2 <- sample_data(SecDataGutFilConv)

# For consistency with the other elements of phyloseq object, use "sampleID"
# as sample_names
sample_names(sam_Gut2) <- as.factor(sam_Gut2$SampleID)
head(sample_names(sam_Gut2))


# Create new phyloseq object with sam_Gut2
physeqGut2 <- phyloseq(OTU_Gut1, tax_Gut1, sam_Gut2, refseq_Gut1, tree_Gut1)
head(sample_names(physeqGut2))
sample_variables(physeqGut2)

# Correct the order of levels for different variables
levels(get_variable(physeqGut2, "Filtered_Data_AK_Dec2020"))

Correct_order_Incl_Inf <- c("Inclusion", "Infect1D1", "Infect1D5",
                            "Extub", "Discharge", "Excluded", "NegCtrl")

# Edit names of the levels
Correct_labels <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation",
                    "Discharge", "Excluded", "Neg control")

sample_data(physeqGut2)$Filtered_Data_AK_Dec2020 <- factor(sample_data(physeqGut2)$Filtered_Data_AK_Dec2020,
                                                           levels = Correct_order_Incl_Inf,
                                                           labels = Correct_labels)
levels(get_variable(physeqGut2, "Filtered_Data_AK_Dec2020"))

# Set a manual scale of colors
InclInfectColors5 <- c("Inclusion"="gold2", "Infection_D1"="coral2",
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3", "Discharge"="orchid1")

# Remove samples with Inclusion-Infection group = "Excluded"
physeqGut2NoExcl <- subset_samples(physeqGut2,
                                   !Filtered_Data_AK_Dec2020 %in% c("Excluded"))

sample_names(physeqGut2NoExcl)

# If present remove any samples with a total number of reads below 1000
physeqGut3 <- prune_samples(sample_sums(physeqGut2NoExcl) > 1000, physeqGut2NoExcl)

sort(sample_sums(physeqGut3), decreasing = F)


## Keep only taxa with total abundance at least 2 in at least 1 samples
keep_function <- function(x) {x >= 2}
TaxaToKeep <- genefilter_sample(physeqGut3, keep_function, A = 1)
TaxaToKeep[1:20]

physeqGut4 <- prune_taxa(TaxaToKeep, physeqGut3)

sort(taxa_sums(physeqGut4), decreasing = F)[1:50]



############################################################################################################## 

############################################################################################################# 

##THE ABOVE STEPS HAVE BEEN PERFORMED, MEANING DATA IMPORT AND PREPROCESSING HAS ALREADY BEEN PERFORMED /# THE FOLLOWING STEPS HAVE TO BEEN PERFORMED. THIS ANALYSIS IS FOCUSING ON GUT MICROBIOTA SO GO DIRECTLY TO LINE 4333 #########
### Lung data analysis / Community abundance manipulation and exploration ----


get_taxa_unique(physeqLung2, "Kingdom")
get_taxa_unique(physeqLung2, "Phylum")
get_taxa_unique(physeqLung2, "Order")
get_taxa_unique(physeqLung2, "Family")
get_taxa_unique(physeqLung2, "Genus")
get_taxa_unique(physeqLung2, "Species")

# Top100 most abundant ASVs
Top100AbundantASVCountsLung <- sort(taxa_sums(physeqLung2), decreasing=TRUE)[1:100]
Top100AbundantASVNamesLung  <-  names(Top100AbundantASVCountsLung)
physeqLung2Top100  <-  prune_taxa(Top100AbundantASVNamesLung, physeqLung2)
physeqLung2Top100

# Visualize taxonomy of most abundant ASVs
Top100AbundantTaxaLung <- tax_table(physeqLung2Top100)[ , c("Phylum", "Genus" ,"Species")]
Top100AbundantTaxaLung

# Number of total ASVs per phylum 
TablePhylaTaxa <- table(tax_table(physeqLung4)[ , 2])
sort(TablePhylaTaxa, decreasing=TRUE)

# Set a palette for Phyla with colours visible to colour-blind people
display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 12, name = "Paired")
PhylumPalLung13 <- c("Firmicutes"="#33A02C", "Fusobacteria"="#CAB2D6",
                     "Actinobacteria"="#E31A1C", "Bacteroidetes"="#FF7F00",
                     "Proteobacteria"="#1F78B4", "Saccharibacteria_TM7"="#FB9A99",
                     "Spirochaetes"="#FDBF6F", "Synergistetes"="#A6CEE3",
                     "Tenericutes"="#B2DF8A", "SR1"="#FF7F00", "Acidobacteria"="#FFFF99",
                     "Chloroflexi"="#6A3D9A", "Deinococcus-Thermus"="#666666")

# Examining the distribution of poorly represented phyla
Acidobacteria_ps <- subset_taxa(physeqLung4, Phylum == "Acidobacteria")
tax_table(Acidobacteria_ps)
taxa_sums(Acidobacteria_ps)
otu_table(Acidobacteria_ps)
# ASV_2170 "Stenotrophobacter" is detected in only 1 sample in low abundance.
# This genus has never been described in human airways, and is associated with
# a single sample with a low 16S copy number. Therefore, this taxon will be removed.
physeqLung5 <- prune_taxa(taxa_names(physeqLung4) != "ASV_2170", physeqLung4)

Chloroflexi_ps <- subset_taxa(physeqLung5, Phylum == "Chloroflexi")
tax_table(Chloroflexi_ps)
taxa_sums(Chloroflexi_ps)
otu_table(Chloroflexi_ps)
# ASV_2144 "Anaerolineae" is detected in only 1 sample in low abundance.
# This class is rarely reported in human airways, and is associated with
# a single sample with a very low 16S copy number. Therefore, this taxon will be removed.
physeqLung6 <- prune_taxa(taxa_names(physeqLung5) != "ASV_2144", physeqLung5)

DeinococcusThermus_ps <- subset_taxa(physeqLung6, Phylum == "Deinococcus-Thermus")
tax_table(DeinococcusThermus_ps)
taxa_sums(DeinococcusThermus_ps)
otu_table(DeinococcusThermus_ps)
# ASV_2056 "Deinococcus geothermalis" is detected in only 1 sample in low abundance.
# This genus is rarely reported in human airways, and is associated with
# a single sample with a low 16S copy number. Therefore, this taxon will be removed.
physeqLung7 <- prune_taxa(taxa_names(physeqLung6) != "ASV_2056", physeqLung6)

SR1_ps <- subset_taxa(physeqLung7, Phylum == "SR1")
tax_table(SR1_ps)
taxa_sums(SR1_ps)
otu_table(SR1_ps)
# ASV_1652 and ASV_1981 were each detected in only 1 sample in low abundance.
# Phylum SR1 is rarely reported in the human airways (oral cavity).
# Therefore, these taxa will be removed.
physeqLung8 <- prune_taxa(taxa_names(physeqLung7) != "ASV_1652", physeqLung7)
physeqLung9 <- prune_taxa(taxa_names(physeqLung8) != "ASV_1981", physeqLung8)

get_taxa_unique(physeqLung9, "Kingdom")
get_taxa_unique(physeqLung9, "Phylum")
get_taxa_unique(physeqLung9, "Order")
get_taxa_unique(physeqLung9, "Family")
get_taxa_unique(physeqLung9, "Genus")
get_taxa_unique(physeqLung9, "Species")

# Transform counts to relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqLung9Rel <- transform_sample_counts(physeqLung2, count_to_rel_abund)

# Check total = 1 for the first 10 samples
sample_sums(physeqLung9Rel)[1:10]



### Composition analysis - Phylum-level ----


physeqLung9RelPhy <- tax_glom(physeqLung9Rel, "Phylum", NArm = TRUE) # Agglomeration

physeqLung9RelPhydf <- psmelt(physeqLung9RelPhy) # Obtain a data frame from phyloseq object.
str(physeqLung9RelPhydf)

physeqLung9RelPhydf %<>% 
  select(OTU, Sample, Abundance, Phylum, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(OTU, Sample, Phylum, Filtered_Data_AK_Dec2020)) %>% arrange(-(Abundance))

PhyDescBoxPlot <- physeqLung9RelPhydf %>% 
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill = Phylum), alpha = 0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Phylum" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Phylum_TotalSamples_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InclInfPhyDescBoxPlot <- PhyDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("Phylum_TotalSamples_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 8, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
physeqLung9RelPhydf$Filtered_Data_AK_Dec2020 <- factor(physeqLung9RelPhydf$Filtered_Data_AK_Dec2020,
                                                       levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

physeqLung9RelPhydf %>%
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("Phylum_TotalSamples_InclInfGroup_Lung_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
# Construct a highly agglomerated tree to obtain a reference for further comparisons
tax_table(physeqLung9RelPhy)

# Merge from Bacteroidetes to Firmicutes to get the full row
# of samples with total counts of Bacteria
MergedPhy <- merge_taxa(physeqLung9RelPhy, taxa_names(physeqLung9RelPhy)[c(1:8)], 1)

plotTreeMerged <- plot_tree(MergedPhy, label.tips="Kingdom", color="Filtered_Data_AK_Dec2020", 
                            size="Abundance", sizebase=10, plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreePhyMergedInclInf_Lung.pdf",
       width = 75, height = 15, dpi = 200, units = "cm", limitsize = FALSE)

# Phylum-level tree
plotTreePhy <- plot_tree(physeqLung9RelPhy,label.tips="Phylum",
                         color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10,
                         plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreePhyInclInf_Lung.pdf",
       width = 55, height = 15, dpi = 200, units = "cm", limitsize = FALSE)

plotTreePhyCirc <- plot_tree(physeqLung9RelPhy, label.tips="Phylum", 
                             color="Filtered_Data_AK_Dec2020", sizebase=10, ladderize="left", 
                             plot.margin= 0.9) + coord_polar(theta="y")

```

### Composition analysis - Class-level ----

```{r}
physeqLung9RelCla <- tax_glom(physeqLung9Rel, "Class", NArm = TRUE) # Agglomeration
ntaxa(physeqLung9RelCla)

physeqLung9RelCladf <- psmelt(physeqLung9RelCla) # Obtain a data frame from phyloseq object.
str(physeqLung9RelCladf)

physeqLung9RelCladf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class) %>% 
  convert(fct(OTU, Sample, Phylum, Class)) %>% arrange(-(Abundance))

ClaDescBoxPlot <- physeqLung9RelCladf %>% 
  ggplot(aes(x = reorder(Class, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Class", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Class_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
plotTreeCla <- plot_tree(physeqLung9RelCla, label.tips="Class",
                         color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10,
                         plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeClaInclInf_Lung.pdf",
       width = 55, height = 15, dpi = 200, units = "cm", limitsize = FALSE)

plotTreeClaCirc <- plot_tree(physeqLung9RelCla, label.tips="Class", 
                             color="Filtered_Data_AK_Dec2020",
                             sizebase=10, ladderize="left", 
                             plot.margin= 0.1) +
  scale_color_manual(values = InclInfectColors) +
  coord_polar(theta="y")

ggsave("PtreeClaInclInf_Lung_Circ.pdf",
       width = 25, height = 25, dpi = 200, units = "cm", limitsize = FALSE)

```

### Composition analysis - Order-level ----

```{r}
physeqLung9RelOrd <- tax_glom(physeqLung9Rel, "Order", NArm = TRUE) # Agglomeration
ntaxa(physeqLung9RelOrd)

physeqLung9RelOrddf <- psmelt(physeqLung9RelOrd) # Obtain a data frame from phyloseq object.
str(physeqLung9RelOrddf)

physeqLung9RelOrddf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class, Order) %>% 
  convert(fct(OTU, Sample, Phylum, Class, Order)) %>% arrange(-(Abundance))

OrdDescBoxPlot <- physeqLung9RelOrddf %>% 
  ggplot(aes(x = reorder(Order, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Order", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Order_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
plotTreeOrd <- plot_tree(physeqLung9RelOrd, label.tips="Order",
                         color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10,
                         plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeOrdInclInf_Lung.pdf",
       width = 55, height = 25, dpi = 200, units = "cm", limitsize = FALSE)

plotTreeOrdCirc <- plot_tree(physeqLung9RelOrd, label.tips="Order", 
                             color="Filtered_Data_AK_Dec2020",
                             sizebase=10, ladderize="left", 
                             plot.margin= 0.1) +
  scale_color_manual(values = InclInfectColors) +
  coord_polar(theta="y")

ggsave("PtreeOrdInclInf_Lung_Circ.pdf",
       width = 25, height = 25, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Family-level - Total samples ----

```{r}
physeqLung9RelFam <- tax_glom(physeqLung9Rel, "Family", NArm = TRUE) # Agglomeration
ntaxa(physeqLung9RelFam)

physeqLung9RelFamdf <- psmelt(physeqLung9RelFam) # Obtain a data frame from phyloseq object.
str(physeqLung9RelFamdf)

physeqLung9RelFamdf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class, Order, Family, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(OTU, Sample, Phylum, Class, Order, Family)) %>% arrange(-(Abundance))

FamDescBoxPlot <- physeqLung9RelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Family_TotalSamples_Lung_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

FamDescBoxPlotInclInf <- FamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("Family_per_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
physeqLung9RelFamdf$Filtered_Data_AK_Dec2020 <- factor(physeqLung9RelFamdf$Filtered_Data_AK_Dec2020,
                                                       levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))
# Linear scale
physeqLung9RelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

# Log2 scale
physeqLung9RelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  ylim(c(0, 1.0)) +
  scale_y_continuous(trans = "log2") +
  coord_flip() +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("2Family_TotalSamples_InclInfGroup_Lung_Boxplot.pdf",
       width = 15, height = 45, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
plotTreeFam <- plot_tree(physeqLung9RelFam, label.tips="Family",
                         color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10,
                         plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeFamInclInf_Lung.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)

# plot_tree with shape corresponding to No Infect/Pulm Infect/Extra-pulm Infect
plotTreeFam <- plot_tree(physeqLung9RelFam, nodelabf=nodeplotblank, label.tips="Family",
                         color="Filtered_Data_AK_Dec2020", size="Abundance", sizebase=10,
                         shape="POCGroups20200115", base.spacing=0.04, plot.margin = 0.1,
                         ladderize="right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeFamInclInf_shapeNoInfPulmExtraPulmInf_Lung.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)
```

### Plot tree - Families - Total samples ----

```{r}
# Display prevalence and abundance according to the Inclusion-Infection group
plotTreeFam <- plot_tree(physeqLung9RelFam, label.tips="Family", color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10, plot.margin = 0.1, ladderize="right")

plotTreeFam + scale_color_manual(values = InclInfectColors)

ggsave("PtreeFamInclInf.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)

# plot_tree with shape corresponding to No Infect/Pulm Infect/Extra-pulm Infect
plotTreeFam <- plot_tree(physeqLung9RelFam, nodelabf=nodeplotblank, label.tips="Family",
                         color="Filtered_Data_AK_Dec2020", size="Abundance", sizebase=10,
                         shape="POCGroups20200115", base.spacing=0.04, plot.margin = 0.1, ladderize="right")

plotTreeFam + scale_color_manual(values = c("Inclusion"="gold2", "Infection_D1"="coral2",
                                            "Infection_D5"="chartreuse4",
                                            "Extubation"="cadetblue3"))

ggsave("PtreeFamInclInf_shapeNoInfPulmExtraPulmInf.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)
```

### Taxa representation across total samples ----

```{r}
physeqLung9RelFam <- tax_glom(physeqLung9Rel, "Family", NArm = TRUE) # Agglomeration

# Lachnospiraceae
physeqSubset <- subset_taxa(physeqLung9Rel, Family == "Lachnospiraceae")
physeqGlom <- tax_glom(physeqSubset, "Family", NArm = TRUE)
TaxName <- taxa_names(physeqGlom)
physeqGlomPresent <- prune_taxa(taxa_sums(physeqGlom) > 0, physeqGlom)

Lachnospiraceae <- get_sample(physeqGlomPresent, TaxName)
LachnospiraceaePresent <- Lachnospiraceae[Lachnospiraceae>0]
length(LachnospiraceaePresent)
```

### Composition analysis - Family-level - Per Phylum: Firmicutes -----

```{r}
psFirmiRel <- subset_taxa(physeqLung9Rel, Phylum == "Firmicutes", NArm = TRUE)

psFirmiRelFam <- tax_glom(psFirmiRel, "Family", NArm = TRUE) # Agglomeration
ntaxa(psFirmiRelFam)

psFirmiRelFamdf <- psmelt(psFirmiRelFam) # Obtain a data frame from phyloseq object
str(psFirmiRelFamdf)

psFirmiRelFamdf %<>% 
  select(Phylum, Family, OTU, AKSampleID, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, OTU, AKSampleID)) %>% arrange(-(Abundance))

FirmiFamDescBoxPlot <- psFirmiRelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("FirmiFamily_TotalSamples_Lung_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

FirmiFamDescBoxPlotInclInf <- FirmiFamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("FirmiFamily_per_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
psFirmiRelFamdf$Filtered_Data_AK_Dec2020 <- factor(psFirmiRelFamdf$Filtered_Data_AK_Dec2020,
                                                   levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))
# Linear scale
psFirmiRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

# Log2 scale
psFirmiRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  ylim(c(0, 1.0)) +
  scale_y_continuous(trans = "log2") +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("FirmiFamily_TotalSamples_InclInfGroup_log2_Lung_Boxplot.pdf",
       width = 15, height = 25, dpi = 200, units = "cm")

# Limit to the Top 10 abundant Firmicutes families
Top10Fam <-  sort(tapply(taxa_sums(psFirmiRelFam),
                         tax_table(psFirmiRelFam)[, "Family"], sum), decreasing = TRUE)[1:10]

# Top 10 ordered family names
Top10ListForGraph <- as.data.frame(Top10Fam) %>% rownames()
# Reverse order for vertical graph
Top10ListForGraphRev <- rev(Top10ListForGraph)

# Linear scale
psFirmiRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  xlim(Top10ListForGraphRev) +
  ylim(c(0,1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("FirmiFamily_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")

# Zoom on Firmicutes families 2 to 10
Top2to10Fam <-  sort(tapply(taxa_sums(psFirmiRelFam),
                            tax_table(psFirmiRelFam)[, "Family"], sum), decreasing = TRUE)[2:10]

# Top 2 to 10 ordered family names
Top2to10ListForGraph <- as.data.frame(Top2to10Fam) %>% rownames()
# Reverse order for vertical graph
Top2to10ListForGraphRev <- rev(Top2to10ListForGraph)

# Linear scale
psFirmiRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  xlim (Top2to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.2)) +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("Firmi2to10Family_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")
```

### Composition analysis - Family-level - Per Phylum: Proteobacteria -----

```{r}
psProteoRel <- subset_taxa(physeqLung9Rel, Phylum == "Proteobacteria", NArm = TRUE)

psProteoRelFam <- tax_glom(psProteoRel, "Family", NArm = TRUE) # Agglomeration
ntaxa(psProteoRelFam)

psProteoRelFamdf <- psmelt(psProteoRelFam) # Obtain a data frame from phyloseq object
str(psProteoRelFamdf)

psProteoRelFamdf %<>% 
  select(Phylum, Family, OTU, AKSampleID, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, OTU, AKSampleID)) %>% arrange(-(Abundance))

ProteoFamDescBoxPlot <- psProteoRelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("ProteoFamily_TotalSamples_Lung_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

ProteoFamDescBoxPlotInclInf <- ProteoFamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("ProteoFamily_per_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
psProteoRelFamdf$Filtered_Data_AK_Dec2020 <- factor(psProteoRelFamdf$Filtered_Data_AK_Dec2020,
                                                    levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))
# Linear scale - Order families by sum to show those with tipical pathogens
psProteoRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=sum), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("ProteoFamily_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

# Limit to the Top 10 abundant Proteobacteria families
Top10Fam <-  sort(tapply(taxa_sums(psProteoRelFam),
                         tax_table(psProteoRelFam)[, "Family"], sum), decreasing = TRUE)[1:10]

# Top 10 ordered family names
Top10ListForGraph <- as.data.frame(Top10Fam) %>% rownames()
# Reverse order for vertical graph
Top10ListForGraphRev <- rev(Top10ListForGraph)

# Linear scale
psProteoRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  xlim(Top10ListForGraphRev) +
  ylim(c(0,1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("ProteoFamily1to2_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")

# Zoom on Proteobacteria families 3 to 10
Top3to10Fam <-  sort(tapply(taxa_sums(psProteoRelFam),
                            tax_table(psProteoRelFam)[, "Family"], sum), decreasing = TRUE)[3:10]

# Top 3 to 10 ordered family names
Top3to10ListForGraph <- as.data.frame(Top3to10Fam) %>% rownames()
# Reverse order for vertical graph
Top3to10ListForGraphRev <- rev(Top3to10ListForGraph)

# Linear scale
psProteoRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  xlim (Top3to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.2)) +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("ProteoFamily3to10_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")
```

### Composition analysis - Family-level - Per Phylum: Bacteroidetes -----

```{r}
psBactRel <- subset_taxa(physeqLung9Rel, Phylum == "Bacteroidetes", NArm = TRUE)

psBactRelFam <- tax_glom(psBactRel, "Family", NArm = TRUE) # Agglomeration
ntaxa(psBactRelFam)

psBactRelFamdf <- psmelt(psBactRelFam) # Obtain a data frame from phyloseq object
str(psBactRelFamdf)

psBactRelFamdf %<>% 
  select(Phylum, Family, OTU, AKSampleID, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, OTU, AKSampleID)) %>% arrange(-(Abundance))

BactFamDescBoxPlot <- psBactRelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("BactFamily_TotalSamples_Lung_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

BactFamDescBoxPlotInclInf <- BactFamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("BactFamily_per_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
psBactRelFamdf$Filtered_Data_AK_Dec2020 <- factor(psBactRelFamdf$Filtered_Data_AK_Dec2020,
                                                  levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))
# Linear scale - Order families by sum to show those with tipical pathogens
psBactRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=sum), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("BactFamily_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

# Limit to the Top 10 abundant Bacteroidetes families
Top10Fam <-  sort(tapply(taxa_sums(psBactRelFam),
                         tax_table(psBactRelFam)[, "Family"], sum), decreasing = TRUE)[1:10]

# Top 10 ordered family names
Top10ListForGraph <- as.data.frame(Top10Fam) %>% rownames()
# Reverse order for vertical graph
Top10ListForGraphRev <- rev(Top10ListForGraph)

# Linear scale
psBactRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  xlim(Top10ListForGraphRev) +
  ylim(c(0,1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("BactFamily1to2_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")

# Zoom on Bacteroidetes families 3 to 10
Top3to10Fam <-  sort(tapply(taxa_sums(psBactRelFam),
                            tax_table(psBactRelFam)[, "Family"], sum), decreasing = TRUE)[3:10]

# Top 3 to 10 ordered family names
Top3to10ListForGraph <- as.data.frame(Top3to10Fam) %>% rownames()
# Reverse order for vertical graph
Top3to10ListForGraphRev <- rev(Top3to10ListForGraph)

# Linear scale
psBactRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors) +
  xlim (Top3to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.1)) +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("BactFamily3to10_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")
```

### Composition analysis - Genus-level - Total samples ----

```{r}
physeqLung9RelGen <- tax_glom(physeqLung9Rel, "Genus", NArm = TRUE) # Agglomeration
ntaxa(physeqLung9RelGen)

physeqLung9RelGendf <- psmelt(physeqLung9RelGen) # Obtain a data frame from phyloseq object.
str(physeqLung9RelGendf)

physeqLung9RelGendf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class, Order, Family, Genus, Filtered_Data_AK_Dec2020) %>%
  convert(fct(OTU, Sample, Phylum, Class, Order, Family, Genus)) %>% arrange(-(Abundance))

GenDescBoxPlot <- physeqLung9RelGendf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Genus_TotalSamples_Lung_Boxplot.pdf",
       width = 25, height = 45, dpi = 200, units = "cm")

GenDescBoxPlotInclInf <- GenDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("Genus_per_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 45, dpi = 200, units = "cm")

# plot_tree with shape corresponding to No Infect/Pulm Infect/Extra-pulm Infect
plotTreeGen <- plot_tree(physeqLung9RelGen, nodelabf=nodeplotblank, label.tips="Genus",
                         color="Filtered_Data_AK_Dec2020", size="Abundance", sizebase=10,
                         shape="POCGroups20200115", base.spacing=0.2, plot.margin = 0.2,
                         text.size=5, ladderize="right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeGenInclInf_shapeNoInfPulmExtraPulmInf_Lung.pdf",
       width = 55, height = 80, dpi = 200, units = "cm", limitsize = FALSE)


plotTreeGenCirc <- plot_tree(physeqLung9RelGen, label.tips="Genus", 
                             color="Filtered_Data_AK_Dec2020", sizebase=2, ladderize="left", 
                             plot.margin= 0) + coord_polar(theta="y")

ggsave("PtreeGenInclInf_Lung_Circ.pdf",
       width = 40, height = 40, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Species-level - Per Family: Streptococcaceae ----

```{r}
physeqSub <- subset_taxa(physeqLung9Rel, Family == "Streptococcaceae")

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeGen <- plot_tree(physeqSub, label.tips="Species",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         nodelabf = nodeplotblank,
                         plot.margin = 0.1,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeStreptoSpecInclInf_Lung.pdf",
       width = 35, height = 90, dpi = 200, units = "cm", limitsize = FALSE)

### Composition analysis - Genus-level - Per Family: Lachnospiraceae ----
physeqSub <- subset_taxa(physeqLung9RelGen, Family == "Lachnospiraceae")

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeGen <- plot_tree(physeqSub, label.tips="Genus",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         plot.margin = 0.1,
                         nodelabf = nodeplotblank,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeLachnoGenInclInf_Lung.pdf",
       width = 37, height = 6, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Species-level - Per Family: Staphylococcaceae ----

```{r}
physeqSub <- subset_taxa(physeqLung9Rel, Family == "Staphylococcaceae")
tax_table(physeqSub)
physeqSub <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSub)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSub, label.tips="Species",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeStaphyloSpecInclInf_Lung.pdf",
       width = 40, height = 6, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Species-level - Per Family: Pasteurellaceae ----

```{r}
physeqSub <- subset_taxa(physeqLung9Rel, Family == "Pasteurellaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

physeqSubHaemo <- subset_taxa(physeqSubSpe, Genus == "Haemophilus")

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSubHaemo, label.tips="Species",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeHaemophilusInclInf_Lung.pdf",
       width = 43, height = 7, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Genus-level - Per Family: Enterobacteraceae ----

```{r}
physeqSub <- subset_taxa(physeqLung9Rel, Family == "Enterobacteriaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSubSpe, label.tips="Species",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeEnterobactSpecInclInf_Lung.pdf",
       width = 40, height = 6, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Genus-level - Per Family: Prevotelleceae ----

```{r}
physeqSub <- subset_taxa(physeqLung9Rel, Family == "Prevotellaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSubSpe, label.tips="Genus",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreePrevotSpecInclInf_Lung.pdf",
       width = 45, height = 65, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Patients with Pulm. infection AND (InfD1 + InfD5) available ----

```{r}
PatPulmInf_InfD1InfD5 <- c("2", "4", "8", "28", "30", "36", "37", "41", "43")

physPulmInf_InfD1InfD5 <- subset_samples(physeqLung9, record_ID %in% PatPulmInf_InfD1InfD5)

# Focus on Proteobacteria genera
InfD1InfD5Prot <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Proteobacteria")

InfD1InfD5ProtGen <- tax_glom(InfD1InfD5Prot, "Genus", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5ProtPresent <- prune_taxa(taxa_sums(InfD1InfD5ProtGen) > 0, 
                                    InfD1InfD5ProtGen)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5ProtGenRel <- transform_sample_counts(InfD1InfD5ProtPresent, count_to_rel_abund)

InfD1InfD5ProtGenReldf <- psmelt(InfD1InfD5ProtGenRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5ProtGenReldf)

InfD1InfD5ProtGenReldf %<>% 
  select(Sample, Phylum, Genus, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Genus, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5ProtGenReldf)
tail(InfD1InfD5ProtGenReldf)

InfD1InfD5ProtGenBoxPlot <- InfD1InfD5ProtGenReldf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#1F78B4", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_ProteoGenera_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5ProtGenBoxPlot_perInfGroup <- InfD1InfD5ProtGenBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_ProteoGenera_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 15, dpi = 200, units = "cm")

# Focus on Firmicutes families
InfD1InfD5Firm <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Firmicutes")

InfD1InfD5FirmFam <- tax_glom(InfD1InfD5Firm, "Family", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5FirmFamPresent <- prune_taxa(taxa_sums(InfD1InfD5FirmFam) > 0, 
                                       InfD1InfD5FirmFam)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5FirmFamRel <- transform_sample_counts(InfD1InfD5FirmFamPresent, count_to_rel_abund)

InfD1InfD5FirmFamReldf <- psmelt(InfD1InfD5FirmFamRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5FirmFamReldf)

InfD1InfD5FirmFamReldf %<>% 
  select(Sample, Phylum, Family, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Family, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5FirmFamReldf)
tail(InfD1InfD5FirmFamReldf)

InfD1InfD5FirmFamBoxPlot <- InfD1InfD5FirmFamReldf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#33A02C", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Family" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_FirmFamily_PulmInf_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5FirmFamBoxPlot_perInfGroup <- InfD1InfD5FirmFamBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_FirmFamilies_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Focus on Firmicutes genera
InfD1InfD5Firm <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Firmicutes")

InfD1InfD5FirmGen <- tax_glom(InfD1InfD5Firm, "Genus", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5FirmGenPresent <- prune_taxa(taxa_sums(InfD1InfD5FirmGen) > 0, 
                                       InfD1InfD5FirmGen)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5FirmGenRel <- transform_sample_counts(InfD1InfD5FirmGenPresent, count_to_rel_abund)

InfD1InfD5FirmGenReldf <- psmelt(InfD1InfD5FirmGenRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5FirmGenReldf)

InfD1InfD5FirmGenReldf %<>% 
  select(Sample, Phylum, Family, Genus, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Family, Genus, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5FirmGenReldf)
tail(InfD1InfD5FirmGenReldf)

InfD1InfD5FirmGenBoxPlot <- InfD1InfD5FirmGenReldf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#33A02C", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_FirmGenera_PulmInf_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5FirmGenBoxPlot_perInfGroup <- InfD1InfD5FirmGenBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_FirmFamilies_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Focus on Bacteroidetes families
InfD1InfD5Bact <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Bacteroidetes")

InfD1InfD5BactFam <- tax_glom(InfD1InfD5Bact, "Family", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5BactFamPresent <- prune_taxa(taxa_sums(InfD1InfD5BactFam) > 0, 
                                       InfD1InfD5BactFam)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5BactFamRel <- transform_sample_counts(InfD1InfD5BactFamPresent, count_to_rel_abund)

InfD1InfD5BactFamReldf <- psmelt(InfD1InfD5BactFamRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5BactFamReldf)

InfD1InfD5BactFamReldf %<>% 
  select(Sample, Phylum, Family, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Family, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5BactFamReldf)
tail(InfD1InfD5BactFamReldf)

InfD1InfD5BactFamBoxPlot <- InfD1InfD5BactFamReldf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#FF7F00", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Family" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_BactFamily_PulmInf_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5BactFamBoxPlot_perInfGroup <- InfD1InfD5BactFamBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_BactFamilies_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 7, dpi = 200, units = "cm")

# Focus on Bacteroidetes genera
InfD1InfD5Bact <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Bacteroidetes")

InfD1InfD5BactGen <- tax_glom(InfD1InfD5Bact, "Genus", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5BactPresent <- prune_taxa(taxa_sums(InfD1InfD5BactGen) > 0, 
                                    InfD1InfD5BactGen)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5BactGenRel <- transform_sample_counts(InfD1InfD5BactPresent, count_to_rel_abund)

InfD1InfD5BactGenReldf <- psmelt(InfD1InfD5BactGenRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5BactGenReldf)

InfD1InfD5BactGenReldf %<>% 
  select(Sample, Phylum, Genus, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Genus, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5BactGenReldf)
tail(InfD1InfD5BactGenReldf)

InfD1InfD5BactGenBoxPlot <- InfD1InfD5BactGenReldf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#FF7F00", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_BactGenera_InclInfGroups_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5BactGenBoxPlot_perInfGroup <- InfD1InfD5BactGenBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_BactGenera_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 15, dpi = 200, units = "cm")

### Inclusion group: Phylum-level analysis ----
Inclusion <- subset_samples(physeqLung9Rel, Filtered_Data_AK_Dec2020=="Inclusion")
InclusionPhy <- tax_glom(Inclusion, "Phylum", NArm = TRUE) # Agglomeration

InclusionPhydf <- psmelt(InclusionPhy) # Obtain a data frame from phyloseq object.
str(InclusionPhydf)

InclusionPhydf %<>% 
  select(OTU, Sample, Abundance, Phylum) %>% 
  convert(fct(OTU, Sample, Phylum)) %>% arrange(-(Abundance))

InclPhyDescBoxPlot <- InclusionPhydf %>% 
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum, alpha=0.9)) +
  scale_fill_manual(values = PhylumPalLung13) +
  ylim(c(0, 1.0)) +
  labs(x = "Phyla", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Inclusion_Phyla_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")
```

### Phyloseq heatmaps ----

```{r}
# For plot_heatmap analysis, the sample and taxa names in secondary data, 
# OTU table and taxa table must be changed.
# Include a column in sec data with PatientID_Sample_InclInfGr for 
# labeling the x axis of heatmap
physeqLung9Dat_df <- data.frame(sample_data(physeqLung9))
head(physeqLung9Dat_df)
rownames(physeqLung9Dat_df)

ForUnXlabeldf <- physeqLung9Dat_df %>% 
  unite(UnXlabel, SeqSampleID, AKSampleID, record_ID, Filtered_Data_AK_Dec2020, sep = "_")

physeqLung9Dat_df$UnXlabel <- ForUnXlabeldf$UnXlabel
head(physeqLung9Dat_df)

# For labeling the y axis of heatmap: 
# Attach taxa names to OTU names in taxa table
physeqLung9Tax_df <- data.frame(tax_table(physeqLung9))
head(physeqLung9Tax_df)
rownames(physeqLung9Tax_df)
physeqLung9Tax2_df <- physeqLung9Tax_df %>% rownames_to_column("SpID")
head(physeqLung9Tax2_df)

ForUnYlabeldf <- physeqLung9Tax2_df %>% 
  unite(UnYlabel, SpID, Phylum, Family, Genus, Species, sep = "_") %>% 
  arrange(UnYlabel)
head(ForUnYlabeldf)

physeqLung9Tax2Sor_df <- physeqLung9Tax2_df %>% arrange(SpID)
head(physeqLung9Tax2Sor_df)

physeqLung9Tax2Sor_df$UnYlabel <- ForUnYlabeldf$UnYlabel
head(physeqLung9Tax2Sor_df)

rownames(physeqLung9Tax2Sor_df) <- physeqLung9Tax2Sor_df$UnYlabel
head(physeqLung9Tax2Sor_df)

# Change rownames in OTU table to match taxa table
physeqLung9OTU_df <- data.frame(otu_table(physeqLung9))
physeqLung9OTU2Sor_df <- physeqLung9OTU_df %>% 
  rownames_to_column("SpID") %>% 
  arrange(SpID)
head(physeqLung9OTU2Sor_df)
dim(physeqLung9OTU2Sor_df)
dim(physeqLung9Tax2Sor_df)

rownames(physeqLung9OTU2Sor_df) <- rownames(physeqLung9Tax2Sor_df)

# Remove SpID column to obtain further below a matrix
physeqLung9OTU2Sor_df %<>% select(-1)
str(physeqLung9OTU2Sor_df)

# Prepare phyloseq object for heatmap plotting
physeqLung9UnDat <- sample_data(physeqLung9Dat_df)
class(physeqLung9UnDat)
str(physeqLung9UnDat)
sample_names(physeqLung9UnDat)

physeqLung9OTU2 <- otu_table(physeqLung9OTU2Sor_df, taxa_are_rows = TRUE)
class(physeqLung9OTU2)
sample_names(physeqLung9OTU2)

physeqLung9Tax2 <- tax_table(as.matrix(physeqLung9Tax2Sor_df))
class(physeqLung9Tax2)

physeqLung9Un <- phyloseq(physeqLung9UnDat, physeqLung9OTU2, physeqLung9Tax2)

# Plot for figure - 2021.01.21
PHeat1_20210121 <- plot_heatmap(physeqLung9Un, "NMDS", "bray", sample.label="UnXlabel", 
                                taxa.label="UnYlabel", max.label = 1989)

ggsave(filename = "PHeat1_20210121_high_tax_names.pdf",
       width = 40, height = 450, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

ggsave(filename = "PHeat1_20210121.pdf",
       width = 25, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

PHeat1_20210121[["plot_env"]][["sample.order"]]
PHeat1_20210121[["plot_env"]][["taxa.order"]]
tail(PHeat1_20210121[["plot_env"]][["taxa.order"]])

PHeat2 <- plot_heatmap(physeqLung9Un, "NMDS", "jaccard", sample.label="UnXlabel", 
                       species.label="Genus", max.label = 1989)

ggsave(filename = "PHeat1_20210121.pdf",
       width = 25, height = 600, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

PHeat1facetFiltered_Data_AK_Dec2020 <- PHeat1_20210121 + facet_grid(~Filtered_Data_AK_Dec2020,
                                                                    scales = "free_x")

plot(PHeat1facetFiltered_Data_AK_Dec2020)

ggsave(filename = "PHeat1facetFiltered_Data_AK_Dec2020.pdf",
       width = 25, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

PHeat1facetPatient <- PHeat1_20210121 + facet_grid(~record_ID, scales = "free_x")

plot(PHeat1facetPatient)

ggsave(filename = "PHeat1facetPatient.pdf",
       width = 50, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)
```

### Venn diagram - Total samples at Inclusion versus Infection_D5 ----

```{r}
# Get taxa names at Inclusion
sample_data(physeqLung9)
physeqLung9Incl <- subset_samples(physeqLung9, Filtered_Data_AK_Dec2020 %in% 
                                    "Inclusion")
ntaxa(physeqLung9Incl)

# Keep only represented taxa
physeqLung9InclPresent <- prune_taxa(taxa_sums(physeqLung9Incl) > 0, 
                                     physeqLung9Incl)

ntaxa(physeqLung9InclPresent)

physeqLung9InclPresent_df <- data.frame(tax_table(physeqLung9InclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(physeqLung9InclPresent_df) # enter in area 1 in venn diagram
head(physeqLung9InclPresent_df)

# Get taxa names at Infection_D5
physeqLung9InfD5 <- subset_samples(physeqLung9, Filtered_Data_AK_Dec2020 %in% 
                                     "Infection_D5")
ntaxa(physeqLung9InfD5)

# Keep only represented taxa
physeqLung9InfD5Present <- prune_taxa(taxa_sums(physeqLung9InfD5) > 0, 
                                      physeqLung9InfD5)

ntaxa(physeqLung9InfD5Present)

physeqLung9InfD5Present_df <- data.frame(tax_table(physeqLung9InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(physeqLung9InfD5Present_df) # enter in area 2 in venn diagram
head(physeqLung9InfD5Present_df)

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physeqLung9InclPresent_df, 
                                 physeqLung9InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
head(taxInclusionOnly_df)

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physeqLung9InfD5Present_df, 
                             physeqLung9InclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
head(taxInfD5Only_df)

# Taxa present at Inclusion and Infection_D5
taxIncl_AND_InfD5 <- intersect(physeqLung9InclPresent_df$ASV, physeqLung9InfD5Present_df$ASV)
length(taxIncl_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create Venn diagram
VennDiagram::draw.pairwise.venn(area1 = 1012,area2 = 544, cross.area = 247,
                                fill = c("Inclusion"="gold2", "Infection D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Inclusion versus Infection_D1 ----

```{r}
Pat3 <- subset_samples(physeqLung9Un, record_ID %in% "3")

# Get taxa names at Inclusion
sample_data(Pat3)
Pat3Incl <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                             "Inclusion")
ntaxa(Pat3Incl)

# Keep only represented taxa
Pat3InclPresent <- prune_taxa(taxa_sums(Pat3Incl) > 0, 
                              Pat3Incl)
# n taxa at Inclusion + Inclusion AND Infection_D1
ntaxa(Pat3InclPresent)

Pat3InclPresent_df <- data.frame(tax_table(Pat3InclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InclPresent_df) # enter in area 1 in venn diagram
head(Pat3InclPresent_df)

# Get taxa names at Infection_D1
Pat3InfD1 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D1")
ntaxa(Pat3InfD1)

# Keep only represented taxa
Pat3InfD1Present <- prune_taxa(taxa_sums(Pat3InfD1) > 0, 
                               Pat3InfD1)
# n taxa at Infection_D1 + Inclusion AND Infection_D1
ntaxa(Pat3InfD1Present)

Pat3InfD1Present_df <- data.frame(tax_table(Pat3InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD1Present_df) # enter in area 2 in venn diagram
head(Pat3InfD1Present_df)

# Get taxa names present at Inclusion ONLY
Pat3taxInclusionOnly_df <- anti_join(Pat3InclPresent_df, 
                                     Pat3InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInclusionOnly_df)
head(Pat3taxInclusionOnly_df)

# Get taxa names present at Infection_D1 ONLY
Pat3taxInfD1Only_df <- anti_join(Pat3InfD1Present_df, 
                                 Pat3InclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD1Only_df)
head(Pat3taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
Pat3taxIncl_AND_InfD1 <- intersect(Pat3InclPresent_df$ASV, Pat3InfD1Present_df$ASV)
length(Pat3taxIncl_AND_InfD1) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 152,area2 = 144, cross.area = 94,
                                fill = c("Inclusion"="gold2", "Infection_D1"="coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D1"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD1_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Infection_D1 vs Infection_D5 ----

```{r}
# Get taxa names at Infection_D1
sample_data(Pat3)
Pat3InfD1 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D1")
ntaxa(Pat3InfD1)

# Keep only represented taxa
Pat3InfD1Present <- prune_taxa(taxa_sums(Pat3InfD1) > 0, 
                               Pat3InfD1)
# n taxa at Infection_D1 + Infection_D1 AND Infection_D5
ntaxa(Pat3InfD1Present)

Pat3InfD1Present_df <- data.frame(tax_table(Pat3InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD1Present_df) # enter in area 1 in venn diagram
head(Pat3InfD1Present_df)

# Get taxa names at Infection_D5
Pat3InfD5 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D5")
ntaxa(Pat3InfD5)

# Keep only represented taxa
Pat3InfD5Present <- prune_taxa(taxa_sums(Pat3InfD5) > 0, 
                               Pat3InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Infection_D1
ntaxa(Pat3InfD5Present)

Pat3InfD5Present_df <- data.frame(tax_table(Pat3InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD5Present_df) # enter in area 2 in venn diagram
head(Pat3InfD5Present_df)

# Get taxa names present at Infection_D1 ONLY
Pat3taxInfD1Only_df <- anti_join(Pat3InfD1Present_df, 
                                 Pat3InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD1Only_df)
head(Pat3taxInfD1Only_df)

# Get taxa names present at Infection_D5 ONLY
Pat3taxInfD5Only_df <- anti_join(Pat3InfD5Present_df, 
                                 Pat3InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD5Only_df)
head(Pat3taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
Pat3taxInfD1_AND_InfD5 <- intersect(Pat3InfD1Present_df$ASV, Pat3InfD5Present_df$ASV)
length(Pat3taxInfD1_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 144,area2 = 46, cross.area = 25,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Infection_D5 vs Extubation ----

```{r}
# Get taxa names at Infection_D5
sample_data(Pat3)
Pat3InfD5 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D5")
ntaxa(Pat3InfD5)

# Keep only represented taxa
Pat3InfD5Present <- prune_taxa(taxa_sums(Pat3InfD5) > 0, 
                               Pat3InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Extubation
ntaxa(Pat3InfD5Present)

Pat3InfD5Present_df <- data.frame(tax_table(Pat3InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD5Present_df) # enter in area 1 in venn diagram
head(Pat3InfD5Present_df)

# Get taxa names at Extubation
Pat3Extub <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Extubation")
ntaxa(Pat3Extub)

# Keep only represented taxa
Pat3ExtubPresent <- prune_taxa(taxa_sums(Pat3Extub) > 0, 
                               Pat3Extub)
# n taxa at Extubation + Extubation AND Infection_D5
ntaxa(Pat3ExtubPresent)

Pat3ExtubPresent_df <- data.frame(tax_table(Pat3ExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3ExtubPresent_df) # enter in area 2 in venn diagram
head(Pat3ExtubPresent_df)

# Get taxa names present at Infection_D5 ONLY
Pat3taxInfD5Only_df <- anti_join(Pat3InfD5Present_df, 
                                 Pat3ExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD5Only_df)
head(Pat3taxInfD5Only_df)

# Get taxa names present at Extubation ONLY
Pat3taxExtubOnly_df <- anti_join(Pat3ExtubPresent_df, 
                                 Pat3InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxExtubOnly_df)
head(Pat3taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
Pat3taxInfD5_AND_Extub <- intersect(Pat3InfD5Present_df$ASV, Pat3ExtubPresent_df$ASV)
length(Pat3taxInfD5_AND_Extub) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 46,area2 = 45, cross.area = 18,
                                fill = c("Infection_D5"="chartreuse4", "Extubation"="cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D5", "Extubation"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD5Extub_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Inclusion vs Infection_D5 ----

```{r}
# Get taxa names at Inclusion
sample_data(Pat3)
Pat3Incl <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                             "Inclusion")
ntaxa(Pat3Incl)

# Keep only represented taxa
Pat3InclPresent <- prune_taxa(taxa_sums(Pat3Incl) > 0, 
                              Pat3Incl)
# n taxa at Inclusion + Inclusion AND Infection_D5
ntaxa(Pat3InclPresent)

Pat3InclPresent_df <- data.frame(tax_table(Pat3InclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InclPresent_df) # enter in area 1 in venn diagram
head(Pat3InclPresent_df)

# Get taxa names at Infection_D5
Pat3InfD5 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D5")
ntaxa(Pat3InfD5)

# Keep only represented taxa
Pat3InfD5Present <- prune_taxa(taxa_sums(Pat3InfD5) > 0, 
                               Pat3InfD5)
# n taxa at Infection_D5 + Inclusion AND Infection_D5
ntaxa(Pat3InfD5Present)

Pat3InfD5Present_df <- data.frame(tax_table(Pat3InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD5Present_df) # enter in area 2 in venn diagram
head(Pat3InfD5Present_df)

# Get taxa names present at Inclusion ONLY
Pat3taxInclusionOnly_df <- anti_join(Pat3InclPresent_df, 
                                     Pat3InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInclusionOnly_df)
head(Pat3taxInclusionOnly_df)

# Get taxa names present at Infection_D5 ONLY
Pat3taxInfD5Only_df <- anti_join(Pat3InfD5Present_df, 
                                 Pat3InclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD5Only_df)
head(Pat3taxInfD5Only_df)

# Taxa present at Inclusion and Infection_D5
Pat3taxIncl_AND_InfD5 <- intersect(Pat3InclPresent_df$ASV, Pat3InfD5Present_df$ASV)
length(Pat3taxIncl_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 152,area2 = 46, cross.area = 20,
                                fill = c("Inclusion"="gold2", "Infection D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Pulmonary infection Patient 28 - Inclusion versus Infection_D1 ----

```{r}
Pat28 <- subset_samples(physeqLung9Un, record_ID %in% "28")

# Get taxa names at Inclusion
sample_data(Pat28)
Pat28Incl <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                              "Inclusion")
ntaxa(Pat28Incl)

# Keep only represented taxa
Pat28InclPresent <- prune_taxa(taxa_sums(Pat28Incl) > 0, 
                               Pat28Incl)
# n taxa at Inclusion + Inclusion AND Infection_D1
ntaxa(Pat28InclPresent)

Pat28InclPresent_df <- data.frame(tax_table(Pat28InclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InclPresent_df) # enter in area 1 in venn diagram
head(Pat28InclPresent_df)

# Get taxa names at Infection_D1
Pat28InfD1 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D1")
ntaxa(Pat28InfD1)

# Keep only represented taxa
Pat28InfD1Present <- prune_taxa(taxa_sums(Pat28InfD1) > 0, 
                                Pat28InfD1)
# n taxa at Infection_D1 + Inclusion AND Infection_D1
ntaxa(Pat28InfD1Present)

Pat28InfD1Present_df <- data.frame(tax_table(Pat28InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD1Present_df) # enter in area 2 in venn diagram
head(Pat28InfD1Present_df)

# Get taxa names present at Inclusion ONLY
Pat28taxInclusionOnly_df <- anti_join(Pat28InclPresent_df, 
                                      Pat28InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInclusionOnly_df)
head(Pat28taxInclusionOnly_df)

# Get taxa names present at Infection_D1 ONLY
Pat28taxInfD1Only_df <- anti_join(Pat28InfD1Present_df, 
                                  Pat28InclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD1Only_df)
head(Pat28taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
Pat28taxIncl_AND_InfD1 <- intersect(Pat28InclPresent_df$ASV, Pat28InfD1Present_df$ASV)
length(Pat28taxIncl_AND_InfD1) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 104,area2 = 39, cross.area = 13,
                                fill = c("Inclusion"="gold2", "Infection_D1"="coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D1"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD1_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Pulmonary infection Patient 28 - Infection_D1 versus Infection_D5 ----

```{r}
# Get taxa names at Infection_D1
sample_data(Pat28)
Pat28InfD1 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D1")
ntaxa(Pat28InfD1)

# Keep only represented taxa
Pat28InfD1Present <- prune_taxa(taxa_sums(Pat28InfD1) > 0, 
                                Pat28InfD1)
# n taxa at Infection_D1 + Infection_D1 AND Infection_D5
ntaxa(Pat28InfD1Present)

Pat28InfD1Present_df <- data.frame(tax_table(Pat28InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD1Present_df) # enter in area 1 in venn diagram
head(Pat28InfD1Present_df)

# Get taxa names at Infection_D5
Pat28InfD5 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D5")
ntaxa(Pat28InfD5)

# Keep only represented taxa
Pat28InfD5Present <- prune_taxa(taxa_sums(Pat28InfD5) > 0, 
                                Pat28InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Infection_D1
ntaxa(Pat28InfD5Present)

Pat28InfD5Present_df <- data.frame(tax_table(Pat28InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD5Present_df) # enter in area 2 in venn diagram
head(Pat28InfD5Present_df)

# Get taxa names present at Infection_D1 ONLY
Pat28taxInfD1Only_df <- anti_join(Pat28InfD1Present_df, 
                                  Pat28InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD1Only_df)
head(Pat28taxInfD1Only_df)

# Get taxa names present at Infection_D5 ONLY
Pat28taxInfD5Only_df <- anti_join(Pat28InfD5Present_df, 
                                  Pat28InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD5Only_df)
head(Pat28taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
Pat28taxInfD1_AND_InfD5 <- intersect(Pat28InfD1Present_df$ASV, Pat28InfD5Present_df$ASV)
length(Pat28taxInfD1_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 39,area2 = 38, cross.area = 8,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Pulmonary infection Patient 28 - Infection_D5 versus Extubation ----

```{r}
# Get taxa names at Infection_D5
sample_data(Pat28)
Pat28InfD5 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D5")
ntaxa(Pat28InfD5)

# Keep only represented taxa
Pat28InfD5Present <- prune_taxa(taxa_sums(Pat28InfD5) > 0, 
                                Pat28InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Extubation
ntaxa(Pat28InfD5Present)

Pat28InfD5Present_df <- data.frame(tax_table(Pat28InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD5Present_df) # enter in area 1 in venn diagram
head(Pat28InfD5Present_df)

# Get taxa names at Extubation
Pat28Extub <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Extubation")
ntaxa(Pat28Extub)

# Keep only represented taxa
Pat28ExtubPresent <- prune_taxa(taxa_sums(Pat28Extub) > 0, 
                                Pat28Extub)
# n taxa at Extubation + Extubation AND Infection_D5
ntaxa(Pat28ExtubPresent)

Pat28ExtubPresent_df <- data.frame(tax_table(Pat28ExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28ExtubPresent_df) # enter in area 2 in venn diagram
head(Pat28ExtubPresent_df)

# Get taxa names present at Infection_D5 ONLY
Pat28taxInfD5Only_df <- anti_join(Pat28InfD5Present_df, 
                                  Pat28ExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD5Only_df)
head(Pat28taxInfD5Only_df)

# Get taxa names present at Extubation ONLY
Pat28taxExtubOnly_df <- anti_join(Pat28ExtubPresent_df, 
                                  Pat28InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxExtubOnly_df)
head(Pat28taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
Pat28taxInfD5_AND_Extub <- intersect(Pat28InfD5Present_df$ASV, Pat28ExtubPresent_df$ASV)
length(Pat28taxInfD5_AND_Extub) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 38,area2 = 24, cross.area = 2,
                                fill = c("Infection_D5"="chartreuse4", "Extubation"="cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D5", "Extubation"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with No infection AND (Inclusion + Extub) available ----

```{r}
# Compare Inclusion with Extubation
PatNoInf_InclExtub <- c("7","20","24","39")

physNoInf_InclExtub <- subset_samples(physeqLung9, record_ID %in% PatNoInf_InclExtub)

# Keep only represented taxa
physNoInf_InclExtubPresent <- prune_taxa(taxa_sums(physNoInf_InclExtub) > 0, 
                                         physNoInf_InclExtub)

TotalTaxaPresent <- ntaxa(physNoInf_InclExtubPresent)

# Subset taxa present either at Inclusion or (Inclusion AND Extubation)
physNoInfIncl <- subset_samples(physNoInf_InclExtubPresent,
                                Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physNoInfInclPresent <- prune_taxa(taxa_sums(physNoInfIncl) > 0, 
                                   physNoInfIncl)

physNoInfInclPresent_df <- data.frame(tax_table(physNoInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physNoInfInclPresent_df)[1] # for area 1 in venn diagram
head(physNoInfInclPresent_df)

# Subset taxa present either at Extubation or (Extubation AND Inclusion)
physNoInfExtub <- subset_samples(physNoInf_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                   "Extubation")

# Keep only represented taxa
physNoInfExtubPresent <- prune_taxa(taxa_sums(physNoInfExtub) > 0, 
                                    physNoInfExtub)

physNoInfExtubPresent_df <- data.frame(tax_table(physNoInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physNoInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physNoInfInclPresent_df, 
                                 physNoInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physNoInfExtubPresent_df, 
                             physNoInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physNoInfInclPresent_df$ASV, physNoInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                fill = c("gold2", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInclExtubAbs_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Inclusion + Extub) available ----

```{r}
# Compare Inclusion with Extubation
PatPulmInf_InclExtub <- c("28")

physPulmInf_InclExtub <- subset_samples(physeqLung9, record_ID %in% PatPulmInf_InclExtub)

# Keep only represented taxa
physPulmInf_InclExtubPresent <- prune_taxa(taxa_sums(physPulmInf_InclExtub) > 0, 
                                           physPulmInf_InclExtub)

TotalTaxaPresent <- ntaxa(physPulmInf_InclExtubPresent)

# Subset taxa present either in Inclusion or Inclusion AND Extubation
physPulmInfIncl <- subset_samples(physPulmInf_InclExtubPresent,
                                  Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physPulmInfInclPresent <- prune_taxa(taxa_sums(physPulmInfIncl) > 0, 
                                     physPulmInfIncl)

physPulmInfInclPresent_df <- data.frame(tax_table(physPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physPulmInfInclPresent_df)

# Subset taxa present either in Extubation or Extubation AND Inclusion
physPulmInfExtub <- subset_samples(physPulmInf_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                     "Extubation")

# Keep only represented taxa
physPulmInfExtubPresent <- prune_taxa(taxa_sums(physPulmInfExtub) > 0, 
                                      physPulmInfExtub)

physPulmInfExtubPresent_df <- data.frame(tax_table(physPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physPulmInfInclPresent_df, 
                                 physPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physPulmInfExtubPresent_df, 
                             physPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physPulmInfInclPresent_df$ASV, physPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                fill = c("gold2", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInclExtubAbs_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Inclusion + Extub) available ----

```{r}
# Compare Inclusion with Extubation
PatExtraPulmInf_InclExtub <- c("3", "16", "18", "32")

physExtraPulmInf_InclExtub <- subset_samples(physeqLung9, record_ID %in% PatExtraPulmInf_InclExtub)

# Keep only represented taxa
physExtraPulmInf_InclExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInf_InclExtub) > 0, 
                                                physExtraPulmInf_InclExtub)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InclExtubPresent)

# Subset taxa present either in Inclusion or Inclusion AND Extubation
physExtraPulmInfIncl <- subset_samples(physExtraPulmInf_InclExtubPresent,
                                       Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physExtraPulmInfInclPresent <- prune_taxa(taxa_sums(physExtraPulmInfIncl) > 0, 
                                          physExtraPulmInfIncl)

physExtraPulmInfInclPresent_df <- data.frame(tax_table(physExtraPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInclPresent_df)

# Subset taxa present either in Extubation or Extubation AND Inclusion
physExtraPulmInfExtub <- subset_samples(physExtraPulmInf_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                          "Extubation")

# Keep only represented taxa
physExtraPulmInfExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInfExtub) > 0, 
                                           physExtraPulmInfExtub)

physExtraPulmInfExtubPresent_df <- data.frame(tax_table(physExtraPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physExtraPulmInfInclPresent_df, 
                                 physExtraPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physExtraPulmInfExtubPresent_df, 
                             physExtraPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physExtraPulmInfInclPresent_df$ASV, physExtraPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                fill = c("gold2", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInclExtubAbs_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Infection_D1 + Infection_D5) available ----

```{r}
PatPulmInf_InfD1InfD5 <- c("2", "4", "8", "28", "30", "36", "37", "41", "43")

physPulmInf_InfD1InfD5 <- subset_samples(physeqLung9, record_ID %in% PatPulmInf_InfD1InfD5)

# Keep only represented taxa
physPulmInf_InfD1InfD5Present <- prune_taxa(taxa_sums(physPulmInf_InfD1InfD5) > 0, 
                                            physPulmInf_InfD1InfD5)

TotalTaxaPresent <- ntaxa(physPulmInf_InfD1InfD5Present)

# Subset taxa present either at Infection_D1 or (Infection_D1 AND Infection_D5)
physPulmInfInfD1 <- subset_samples(physPulmInf_InfD1InfD5Present,
                                   Filtered_Data_AK_Dec2020 %in% "Infection_D1")

# Keep only represented taxa
physPulmInfInfD1Present <- prune_taxa(taxa_sums(physPulmInfInfD1) > 0, 
                                      physPulmInfInfD1)

physPulmInfInfD1Present_df <- data.frame(tax_table(physPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInfD1Present_df)[1] # for area 1 in venn diagram
head(physPulmInfInfD1Present_df)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Infection_D1)
physPulmInfInfD5 <- subset_samples(physPulmInf_InfD1InfD5Present, Filtered_Data_AK_Dec2020 %in% 
                                     "Infection_D5")

# Keep only represented taxa
physPulmInfInfD5Present <- prune_taxa(taxa_sums(physPulmInfInfD5) > 0, 
                                      physPulmInfInfD5)

physPulmInfInfD5Present_df <- data.frame(tax_table(physPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfInfD5Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physPulmInfInfD1Present_df, 
                             physPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
taxInfD1Only_df[1:10,]

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physPulmInfInfD5Present_df, 
                             physPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
head(taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
taxInfD1_AND_InfD5 <- intersect(physPulmInfInfD1Present_df$ASV, physPulmInfInfD5Present_df$ASV)
Venn_cross_area <- length(taxInfD1_AND_InfD5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Infection_D1 + Infection_D5) available ----

```{r}
PatExtraPulmInf_InfD1InfD5 <- c("3", "6", "13", "16", "18", "23", "29", "32", "34")

physExtraPulmInf_InfD1InfD5 <- subset_samples(physeqLung9, record_ID %in% PatExtraPulmInf_InfD1InfD5)

# Keep only represented taxa
physExtraPulmInf_InfD1InfD5Present <- prune_taxa(taxa_sums(physExtraPulmInf_InfD1InfD5) > 0, 
                                                 physExtraPulmInf_InfD1InfD5)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InfD1InfD5Present)

# Subset taxa present either at Infection_D1 or (Infection_D1 AND Infection_D5)
physExtraPulmInfInfD1 <- subset_samples(physExtraPulmInf_InfD1InfD5Present,
                                        Filtered_Data_AK_Dec2020 %in% "Infection_D1")

# Keep only represented taxa
physExtraPulmInfInfD1Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD1) > 0, 
                                           physExtraPulmInfInfD1)

physExtraPulmInfInfD1Present_df <- data.frame(tax_table(physExtraPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInfD1Present_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInfD1Present_df)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Infection_D1)
physExtraPulmInfInfD5 <- subset_samples(physExtraPulmInf_InfD1InfD5Present, Filtered_Data_AK_Dec2020 %in% 
                                          "Infection_D5")

# Keep only represented taxa
physExtraPulmInfInfD5Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD5) > 0, 
                                           physExtraPulmInfInfD5)

physExtraPulmInfInfD5Present_df <- data.frame(tax_table(physExtraPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfInfD5Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physExtraPulmInfInfD1Present_df, 
                             physExtraPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
taxInfD1Only_df[1:10,]

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physExtraPulmInfInfD5Present_df, 
                             physExtraPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
head(taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
taxInfD1_AND_InfD5 <- intersect(physExtraPulmInfInfD1Present_df$ASV, physExtraPulmInfInfD5Present_df$ASV)
Venn_cross_area <- length(taxInfD1_AND_InfD5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Inclusion + Infection_D1) available ----

```{r}
PatPulmInf_InclInfD1 <- c("1", "9", "28", "33", "40", "41")

physPulmInf_InclInfD1 <- subset_samples(physeqLung9, record_ID %in% PatPulmInf_InclInfD1)

# Keep only represented taxa
physPulmInf_InclInfD1Present <- prune_taxa(taxa_sums(physPulmInf_InclInfD1) > 0, 
                                           physPulmInf_InclInfD1)

TotalTaxaPresent <- ntaxa(physPulmInf_InclInfD1Present)

# Subset taxa present either at Inclusion or (Inclusion AND Infection_D1)
physPulmInfIncl <- subset_samples(physPulmInf_InclInfD1Present,
                                  Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physPulmInfInclPresent <- prune_taxa(taxa_sums(physPulmInfIncl) > 0, 
                                     physPulmInfIncl)

physPulmInfInclPresent_df <- data.frame(tax_table(physPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physPulmInfInclPresent_df)

# Subset taxa present either at Infection_D1 or (Infection_D1 AND Inclusion)
physPulmInfInfD1 <- subset_samples(physPulmInf_InclInfD1Present, Filtered_Data_AK_Dec2020 %in% 
                                     "Infection_D1")

# Keep only represented taxa
physPulmInfInfD1Present <- prune_taxa(taxa_sums(physPulmInfInfD1) > 0, 
                                      physPulmInfInfD1)

physPulmInfInfD1Present_df <- data.frame(tax_table(physPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfInfD1Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclOnly_df <- anti_join(physPulmInfInclPresent_df, 
                            physPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclOnly_df)
taxInclOnly_df[1:10,]

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physPulmInfInfD1Present_df, 
                             physPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
head(taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
taxIncl_AND_InfD1 <- intersect(physPulmInfInclPresent_df$ASV, physPulmInfInfD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_InfD1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Inclusion" = "gold2", "Infection_D1" = "coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection_D1"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Inclusion + Infection_D1) available ----

```{r}
PatExtraPulmInf_InclInfD1 <- c("3", "16", "18", "32")

physExtraPulmInf_InclInfD1 <- subset_samples(physeqLung9, record_ID %in% PatExtraPulmInf_InclInfD1)

# Keep only represented taxa
physExtraPulmInf_InclInfD1Present <- prune_taxa(taxa_sums(physExtraPulmInf_InclInfD1) > 0, 
                                                physExtraPulmInf_InclInfD1)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InclInfD1Present)

# Subset taxa present either at Inclusion or (Inclusion AND Infection_D1)
physExtraPulmInfIncl <- subset_samples(physExtraPulmInf_InclInfD1Present,
                                       Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physExtraPulmInfInclPresent <- prune_taxa(taxa_sums(physExtraPulmInfIncl) > 0, 
                                          physExtraPulmInfIncl)

physExtraPulmInfInclPresent_df <- data.frame(tax_table(physExtraPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInclPresent_df)

# Subset taxa present either at Infection_D1 or (Infection_D1 AND Inclusion)
physExtraPulmInfInfD1 <- subset_samples(physExtraPulmInf_InclInfD1Present, Filtered_Data_AK_Dec2020 %in% 
                                          "Infection_D1")

# Keep only represented taxa
physExtraPulmInfInfD1Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD1) > 0, 
                                           physExtraPulmInfInfD1)

physExtraPulmInfInfD1Present_df <- data.frame(tax_table(physExtraPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfInfD1Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclOnly_df <- anti_join(physExtraPulmInfInclPresent_df, 
                            physExtraPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclOnly_df)
taxInclOnly_df[1:10,]

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physExtraPulmInfInfD1Present_df, 
                             physExtraPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
head(taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
taxIncl_AND_InfD1 <- intersect(physExtraPulmInfInclPresent_df$ASV, physExtraPulmInfInfD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_InfD1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Inclusion" = "gold2", "Infection_D1" = "coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection_D1"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Infection_D5 + Extubation) available ----

```{r}
PatPulmInf_InfD5Extub <- c("2", "4", "8", "28", "30", "37")

physPulmInf_InfD5Extub <- subset_samples(physeqLung9, record_ID %in% PatPulmInf_InfD5Extub)

# Keep only represented taxa
physPulmInf_InfD5ExtubPresent <- prune_taxa(taxa_sums(physPulmInf_InfD5Extub) > 0, 
                                            physPulmInf_InfD5Extub)

TotalTaxaPresent <- ntaxa(physPulmInf_InfD5ExtubPresent)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Extubation)
physPulmInfInfD5 <- subset_samples(physPulmInf_InfD5ExtubPresent,
                                   Filtered_Data_AK_Dec2020 %in% "Infection_D5")

# Keep only represented taxa
physPulmInfInfD5Present <- prune_taxa(taxa_sums(physPulmInfInfD5) > 0, 
                                      physPulmInfInfD5)

physPulmInfInfD5Present_df <- data.frame(tax_table(physPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInfD5Present_df)[1] # for area 1 in venn diagram
head(physPulmInfInfD5Present_df)

# Subset taxa present either at Extubation or (Extubation AND Infection_D5)
physPulmInfExtub <- subset_samples(physPulmInf_InfD5ExtubPresent, Filtered_Data_AK_Dec2020 %in% 
                                     "Extubation")

# Keep only represented taxa
physPulmInfExtubPresent <- prune_taxa(taxa_sums(physPulmInfExtub) > 0, 
                                      physPulmInfExtub)

physPulmInfExtubPresent_df <- data.frame(tax_table(physPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physPulmInfInfD5Present_df, 
                             physPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
taxInfD5Only_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physPulmInfExtubPresent_df, 
                             physPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
taxInfD5_AND_Extub <- intersect(physPulmInfInfD5Present_df$ASV, physPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxInfD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Infection_D5" = "chartreuse4", "Extubation" = "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection_D5", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD5Extub_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Infection_D5 + Extubation) available ----

```{r}
PatExtraPulmInf_InfD5Extub <- c("3", "13", "16", "18", "25", "32") 

physExtraPulmInf_InfD5Extub <- subset_samples(physeqLung9, record_ID %in% PatExtraPulmInf_InfD5Extub)

# Keep only represented taxa
physExtraPulmInf_InfD5ExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInf_InfD5Extub) > 0, 
                                                 physExtraPulmInf_InfD5Extub)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InfD5ExtubPresent)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Extubation)
physExtraPulmInfInfD5 <- subset_samples(physExtraPulmInf_InfD5ExtubPresent,
                                        Filtered_Data_AK_Dec2020 %in% "Infection_D5")

# Keep only represented taxa
physExtraPulmInfInfD5Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD5) > 0, 
                                           physExtraPulmInfInfD5)

physExtraPulmInfInfD5Present_df <- data.frame(tax_table(physExtraPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInfD5Present_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInfD5Present_df)

# Subset taxa present either at Extubation or (Extubation AND Infection_D5)
physExtraPulmInfExtub <- subset_samples(physExtraPulmInf_InfD5ExtubPresent,
                                        Filtered_Data_AK_Dec2020 %in% "Extubation")

# Keep only represented taxa
physExtraPulmInfExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInfExtub) > 0, 
                                           physExtraPulmInfExtub)

physExtraPulmInfExtubPresent_df <- data.frame(tax_table(physExtraPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physExtraPulmInfInfD5Present_df, 
                             physExtraPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
taxInfD5Only_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physExtraPulmInfExtubPresent_df, 
                             physExtraPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
taxInfD5_AND_Extub <- intersect(physExtraPulmInfInfD5Present_df$ASV, physExtraPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxInfD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Infection_D5" = "chartreuse4", "Extubation" = "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection_D5", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD5Extub_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Heatmap Patient3 with represented taxa only ----

```{r}
Pat3Present <- prune_taxa(taxa_sums(Pat3) > 0, Pat3)
Pat3Present_PHeat1_20210125 <- plot_heatmap(Pat3Present, "NMDS", "bray", sample.label="UnXlabel", 
                                            taxa.label="UnYlabel", max.label = 240)

ggsave(filename = "Pat3Present_PHeat1_20210125_high_tax_names.pdf",
       width = 15, height = 45, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

ggsave(filename = "Pat3_PHeat1_20210125.pdf",
       width = 12, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

Pat3Present_PHeat1_20210125[["plot_env"]][["sample.order"]]
Pat3Present_PHeat1_20210125[["plot_env"]][["taxa.order"]]
tail(Pat3Present_PHeat1_20210125[["plot_env"]][["taxa.order"]])
```

### Heatmap Patient28 with represented taxa only ----

```{r}
Pat28Present <- prune_taxa(taxa_sums(Pat28) > 0, Pat28)
Pat28Present_PHeat1_20210127 <- plot_heatmap(Pat28Present, "NMDS", "bray", sample.label="UnXlabel", 
                                             taxa.label="UnYlabel", max.label = 166)

ggsave(filename = "Pat28Present_PHeat1_20210127_high_tax_names.pdf",
       width = 15, height = 45, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

ggsave(filename = "Pat28_PHeat1_20210127.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

Pat28Present_PHeat1_20210127[["plot_env"]][["sample.order"]]
Pat28Present_PHeat1_20210127[["plot_env"]][["taxa.order"]]
tail(Pat28Present_PHeat1_20210127[["plot_env"]][["taxa.order"]])
```

### Alpha diversity - Richness ----

/#/`/`/`{r} psLungAlphaRich_df /<- psmelt(physeqLung9Rel) View(psLungAlphaRich_df)

psLungAlphaRich_df %/</>% select(OTU, Sample, Abundance, Phylum, Filtered_Data_AK_Dec2020, Shannon, InvSimpson,Chao1, Observed,POCGroups20200115) %/>% convert(fct(OTU, Sample, Phylum, Filtered_Data_AK_Dec2020)) %/>% arrange(-(Abundance))

Shannon_lung/<-psLungAlphaRich_df %/>% ggplot(aes(x =POCGroups20200115 , y = Shannon)) + geom_boxplot( aes(fill =Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Group of patients", y = "Shannon") + theme_classic()

ggsave("Shannon_Lung_Groups.tiff", width = 14, height = 9, dpi = 300, units = "cm")

psLungAlphaRich_df %/>% ggplot(aes(x = Filtered_Data_AK_Dec2020, y = InvSimpson)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Sample collection phase", y = "Inverse Simpson") + theme_classic()

ggsave("InvSimpson_InclInfGroups.pdf", width = 9, height = 9, dpi = 200, units = "cm")

Chao1_lung/<-psLungAlphaRich_df %/>% ggplot(aes(x = POCGroups20200115, y = Chao1)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Patient group", y = "Chao1") + theme_classic()

ggsave("Chao1_Lung_Groups.tiff", width = 14, height = 9, dpi = 300, units = "cm")

ggarrange(Shannon_lung,Chao1_lung, labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )

ggsave("Alpha_diversity_combined_lung.tiff", width = 18, height = 9, dpi = 300, units = "cm")

# Ploting of alpha diversity directly as calculated by the pipeline output.

# Plot first for lung microbiome

SecDataLungFil$Filtered_Data_AK_Dec2020<-as.factor(SecDataLungFil$Filtered_Data_AK_Dec2020) levels(SecDataLungFil/$Filtered_Data_AK_Dec2020)

levels(SecDataLungFil/$Filtered_Data_AK_Dec2020)/<-c("Extubation","Inclusion","Infection_D1","Infection_D5") /#rename factor levels

SecDataLungFil$Filtered_Data_AK_Dec2020<-factor(SecDataLungFil$Filtered_Data_AK_Dec2020,levels=c("Inclusion","Infection_D1","Infection_D5","Extubation")) /#reorder levels of the factor by chronological order

levels(SecDataLungFil/$Filtered_Data_AK_Dec2020)

/#Check alpha diversity only on inclusion among groups a_inclusion_shannon/<-SecDataLungFil %/>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%/>% ggplot(aes(x = POCGroups20200115, y = Shannon,label=Sample)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + geom_jitter(position = position_jitter(seed = 1))+ geom_text(position = position_jitter(seed = 1),hjust=0,vjust=0)+ scale_fill_manual(values = InclInfectColors4) + labs(x = "Group ofpatients", y = "Shannon",fill="Timepoint") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_inclusion_shannon

ggsave("a_inclusion_gut.tiff", width = 25, height = 25, dpi = 300, units = "cm")

a_inclusion_chao1/<-SecDataLungFil %/>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%/>% ggplot(aes(x = POCGroups20200115, y = Chao1,label=Sample)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + geom_jitter(position = position_jitter(seed = 1))+ geom_text(position = position_jitter(seed = 1),hjust=0,vjust=0)+ scale_fill_manual(values = InclInfectColors4) + labs(x = "Group of patients", y = "Chao1",fill="Timepoint") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y =180)

a_inclusion_chao1

ggsave("a_inclusion_gut_chao.tiff", width = 25, height = 25, dpi = 300, units = "cm")

/#Check alpha diversity only on extubation among groups

a_extubation_shannon/<-SecDataLungFil %/>%subset(Filtered_Data_AK_Dec2020=="Extubation")%/>% ggplot(aes(x = POCGroups20200115, y = Shannon)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Extubation", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_extubation_shannon

a_extubation_chao1/<-SecDataLungFil %/>%subset(Filtered_Data_AK_Dec2020=="Extubation")%/>% ggplot(aes(x = POCGroups20200115, y = Chao1)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Extubation", y = "Chao1") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 180)

a_extubation_chao1

ggarrange(a_inclusion_shannon,a_inclusion_chao1,a_extubation_shannon,a_extubation_chao1, labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )

ggsave("a_inclusion and extubation_lung.tiff", width = 25, height = 25, dpi = 300, units = "cm")

# Plot alpha diversity upon time for all groups with Shannon and Chao1 indexes

a_shannon_lung/<-SecDataLungFil %/>% ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Sample collection phase", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_shannon_lung+facet_wrap(/~POCGroups20200115)

ggsave("Alpha_diversity_combined_lung_shannon.tiff", width = 25, height = 15, dpi = 300, units = "cm")

a_chao1_lung/<-SecDataLungFil %/>% ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Sample collection phase", y = "Chao1") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 200)

a_chao1_lung

a_chao1_lung+facet_wrap(/~POCGroups20200115)

ggsave("Alpha_diversity_combined_lung_Chao1.tiff", width = 25, height = 15, dpi = 300, units = "cm")

# Repeat measurements of alpha-diversity as above for Gut microbiome

SecDataGutFil$Filtered_Data_AK_Dec2020<-as.factor(SecDataGutFil$Filtered_Data_AK_Dec2020)

levels(SecDataGutFil/$Filtered_Data_AK_Dec2020)

levels(SecDataGutFil/$Filtered_Data_AK_Dec2020)/<-c("Discharge","Extubation","Inclusion","Infection_D1","Infection_D5") /#rename factor levels

SecDataGutFil$Filtered_Data_AK_Dec2020<-factor(SecDataGutFil$Filtered_Data_AK_Dec2020,levels=c("Inclusion","Infection_D1","Infection_D5","Extubation","Discharge")) /#reorder levels of the factor by chronological order

levels(SecDataGutFil/$Filtered_Data_AK_Dec2020)

View(SecDataGutFil)

/#Check alpha diversity only on inclusion among groups metadata_310122 %/>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%/>% ggplot(aes(x = POCGroups20200115, y = Shannon,label=Sample)) + geom_boxplot( aes(fill = POCGroups20200115), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + geom_jitter(position = position_jitter(seed = 1))+ labs(x = "Inclusion", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

metadata_310122 %/>%subset(Subgroup=="Inclusion")%/>% ggplot(aes(x = POCGroups20200115, y = Shannon,label=Sample)) + geom_boxplot( aes(fill = POCGroups20200115), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + geom_jitter(position = position_jitter(seed = 1))+ labs(x = "Inclusion", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

metadata_310122 %/>% ggplot(aes(x = Group, y = Shannon,label=Sample)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + geom_jitter(position = position_jitter(seed = 1))+ scale_fill_manual(values = InclInfectColors5)+ labs(x = "Inclusion", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_inclusion_gut_chao1/<-SecDataGutFil %/>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%/>% ggplot(aes(x = POCGroups20200115, y = Shannon)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + geom_jitter(position = position_jitter(seed = 1))+ scale_fill_manual(values = InclInfectColors4) + labs(x = "Inclusion", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y =500)

a_inclusion_gut_chao1

/#Check alpha diversity only on extubation among groups and/or discharge

a_extubation_gut_shannon/<-SecDataGutFil %/>%subset(Filtered_Data_AK_Dec2020=="Extubation")%/>% ggplot(aes(x = POCGroups20200115, y = Shannon)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Extubation", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_extubation_gut_shannon

a_discharge_gut_shannon/<-SecDataGutFil %/>%subset(Filtered_Data_AK_Dec2020=="Discharge")%/>% ggplot(aes(x = POCGroups20200115, y = Shannon)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Discharge", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_discharge_gut_shannon

/#Combined extubation and discharge a_last_samples_gut_shannon/<-SecDataGutFil %/>% subset((Filtered_Data_AK_Dec2020=="Extubation")/|(Filtered_Data_AK_Dec2020=="Discharge"))%/>% ggplot(aes(x = POCGroups20200115, y = Shannon)) + geom_boxplot( aes(fill = POCGroups20200115), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Discharge", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_last_samples_gut_shannon

a_extubation_gut_chao1/<-SecDataGutFil %/>%subset(Filtered_Data_AK_Dec2020=="Extubation")%/>% ggplot(aes(x = POCGroups20200115, y = Chao1)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors4) + labs(x = "Extubation", y = "Chao1") + theme_classic()+ stat_compare_means(label.x = 2, label.y = 450)

a_extubation_gut_chao1

a_discharge_gut_chao1/<-SecDataGutFil %/>%subset(Filtered_Data_AK_Dec2020=="Discharge")%/>% ggplot(aes(x = POCGroups20200115, y = Chao1)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Discharge", y = "Chao1") + theme_classic()+ stat_compare_means(label.x = 2, label.y = 450)

a_discharge_gut_chao1

/#Combined extubation and discharge a_last_samples_gut_chao1/<-SecDataGutFil %/>% subset((Filtered_Data_AK_Dec2020=="Extubation")/|(Filtered_Data_AK_Dec2020=="Discharge"))%/>% ggplot(aes(x = POCGroups20200115, y = Chao1)) + geom_boxplot( aes(fill = POCGroups20200115), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Discharge", y = "Chao1") + theme_classic()+ stat_compare_means(label.x = 1.5, label.y = 5)

a_last_samples_gut_chao1

ggarrange(a_inclusion_gut_shannon,a_inclusion_gut_chao1,a_extubation_gut_shannon,a_extubation_gut_chao1, labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )

ggsave("a_inclusion and extubation_gut.tiff", width = 25, height = 25, dpi = 300, units = "cm")

ggarrange(a_discharge_gut_shannon,a_discharge_gut_chao1,a_last_samples_gut_shannon,a_last_samples_gut_chao1, labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )

ggsave("a_last_samples(extubation and discharge together)/_gut.tiff", width = 25, height = 25, dpi = 300, units = "cm")

# Plot alpha diversity upon time for all groups with Shannon and Chao1 indexes

a_shannon_gut/<-SecDataGutFil %/>% ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Sample collection phase", y = "Shannon") + theme_classic()+ stat_compare_means(label.x = 2.5, label.y = 5.5)

a_shannon_gut

a_shannon_gut+facet_wrap(/~POCGroups20200115)

ggsave("Alpha_diversity_combined_gut_shannon.tiff", width = 25, height = 15, dpi = 300, units = "cm")

a_chao1_gut/<-SecDataGutFil %/>% ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Sample collection phase", y = "Chao1") + theme_classic()+ stat_compare_means(label.x = 2.5, label.y = 550)

a_chao1_gut

a_chao1_gut+facet_wrap(/~POCGroups20200115)

ggsave("Alpha_diversity_combined_gut_Chao1.tiff", width = 25, height = 15, dpi = 300, units = "cm")

### Distance - beta diversity ----
#```{r}
# Adapted from https://joey711.github.io/phyloseq/distance.html
# Visualize sample distribution according to various distance metrics
dist_methods <- unlist(phyloseq::distanceMethodList)
print(dist_methods)

dist_methods <- dist_methods[-which(dist_methods=="ANY")]

plist <- vector("list", length(dist_methods))
names(plist) <-  dist_methods
names(plist)

for( i in dist_methods ){
  # Calculate distance matrix (distance function requires phyloseq::)
  iDist <- phyloseq::distance(physeqLung9, method=i)
  # Calculate ordination
  iMDS  <- ordinate(physeqLung9, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(physeqLung9, iMDS, color="Filtered_Data_AK_Dec2020")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
multidiv <- ggplot(df, aes(Axis.1, Axis.2, color=Filtered_Data_AK_Dec2020)) + 
  geom_point(size=3, alpha=0.5) + 
  facet_wrap(~distance, scales="free") +
  ggtitle("MDS on various distance metrics")

multidiv

ggsave(filename = "Beta_diversity_multidiv.pdf",
       width = 45, height = 45, dpi = 200, units = "cm", device='pdf')

print(plist[["binary"]])
print(plist[["unifrac"]])

distBin <- phyloseq::distance(physeqLung9, method="binary")

physeqLung9_ord <- ordinate(physeqLung9, "MDS", distBin)
plot_ord_MDS_bin <- plot_ordination(physeqLung9, physeqLung9_ord, type="samples",
                                    color="Filtered_Data_AK_Dec2020") +
  scale_color_manual(values = InclInfectColors4) +
  theme_bw()

ggsave(filename = "plot_ord_MDS_bin_InclInf_Lung.pdf",
       width = 12, height = 6, dpi = 200, units = "cm", device='pdf')

# Assign a different shape to No infection, Pneumonia and Extra-pulmonary infection 
plot_ord_MDS_bin2 <- plot_ordination(physeqLung9, physeqLung9_ord, type="samples",
                                     axes=c(1,2), color="Filtered_Data_AK_Dec2020", shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors4) +
  theme_bw()

ggsave(filename = "plot_ord_MDS_bin2_InclInf_shape_Lung.pdf",
       width = 12, height = 6, dpi = 200, units = "cm", device='pdf')

# wunifrac distance

View(psmelt(physeqLung9))

physeqLung9_ord2 <- ordinate(physeqLung9, "NMDS", "wunifrac")
plot_ord_MDS_wunifrac <- plot_ordination(physeqLung9, physeqLung9_ord2, type="samples",
                                         axes=c(1,2), color="Filtered_Data_AK_Dec2020", shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors4) +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Weighted Unifrac",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_wunifrac 

plot_ord_MDS_wunifrac+stat_ellipse(type="t")  # ad ellipses

wunifraclung<-plot_ord_MDS_wunifrac+stat_ellipse(type="t")+facet_wrap(~POCGroups20200115) # combined graphs

wunifraclung


ggsave(filename = "plot_ord_MDS_unifrac_InclInf_Lung.pdf",
       width = 12, height = 6, dpi = 200, units = "cm", device='pdf')

sample_data(physeqLung9)$record_ID %>% unique() %>% length()


# Bray-Curtis distance
physeqLung9_ord3 <- ordinate(physeqLung9, "NMDS", "bray")
plot_ord_MDS_bray <- plot_ordination(physeqLung9, physeqLung9_ord3, type="samples",
                                     axes=c(1,2), color="Filtered_Data_AK_Dec2020", shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors4) +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Bray-Curtis",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_bray 

plot_ord_MDS_bray+stat_ellipse(type="t")  # ad ellipses

braylung<-plot_ord_MDS_bray+stat_ellipse(type="t")+facet_wrap(~POCGroups20200115) # combined graphs

braylung

# Taxa ordination according to patient group based on Bray-Curtis distance
taxordlung<-plot_ordination(physeqLung9, physeqLung9_ord3, type="split", 
                            color="Phylum",title="Taxa ordination (distance:Bray-Curtis)",shape="POCGroups20200115")

taxordlung+facet_wrap(~Phylum)



# Jaccard distance
physeqLung9_ord4 <- ordinate(physeqLung9, "NMDS", "jaccard")
plot_ord_MDS_jaccard <- plot_ordination(physeqLung9, physeqLung9_ord4, type="samples",
                                        axes=c(1,2), color="Filtered_Data_AK_Dec2020", shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors4) +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Jaccard",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_jaccard 

plot_ord_MDS_jaccard+stat_ellipse(type="t")  # ad ellipses

jaccardlung<-plot_ord_MDS_jaccard+stat_ellipse(type="t")+facet_wrap(~POCGroups20200115) # combined graphs

jaccardlung

ggarrange(wunifraclung,jaccardlung,braylung,labels=c("A","B","C"),common.legend=TRUE,nrow=3,legend="right")

ggsave(filename = "beta_diversity_lung_combined.pdf",
       width =24 , height = 20, dpi = 300, units = "cm", device='pdf')



###  Check for beta-diversity of the 3 groups only at inclusion

physeqLung9incl<-subset_samples(physeqLung9,Filtered_Data_AK_Dec2020=="Inclusion")



# Bray distance on inclusion
physeqLung9_ordincl <- ordinate(physeqLung9incl, "NMDS", "bray")

bray_lung_incl <- plot_ordination(physeqLung9incl, physeqLung9_ordincl, type="samples", axes=c(1,2), color="POCGroups20200115") +
  scale_color_brewer(palette="Dark2") +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Lung samples",subtitle="Bray-Curtis distance only at Inclusion",colour="Group of patients")

bray_lung_incl_elps<-bray_lung_incl+stat_ellipse(type="t")  # ad ellipses
bray_lung_incl_elps

# Repeat the same graph for jaccard

physeqLung9_ordincl <- ordinate(physeqLung9incl, "NMDS", "jaccard")

jaccard_lung_incl <- plot_ordination(physeqLung9incl, physeqLung9_ordincl, type="samples", axes=c(1,2), color="POCGroups20200115") +
  scale_color_brewer(palette="Dark2") +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Lung samples",subtitle="Jaccard distance only at Inclusion",colour="Group of patients")

jaccard_lung_incl_elps<-jaccard_lung_incl+stat_ellipse(type="t")  # ad ellipses

jaccard_lung_incl_elps




#Combine graphs
ggarrange(bray_lung_incl_elps,jaccard_lung_incl_elps,
          labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )

ggsave("Beta_diversity_lung_inclusion.tiff",
       width = 24, height = 9, dpi = 300, units = "cm")

### Differential abundance - DESeq2 - Total patients - Inclusion vs Infection_D5 ----

```{r}
# From https://joey711.github.io/phyloseq-extensions/DESeq2.html
# Subset samples of Inclusion and InfectD5
InclInfD5_ps <- subset_samples(physeqLung9Un,
                               Filtered_Data_AK_Dec2020 %in% c("Inclusion", "Infection_D5"))

InclInfD5_present <- prune_taxa(taxa_sums(InclInfD5_ps)>0, InclInfD5_ps)

diagdds <-  phyloseq_to_deseq2(InclInfD5_present, ~ Filtered_Data_AK_Dec2020)

# calculate geometric means prior to estimate size factors
gm_mean  <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(diagdds), 1, gm_mean)
diagdds  <-  estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds  <-  DESeq(diagdds, fitType="local")

res <-  results(diagdds, cooksCutoff = FALSE)
alpha  <-  0.05
sigtab  <-  res[which(res$padj < alpha), ]
sigtab  <-  cbind(as(sigtab, "data.frame"), as(tax_table(InclInfD5_present)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_classic())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, colour=Phylum)) +
  geom_point(size=4, alpha=0.7) +
  scale_colour_manual(values = PhylumPal) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(filename = "DESeq2_InclInfD5_Lung.pdf",
       width = 15, height = 12, dpi = 200, units = "cm", device='pdf')
```

### Differential abundance - DESeq2 - Patients with Pulmonary infection AND (Infection_D1, Infection_D5 AND/OR Extubation) available ----

```{r}
# From https://joey711.github.io/phyloseq-extensions/DESeq2.html
# Subset samples of Inclusion and InfectD5
PatPulmInf_InfD1InfD5 <- c("2", "4", "8", "28", "30", "36", "37", "41", "43")

physPulmInf_InfD1InfD5 <- subset_samples(physeqLung9, record_ID %in% PatPulmInf_InfD1InfD5)

# Keep only represented taxa
physPulmInf_InfD1InfD5Present <- prune_taxa(taxa_sums(physPulmInf_InfD1InfD5) > 0, 
                                            physPulmInf_InfD1InfD5)

diagdds <-  phyloseq_to_deseq2(physPulmInf_InfD1InfD5Present, ~ Filtered_Data_AK_Dec2020)

# calculate geometric means prior to estimate size factors
gm_mean  <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(diagdds), 1, gm_mean)
diagdds  <-  estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds  <-  DESeq(diagdds, fitType="local")

res <-  results(diagdds, cooksCutoff = FALSE)
alpha  <-  0.05
sigtab  <-  res[which(res$padj < alpha), ]
sigtab  <-  cbind(as(sigtab, "data.frame"), as(tax_table(physPulmInf_InfD1InfD5Present)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_classic())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, colour=Phylum)) +
  geom_point(size=2.5, alpha=0.7) +
  scale_colour_manual(values = PhylumPal) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(filename = "DESeq2_PulmInf_InfD1InfD5_Lung.pdf",
       width = 15, height = 12, dpi = 200, units = "cm", device='pdf')
```

### Differential abundance - DESeq2 - Patients with No infection AND (Inclusion + Extub) available ----

```{r}
# From https://joey711.github.io/phyloseq-extensions/DESeq2.html
# Subset samples
PatNoInf_InclExtub <- c("7","20","24","39")

physNoInf_InclExtub <- subset_samples(physeqLung9, record_ID %in% PatNoInf_InclExtub)

# Keep only represented taxa
physNoInf_InclExtubPresent <- prune_taxa(taxa_sums(physNoInf_InclExtub) > 0, 
                                         physNoInf_InclExtub)

sample_data(physNoInf_InclExtubPresent)

diagdds <-  phyloseq_to_deseq2(physNoInf_InclExtubPresent, ~ Filtered_Data_AK_Dec2020)

# calculate geometric means prior to estimate size factors
gm_mean  <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(diagdds), 1, gm_mean)
diagdds  <-  estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds  <-  DESeq(diagdds, fitType="local")

res <-  results(diagdds, cooksCutoff = FALSE)
alpha  <-  0.05 # significance level
sigtab  <-  res[which(res$padj < alpha), ]
sigtab  <-  cbind(as(sigtab, "data.frame"), as(tax_table(physNoInf_InclExtubPresent)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_classic())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, colour=Phylum)) +
  geom_point(size=3, alpha=0.6) +
  scale_colour_manual(values = PhylumPal) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(filename = "DESeq2_NoInf_InclExtub_Lung.pdf",
       width = 15, height = 12, dpi = 200, units = "cm", device='pdf')
```

### Taxa unique to Patients with No infection AND (Inclusion + Extub) available ----

```{r}
# See Venn diagram - Patients with No infection AND (Inclusion + Extub) available
# Inclusion only
InclOnlyASVNames <- taxInclusionOnly_df$ASV

physInclOnly <- prune_taxa(InclOnlyASVNames, physNoInfIncl)
tax_table(physInclOnly)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physInclOnlyRel <- transform_sample_counts(physInclOnly, count_to_rel_abund)

# Agglomeration
physInclOnlyAgg <- tax_glom(physInclOnlyRel, "Genus", NArm = TRUE)

# Top 20 genera of taxa exclusively represented at Inclusion
TopTaxaNames <-  names(sort(taxa_sums(physInclOnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physInclOnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPalLung13) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusInclOnly_Lung_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Extubation only
ExtubOnlyASVNames <- taxExtubOnly_df$ASV

physExtubOnly <- prune_taxa(ExtubOnlyASVNames, physNoInfExtub)
tax_table(physExtubOnly)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physExtubOnlyRel <- transform_sample_counts(physExtubOnly, count_to_rel_abund)

# Agglomeration
physExtubOnlyAgg <- tax_glom(physExtubOnlyRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa exclusively represented at Extubation
TopTaxaNames <-  names(sort(taxa_sums(physExtubOnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physExtubOnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPalLung13) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusExtubOnly_Lung_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Taxa shared between Inclusion and Extubation
physShared <- prune_taxa(taxIncl_AND_Extub, physNoInf_InclExtubPresent)
tax_table(physShared)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physSharedRel <- transform_sample_counts(physShared, count_to_rel_abund)

# Agglomeration
physSharedAgg <- tax_glom(physSharedRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa shared between Inclusion and Extubation
TopTaxaNames <-  names(sort(taxa_sums(physSharedAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physSharedAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPalLung13) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusSharedInclExtub_Lung_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")
```

### Taxa unique to Patients with Pulm infection AND (Infection_D1 + Infection_D5) available ----

```{r}
# See Venn diagram - Pulm infection AND (Infection_D1 + Infection_D5) available
# Infection_D1 only
InfD1OnlyASVNames <- taxInfD1Only_df$ASV

physInfD1Only <- prune_taxa(InfD1OnlyASVNames, physPulmInfInfD1Present)
tax_table(physInfD1Only)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physInfD1OnlyRel <- transform_sample_counts(physInfD1Only, count_to_rel_abund)

# Agglomeration
physInfD1OnlyAgg <- tax_glom(physInfD1OnlyRel, "Genus", NArm = TRUE)

# Top 20 genera of taxa exclusively represented at Infection_D1
TopTaxaNames <-  names(sort(taxa_sums(physInfD1OnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physInfD1OnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPal) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusPulmInf_InfD1Only_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Infection_D5 only
InfectD5OnlyASVNames <- taxInfD5Only_df$ASV

physInfD5Only <- prune_taxa(InfectD5OnlyASVNames, physPulmInfInfD5Present)
tax_table(physInfD5Only)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physInfD5OnlyRel <- transform_sample_counts(physInfD5Only, count_to_rel_abund)

# Agglomeration
physInfD5OnlyAgg <- tax_glom(physInfD5OnlyRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa exclusively represented at Infection_D5
TopTaxaNames <-  names(sort(taxa_sums(physInfD5OnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physInfD5OnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPal) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusPulmInf_InfD5Only_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Taxa shared between Infection_D1 and Infection_D5
physShared <- prune_taxa(taxInfD1_AND_InfD5, physPulmInf_InfD1InfD5Present)
tax_table(physShared)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physSharedRel <- transform_sample_counts(physShared, count_to_rel_abund)

# Agglomeration
physSharedAgg <- tax_glom(physSharedRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa shared between Infection_D1 and Infection_D5
TopTaxaNames <-  names(sort(taxa_sums(physSharedAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physSharedAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPal) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusPulmInf_SharedInfD1InfD5_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")
```

#### Gut data analysis / Community abundance manipulation and exploration ----

```{r}
physeqGut1

OTU_Gut1 <- otu_table(physeqGut4)
sam_Gut1 <- sample_data(physeqGut4)
tax_Gut1 <- tax_table(physeqGut4)
tree_Gut1 <- phy_tree(physeqGut4)
refseq_Gut1 <- refseq(physeqGut4)

get_taxa_unique(physeqGut4, "Kingdom")
get_taxa_unique(physeqGut4, "Phylum")
get_taxa_unique(physeqGut4, "Order")
get_taxa_unique(physeqGut4, "Family")
get_taxa_unique(physeqGut4, "Genus")
get_taxa_unique(physeqGut4, "Species")

sort(sample_sums(physeqGut4))

# Top100 most abundant ASVs
Top100AbundantASVCountsGut <- sort(taxa_sums(physeqGut4), decreasing=TRUE)[1:100]
Top100AbundantASVNamesGut  <-  names(Top100AbundantASVCountsGut)
physeGut4Top100  <-  prune_taxa(Top100AbundantASVNamesGut, physeqGut4)
# Visualize taxonomy of most abundant ASVs
Top100AbundantTaxaGut <- tax_table(physeGut4Top100)[ , c("Phylum", "Genus" ,"Species")]

# Number of total ASVs per phylum 
TablePhylaTaxa <- table(tax_table(physeqGut4)[ , 2])
sort(TablePhylaTaxa, decreasing=TRUE)

# Set a palette for Phyla with colours visible to colour-blind people
display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 12, name = "Paired")
brewer.pal(n = 8, name = "Dark2")
PhylumPalGut13 <- c("Firmicutes"="#33A02C", "Fusobacteria"="#CAB2D6",
                    "Actinobacteria"="#E31A1C", "Bacteroidetes"="#FF7F00",
                    "Proteobacteria"="#1F78B4", "Saccharibacteria_TM7"="#FB9A99",
                    "Spirochaetes"="#FDBF6F", "Synergistetes"="#A6CEE3",
                    "Tenericutes"="#B2DF8A", "Cyanobacteria"="#FFFF99",
                    "Lentisphaerae"="#6A3D9A", "Verrucomicrobia"="#B15928",
                    "Deinococcus-Thermus"="#666666")

# Transform counts to relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqGut4Rel <- transform_sample_counts(physeqGut4, count_to_rel_abund)
# Check total = 1 for the first 10 samples
sample_sums(physeqGut4Rel)[1:10]
```

#### Composition analysis - Phylum-level ----

```{r}
# Focus on human samples
get_variable(physeqGut4, varName = "Filtered_Data_AK_Dec2020")
physeqGutHum <- subset_samples(physeqGut4, !Filtered_Data_AK_Dec2020 == "Neg control")

# Transform counts to relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqGutHumRel <- transform_sample_counts(physeqGutHum, count_to_rel_abund)

# Phylum agglomeration
physeqGutHumRelPhy <- tax_glom(physeqGutHumRel, "Phylum", NArm = TRUE) # Agglomeration

physeqGutHumRelPhydf <- psmelt(physeqGutHumRelPhy) # Obtain a data frame from phyloseq object.
str(physeqGutHumRelPhydf)

physeqGutHumRelPhydf %<>% 
  select(OTU, Sample, Abundance, Phylum, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(OTU, Sample, Phylum, Filtered_Data_AK_Dec2020)) %>% arrange(-(Abundance))

PhyDescBoxPlot <- physeqGutHumRelPhydf %>% 
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill = Phylum), alpha = 0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Phylum" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Phylum_TotalSamples_Gut_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

InclInfPhyDescBoxPlot <- PhyDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("Phylum_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 30, height = 8, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
physeqGutHumRelPhydf$Filtered_Data_AK_Dec2020 <- factor(physeqGutHumRelPhydf$Filtered_Data_AK_Dec2020,
                                                        levels = c("Discharge", "Extubation",
                                                                   "Infection D5", "Infection D1", "Inclusion"))

physeqGutHumRelPhydf %>%
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("Phylum_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 15, height = 18, dpi = 200, units = "cm")

# Zoom on Phyla families 4 to 10
Top4to10Phy <-  sort(tapply(taxa_sums(physeqGutHumRelPhy),
                            tax_table(physeqGutHumRelPhy)[, "Phylum"], sum), decreasing = TRUE)[4:10]

# Top 4 to 10 ordered phylum names
Top4to10ListForGraph <- as.data.frame(Top4to10Phy) %>% rownames()
# Reverse order for vertical graph
Top4to10ListForGraphRev <- rev(Top4to10ListForGraph)

physeqGutHumRelPhydf %>%
  ggplot(aes(x = Phylum, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim (Top4to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.25)) +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("Phylum4to10_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
# Construct a highly agglomerated tree to obtain a reference for further comparisons
tax_table(physeqGutHumRelPhy)

# Merge all taxa except one poorly represented phylum
MergedPhy <- merge_taxa(physeqGutHumRelPhy, taxa_names(physeqGutHumRelPhy)[-c(2)], 1)
tax_table(MergedPhy)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeGen <- plot_tree(MergedPhy, label.tips="Phylum",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         plot.margin = 0.1,
                         nodelabf = nodeplotblank,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePhyMergedInclInf_Gut.pdf",
       width = 85, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Phylum-level tree
plotTreePhy <- plot_tree(physeqGutHumRelPhy, label.tips="Phylum",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         plot.margin = 0.1,
                         nodelabf = nodeplotblank,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePhyInclInf_Gut.pdf",
       width = 65, height = 15, dpi = 200, units = "cm", limitsize = FALSE)
```

##### Lung-Gut comparisons - Venn diagram - All patients and samples ----

```{r}
# Subset Lung samples with corresponding Gut samples available
psLung_wGutAvai <- subset_samples(physeqLung9, GutSampleAvailable == "Yes")

# Keep only represented taxa
psLung_wGutAvai_Present <- prune_taxa(taxa_sums(psLung_wGutAvai) > 0, 
                                      psLung_wGutAvai)

# Subset Gut samples with corresponding Lung samples available
psGut_wLungAvai <- subset_samples(physeqGutHum, LungSampleAvailable == "Yes")

# Remove replicates
psGut_wLungAvai2 <-  subset_samples(psGut_wLungAvai, Replicate  == "No")

# Keep only represented taxa
psGut_wLungAvai2_Present <- prune_taxa(taxa_sums(psGut_wLungAvai2) > 0, 
                                       psGut_wLungAvai2)

# Species present in Gut and (Gut and Lung) samples
# in patients with (Gut and Lung) samples available
Gut_GutLung_Species <- get_taxa_unique(psGut_wLungAvai2_Present, "Species")
Venn_area1 <- length(Gut_GutLung_Species)

# Species present in Lung and (Lung and Gut) samples
# in patients with Gut and Lung samples available
Lung_LungGut_Species <- get_taxa_unique(psLung_wGutAvai_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples
# in patients with Gut and Lung samples available
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for all samples with shared taxa
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("raw", "percent"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammGutLung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with No infection and Samples at Inclusion ----

```{r}
# Subset samples
psGutIncl <- subset_samples(psGut_wLungAvai2,
                            POCGroups20200115 == "Control" & Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psGutIncl)[, c(7,8)]

psLungIncl <- subset_samples(psLung_wGutAvai,
                             POCGroups20200115 == "Control" &  Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psLungIncl)[, c(4,15)]

# Keep only represented taxa
psGutIncl_Present <- prune_taxa(taxa_sums(psGutIncl) > 0, 
                                psGutIncl)
psLungIncl_Present <- prune_taxa(taxa_sums(psLungIncl) > 0, 
                                 psLungIncl)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutIncl_Present)
otu_table(psGutIncl_Present)[ , "S2104"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungIncl_Present)
otu_table(psLungIncl_Present)[ , "S3905"] %>% unique %>% na.exclude %>% length

# Species present at Inclusion in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutIncl_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Inclusion in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungIncl_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_NoInf_Inclusion.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with No infection and Samples at Extubation ----

```{r}
# Subset samples
psGutExtub <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Control" & Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psGutExtub)[, c(7,8)]

psLungExtub <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Control" &  Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psLungExtub)[, c(4,15)]

# Keep only represented taxa
psGutExtub_Present <- prune_taxa(taxa_sums(psGutExtub) > 0, 
                                 psGutExtub)
psLungExtub_Present <- prune_taxa(taxa_sums(psLungExtub) > 0, 
                                  psLungExtub)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutExtub_Present)
otu_table(psGutExtub_Present)[ , "S3906"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungExtub_Present)
otu_table(psLungExtub_Present)[ , "S3907"] %>% unique %>% na.exclude %>% length

# Species present at Extubation in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutExtub_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Extubation in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungExtub_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_NoInf_Extubation.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Inclusion ----

```{r}
# Subset samples
psGutIncl <- subset_samples(psGut_wLungAvai2,
                            POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psGutIncl)[, c(7,8)]

psLungIncl <- subset_samples(psLung_wGutAvai,
                             POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psLungIncl)[, c(4,15)]

# Keep only represented taxa
psGutIncl_Present <- prune_taxa(taxa_sums(psGutIncl) > 0, 
                                psGutIncl)
psLungIncl_Present <- prune_taxa(taxa_sums(psLungIncl) > 0, 
                                 psLungIncl)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutIncl_Present)
otu_table(psGutIncl_Present)[ , "S3304"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungIncl_Present)
otu_table(psLungIncl_Present)[ , "S4204"] %>% unique %>% na.exclude %>% length

# Species present at Inclusion in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutIncl_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Inclusion in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungIncl_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_Inclusion.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Infection D1 ----

```{r}
# Subset samples
psGutInfD1 <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Infection D1")
sample_data(psGutInfD1)[, c(7,8)]

psLungInfD1 <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Infection_D1")
sample_data(psLungInfD1)[, c(4,15)]

# Keep only represented taxa
psGutInfD1_Present <- prune_taxa(taxa_sums(psGutInfD1) > 0, 
                                 psGutInfD1)
psLungInfD1_Present <- prune_taxa(taxa_sums(psLungInfD1) > 0, 
                                  psLungInfD1)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutInfD1_Present)
otu_table(psGutInfD1_Present)[ , "S3604"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungInfD1_Present)
otu_table(psLungInfD1_Present)[ , "S4305"] %>% unique %>% na.exclude %>% length

# Species present at Infection D1 in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutInfD1_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Infection D1 in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungInfD1_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_InfD1.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Infection D5 ----

```{r}
# Subset samples
psGutInfD5 <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Infection D5")
sample_data(psGutInfD5)[, c(7,8)]

psLungInfD5 <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Infection_D5")
sample_data(psLungInfD5)[, c(4,15)]

# Keep only represented taxa
psGutInfD5_Present <- prune_taxa(taxa_sums(psGutInfD5) > 0, 
                                 psGutInfD5)
psLungInfD5_Present <- prune_taxa(taxa_sums(psLungInfD5) > 0, 
                                  psLungInfD5)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutInfD5_Present)
otu_table(psGutInfD5_Present)[ , "S3610"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungInfD5_Present)
otu_table(psLungInfD5_Present)[ , "S4112"] %>% unique %>% na.exclude %>% length

# Species present at Infection D5 in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutInfD5_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Infection D5 in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungInfD5_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_InfD5.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Extubation ----

```{r}
# Subset samples
psGutExtub <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psGutExtub)[, c(7,8)]

psLungExtub <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psLungExtub)[, c(4,15)]

# Keep only represented taxa
psGutExtub_Present <- prune_taxa(taxa_sums(psGutExtub) > 0, 
                                 psGutExtub)
psLungExtub_Present <- prune_taxa(taxa_sums(psLungExtub) > 0, 
                                  psLungExtub)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutExtub_Present)
otu_table(psGutExtub_Present)[ , "S3012"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungExtub_Present)
otu_table(psLungExtub_Present)[ , "S3714"] %>% unique %>% na.exclude %>% length

# Species present at Extubation in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutExtub_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Extubation in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungExtub_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_Extub.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Analysis limited to taxa shared between Gut and Lung ----

```{r}
# Work with Objects obtained in Lung-Gut comparisons - Venn diagram - ...
psLung_wGutAvai_Present
psGut_wLungAvai2_Present

# Obtain lists of taxa to subset
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)

psLungShared <- subset_taxa(psLung_wGutAvai_Present, Species %in% SharedSpecies)
psGutShared <- subset_taxa(psGut_wLungAvai2_Present, Species %in% SharedSpecies)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

# Agglomeration
psLungSharedRel_Gen <- tax_glom(psLungSharedRel, "Genus", NArm = TRUE)
psGutSharedRel_Gen <- tax_glom(psGutSharedRel, "Genus", NArm = TRUE)

# Top x genera with shared taxa in the Lung
TopGenNamesLung <-  names(sort(taxa_sums(psLungSharedRel_Gen), decreasing = TRUE)[1:20])
psLungTopGen <-  prune_taxa(TopGenNamesLung, psLungSharedRel_Gen)

# Prepare a df for subsequent plot
psLungTopGen_df <- psmelt(psLungTopGen)
str(psLungTopGen_df)

# Convert as required
psLungTopGen_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(psLungTopGen_df)

# Grouped boxplot
# Reverse Group levels for vertical plot
psLungTopGen_df$Filtered_Data_AK_Dec2020 <- factor(psLungTopGen_df$Filtered_Data_AK_Dec2020,
                                                   levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungTopGen_df %>%
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Genus", y = "Relative abundance") +
  theme_classic()

ggsave("TopGenLung_sharedGutLung_Boxplot.pdf",
       width = 25, height = 30, dpi = 200, units = "cm")

# Top x genera with shared taxa in the Gut
TopGenNamesGut <-  names(sort(taxa_sums(psGutSharedRel_Gen), decreasing = TRUE)[1:20])
psGutTopGen <-  prune_taxa(TopGenNamesGut, psGutSharedRel_Gen)

# Prepare a df for subsequent plot
psGutTopGen_df <- psmelt(psGutTopGen)
str(psGutTopGen_df)

# Convert as required
psGutTopGen_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(psGutTopGen_df)

# Grouped boxplot
# Reverse Group levels for vertical plot
psGutTopGen_df$Filtered_Data_AK_Dec2020 <- factor(psGutTopGen_df$Filtered_Data_AK_Dec2020,
                                                  levels = c("Discharge", "Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutTopGen_df %>%
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Genus", y = "Relative abundance") +
  theme_classic()

ggsave("TopGenGut_sharedGutLung_Boxplot.pdf",
       width = 25, height = 30, dpi = 200, units = "cm")

# Analysis of Prevotella species in Lung and Gut
psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

LungPrevotellaShared <- subset_taxa(psLungSharedRel, Genus == "Prevotella")
GutPrevotellaShared <- subset_taxa(psGutSharedRel, Genus == "Prevotella")

# Agglomeration
psLungSharedRel_Spe <- tax_glom(LungPrevotellaShared, "Species", NArm = TRUE)
psGutSharedRel_Spe <- tax_glom(GutPrevotellaShared, "Species", NArm = TRUE)

# Prepare a df for subsequent plot
psLungSharedRel_Spe_df <- psmelt(psLungSharedRel_Spe)
str(psLungSharedRel_Spe_df)
psGutSharedRel_Spe_df <- psmelt(psGutSharedRel_Spe)
str(psGutSharedRel_Spe_df)

# Convert as required
psLungSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psLungSharedRel_Spe_df)

psGutSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psGutSharedRel_Spe_df)

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                          levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("LungPrevotella_sharedGutLung_Boxplot.pdf",
       width = 20, height = 14, dpi = 200, units = "cm")

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                         levels = c("Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("GutPrevotella_sharedGutLung_Boxplot.pdf",
       width = 20, height = 14, dpi = 200, units = "cm")

# Analysis of Enterobacter species in Lung and Gut
psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

LungEnterobacterShared <- subset_taxa(psLungSharedRel, Genus == "Enterobacter")
GutEnterobacterShared <- subset_taxa(psGutSharedRel, Genus == "Enterobacter")

# Agglomeration
psLungSharedRel_Spe <- tax_glom(LungEnterobacterShared, "Species", NArm = TRUE)
psGutSharedRel_Spe <- tax_glom(GutEnterobacterShared, "Species", NArm = TRUE)

# Prepare a df for subsequent plot
psLungSharedRel_Spe_df <- psmelt(psLungSharedRel_Spe)
str(psLungSharedRel_Spe_df)
psGutSharedRel_Spe_df <- psmelt(psGutSharedRel_Spe)
str(psGutSharedRel_Spe_df)

# Convert as required
psLungSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psLungSharedRel_Spe_df)

psGutSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psGutSharedRel_Spe_df)

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                          levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("LungEnterobacter_sharedGutLung_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                         levels = c("Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("GutPrevotella_sharedGutLung_Boxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Analysis of Enterobacteriaceae species in Lung and Gut
psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

LungEnterobacterShared <- subset_taxa(psLungSharedRel, Family == "Enterobacteriaceae")
GutEnterobacterShared <- subset_taxa(psGutSharedRel, Family == "Enterobacteriaceae")

# Agglomeration
psLungSharedRel_Spe <- tax_glom(LungEnterobacterShared, "Species", NArm = TRUE)
psGutSharedRel_Spe <- tax_glom(GutEnterobacterShared, "Species", NArm = TRUE)

# Prepare a df for subsequent plot
psLungSharedRel_Spe_df <- psmelt(psLungSharedRel_Spe)
str(psLungSharedRel_Spe_df)
psGutSharedRel_Spe_df <- psmelt(psGutSharedRel_Spe)
str(psGutSharedRel_Spe_df)

# Convert as required
psLungSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psLungSharedRel_Spe_df)

psGutSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psGutSharedRel_Spe_df)

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                          levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("LungEnterobacteriaceae_sharedGutLung_Boxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                         levels = c("Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("GutPrevotella_sharedGutLung_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")
```

#################################################################################### 

#################################################################################### 

/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#START OF GUT MICROBIOTA ANALYSIS/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#

#################################################################################### 

#################################################################################### 

/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#START OF GUT MICROBIOTA ANALYSIS/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#

#################################################################################### 

#################################################################################### 

/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#START OF GUT MICROBIOTA ANALYSIS/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#

#### Gut data analysis / Community abundance manipulation and exploration ----

```{r}
physeqGut1

OTU_Gut1 <- otu_table(physeqGut4)
sam_Gut1 <- sample_data(physeqGut4)
tax_Gut1 <- tax_table(physeqGut4)
tree_Gut1 <- phy_tree(physeqGut4)
refseq_Gut1 <- refseq(physeqGut4)

get_taxa_unique(physeqGut4, "Kingdom")
get_taxa_unique(physeqGut4, "Phylum")
get_taxa_unique(physeqGut4, "Order")
get_taxa_unique(physeqGut4, "Family")
get_taxa_unique(physeqGut4, "Genus")
get_taxa_unique(physeqGut4, "Species")

sort(sample_sums(physeqGut4))

# Top100 most abundant ASVs
Top100AbundantASVCountsGut <- sort(taxa_sums(physeqGut4), decreasing=TRUE)[1:100]
Top100AbundantASVNamesGut  <-  names(Top100AbundantASVCountsGut)
physeGut4Top100  <-  prune_taxa(Top100AbundantASVNamesGut, physeqGut4)
# Visualize taxonomy of most abundant ASVs
Top100AbundantTaxaGut <- tax_table(physeGut4Top100)[ , c("Phylum", "Genus" ,"Species")]

# Number of total ASVs per phylum 
TablePhylaTaxa <- table(tax_table(physeqGut4)[ , 2])
sort(TablePhylaTaxa, decreasing=TRUE)

# Set a palette for Phyla with colours visible to colour-blind people
display.brewer.all(colorblindFriendly = TRUE)
brewer.pal(n = 12, name = "Paired")
brewer.pal(n = 8, name = "Dark2")
PhylumPalGut13 <- c("Firmicutes"="#33A02C", "Fusobacteria"="#CAB2D6",
                    "Actinobacteria"="#E31A1C", "Bacteroidetes"="#FF7F00",
                    "Proteobacteria"="#1F78B4", "Saccharibacteria_TM7"="#FB9A99",
                    "Spirochaetes"="#FDBF6F", "Synergistetes"="#A6CEE3",
                    "Tenericutes"="#B2DF8A", "Cyanobacteria"="#FFFF99",
                    "Lentisphaerae"="#6A3D9A", "Verrucomicrobia"="#B15928",
                    "Deinococcus-Thermus"="#666666")

# Transform counts to relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqGut4Rel <- transform_sample_counts(physeqGut4, count_to_rel_abund)
# Check total = 1 for the first 10 samples
sample_sums(physeqGut4Rel)[1:10]
```

#### Composition analysis - Phylum-level ----

```{r}
# Focus on human samples
get_variable(physeqGut4, varName = "Filtered_Data_AK_Dec2020")
physeqGutHum <- subset_samples(physeqGut4, !Filtered_Data_AK_Dec2020 == "Neg control")

# Transform counts to relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqGutHumRel <- transform_sample_counts(physeqGutHum, count_to_rel_abund)

# Phylum agglomeration
physeqGutHumRelPhy <- tax_glom(physeqGutHumRel, "Phylum", NArm = TRUE) # Agglomeration

physeqGutHumRelPhydf <- psmelt(physeqGutHumRelPhy) # Obtain a data frame from phyloseq object.
str(physeqGutHumRelPhydf)

physeqGutHumRelPhydf%>% 
  select(OTU, Sample, Abundance, Phylum, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(OTU, Sample, Phylum, Filtered_Data_AK_Dec2020)) %>% arrange(-(Abundance))

PhyDescBoxPlot<-physeqGutHumRelPhydf %>% 
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill = Phylum), alpha = 0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Phylum" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

print(PhyDescBoxPlot)

ggsave("Phylum_TotalSamples_Gut_Boxplot.pdf",
       width = 15, height = 12, dpi = 400, units = "cm")

InclInfPhyDescBoxPlot <- PhyDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("Phylum_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 30, height = 8, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
physeqGutHumRelPhydf$Filtered_Data_AK_Dec2020 <- factor(physeqGutHumRelPhydf$Filtered_Data_AK_Dec2020,
                                                        levels = c("Discharge", "Extubation",
                                                                   "Infection D5", "Infection D1", "Inclusion"))

physeqGutHumRelPhydf %>%
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("Phylum_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 15, height = 18, dpi = 200, units = "cm")

# Zoom on Phyla families 4 to 10
Top4to10Phy <-  sort(tapply(taxa_sums(physeqGutHumRelPhy),
                            tax_table(physeqGutHumRelPhy)[, "Phylum"], sum), decreasing = TRUE)[4:10]

# Top 4 to 10 ordered phylum names
Top4to10ListForGraph <- as.data.frame(Top4to10Phy) %>% rownames()
# Reverse order for vertical graph
Top4to10ListForGraphRev <- rev(Top4to10ListForGraph)

physeqGutHumRelPhydf %>%
  ggplot(aes(x = Phylum, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim (Top4to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.25)) +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("Phylum4to10_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")



# Zoom on Phyla families 1 to 3
Top3 <-  sort(tapply(taxa_sums(physeqGutHumRelPhy),
                     tax_table(physeqGutHumRelPhy)[, "Phylum"], sum), decreasing = TRUE)[1:3]

# Top 3 ordered phylum names
Top3ListForGraph <- as.data.frame(Top3) %>% rownames()
# Reverse order for vertical graph
Top3ListForGraphRev <- rev(Top3ListForGraph)

physeqGutHumRelPhydf %>%
  ggplot(aes(x = Phylum, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim (Top3ListForGraphRev) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("Phylum_Top3_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 25, height = 20, dpi = 200, units = "cm")


# Plot Heatmap of Phyla abundance distribution according to the Inclusion-Infection Group

# Subset phyloseq object for the infected group
physeqgutinfected <- subset_samples(physeqGutHumRelPhy, !POCGroups20200115=="Control")

# Subset phyloseq object for the non infected group
physeqgutnotinfected<-subset_samples(physeqGutHumRelPhy,POCGroups20200115=="Control")

#Plot Heatmap ad "Phylum" level according to sample timepoint (inclusion-D1-D5-Extubation-Discharge)
plot_heatmap(physeqgutinfected, taxa.label="Phylum")+facet_grid(~Filtered_Data_AK_Dec2020, scales = "free_x")
plot_heatmap(physeqgutnotinfected, taxa.label="Phylum")+facet_grid(~Filtered_Data_AK_Dec2020, scales = "free_x")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
# Construct a highly agglomerated tree to obtain a reference for further comparisons
tax_table(physeqGutHumRelPhy)

# Merge all taxa except one poorly represented phylum
MergedPhy <- merge_taxa(physeqGutHumRelPhy, taxa_names(physeqGutHumRelPhy)[-c(2)], 1)
tax_table(MergedPhy)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeGen <- plot_tree(MergedPhy, label.tips="Phylum",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         plot.margin = 0.1,
                         nodelabf = nodeplotblank,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePhyMergedInclInf_Gut.pdf",
       width = 85, height = 10, dpi = 200, units = "cm", limitsize = FALSE)

# Phylum-level tree
plotTreePhy <- plot_tree(physeqGutHumRelPhy, label.tips="Phylum",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         plot.margin = 0.1,
                         nodelabf = nodeplotblank,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePhyInclInf_Gut.pdf",
       width = 65, height = 15, dpi = 200, units = "cm", limitsize = FALSE)



# Phylum-level tree only for infected patients (physeq object was subsetted before for infected patients)
plot_tree(physeqgutinfected, label.tips="Phylum",
          color = "Filtered_Data_AK_Dec2020",
          size = "Abundance",
          sizebase = 10,
          plot.margin = 0.1,
          nodelabf = nodeplotblank,
          ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)


ggsave("PtreePhyonlyinfectedf_Gut.pdf",
       width = 65, height = 15, dpi = 200, units = "cm", limitsize = FALSE)





# Phylum-level tree only for not infected patients (physeq object was subsetted before for infected patients)
plot_tree(physeqgutnotinfected, label.tips="Phylum",
          color = "Filtered_Data_AK_Dec2020",
          size = "Abundance",
          sizebase = 10,
          plot.margin = 0.1,
          nodelabf = nodeplotblank,
          ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)


ggsave("PtreePhynotinfected_Gut.pdf",
       width = 65, height = 15, dpi = 200, units = "cm", limitsize = FALSE)


```

### Composition analysis - Class-level ----

```{r}
physeqGutHumRelCla <- tax_glom(physeqGutHumRel, "Class", NArm = TRUE) # Agglomeration
ntaxa(physeqGutHumRelCla)

physeqGutHumRelCladf <- psmelt(physeqGutHumRelCla) # Obtain a data frame from phyloseq object.
str(physeqGutHumRelCladf)

physeqGutHumRelCladf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class) %>% 
  convert(fct(OTU, Sample, Phylum, Class)) %>% arrange(-(Abundance))

ClaDescBoxPlot <- physeqGutHumRelCladf %>% 
  ggplot(aes(x = reorder(Class, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Class", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Class_Gut_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
plotTreeCla <- plot_tree(physeqGutHumRelCla, label.tips="Class",
                         color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10,
                         plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeClaInclInf_Gut.pdf",
       width = 55, height = 15, dpi = 200, units = "cm", limitsize = FALSE)

plotTreeClaCirc <- plot_tree(physeqGutHumRelCla, label.tips="Class", 
                             color="Filtered_Data_AK_Dec2020",
                             sizebase=10, ladderize="left", 
                             plot.margin= 0.1) +
  scale_color_manual(values = InclInfectColors5) +
  coord_polar(theta="y")

ggsave("PtreeClaInclInf_Gut_Circ.pdf",
       width = 25, height = 25, dpi = 200, units = "cm", limitsize = FALSE)

```

### Composition analysis - Order-level ----

```{r}
physeqGutHumRelOrd <- tax_glom(physeqGutHumRel, "Order", NArm = TRUE) # Agglomeration
ntaxa(physeqGutHumRelOrd)

physeqGutHumRelOrddf <- psmelt(physeqGutHumRelOrd) # Obtain a data frame from phyloseq object.
str(physeqGutHumRelOrddf)

physeqGutHumRelOrddf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class, Order) %>% 
  convert(fct(OTU, Sample, Phylum, Class, Order)) %>% arrange(-(Abundance))

OrdDescBoxPlot <- physeqGutHumRelOrddf %>% 
  ggplot(aes(x = reorder(Order, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Order", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Order_Gut_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
plotTreeOrd <- plot_tree(physeqGutHumRelOrd, label.tips="Order",
                         color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10,
                         plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeOrdInclInf_Gut.pdf",
       width = 55, height = 25, dpi = 200, units = "cm", limitsize = FALSE)

plotTreeOrdCirc <- plot_tree(physeqGutHumRelOrd, label.tips="Order", 
                             color="Filtered_Data_AK_Dec2020",
                             sizebase=10, ladderize="left", 
                             plot.margin= 0.1) +
  scale_color_manual(values = InclInfectColors5) +
  coord_polar(theta="y")

ggsave("PtreeOrdInclInf_Gut_Circ.pdf",
       width = 25, height = 25, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Family-level - Total samples ----

```{r}
physeqGutHumRelFam <- tax_glom(physeqGutHumRel, "Family", NArm = TRUE) # Agglomeration
ntaxa(physeqGutHumRelFam)

physeqGutHumRelFamdf <- psmelt(physeqGutHumRelFam) # Obtain a data frame from phyloseq object.
str(physeqGutHumRelFamdf)

physeqGutHumRelFamdf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class, Order, Family, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(OTU, Sample, Phylum, Class, Order, Family)) %>% arrange(-(Abundance))

FamDescBoxPlot <- physeqGutHumRelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Family_TotalSamples_Gut_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

FamDescBoxPlotInclInf <- FamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("Family_per_InclInfGroup_Gut_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
physeqGutHumRelFamdf$Filtered_Data_AK_Dec2020 <- factor(physeqGutHumRelFamdf$Filtered_Data_AK_Dec2020,
                                                        levels = c("Discharge","Extubation", "Infection D5", "Infection D1", "Inclusion"))
# Linear scale
physeqGutHumRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

# Log2 scale
physeqGutHumRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  scale_y_continuous(trans = "log2") +
  coord_flip() +
  labs(x = "Phylum", y = "Relative abundance") +
  theme_classic()

ggsave("2Family_TotalSamples_InclInfGroup_Gut_Boxplot.pdf",
       width = 15, height = 45, dpi = 200, units = "cm")

# Plot tree and display prevalence and abundance according to the Inclusion-Infection group
plotTreeFam <- plot_tree(physeqGutHumRelFam, label.tips="Family",
                         color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10,
                         plot.margin = 0.1, ladderize="right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeFamInclInf_Gut.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)

# plot_tree with shape corresponding to No Infect/Pulm Infect/Extra-pulm Infect
plotTreeFam <- plot_tree(physeqGutHumRelFam, nodelabf=nodeplotblank, label.tips="Family",
                         color="Filtered_Data_AK_Dec2020", size="Abundance", sizebase=10,
                         shape="POCGroups20200115", base.spacing=0.04, plot.margin = 0.1,
                         ladderize="right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeFamInclInf_shapeNoInfPulmExtraPulmInf_Gut.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)
```

### Plot tree - Families - Total samples ----

```{r}
# Display prevalence and abundance according to the Inclusion-Infection group
plotTreeFam <- plot_tree(physeqGutHumRelFam, label.tips="Family", color="Filtered_Data_AK_Dec2020", 
                         size="Abundance", sizebase=10, plot.margin = 0.1, ladderize="right")

plotTreeFam + scale_color_manual(values = InclInfectColors5)

ggsave("PtreeFamInclInf.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)

# plot_tree with shape corresponding to No Infect/Pulm Infect/Extra-pulm Infect
plotTreeFam <- plot_tree(physeqGutHumRelFam, nodelabf=nodeplotblank, label.tips="Family",
                         color="Filtered_Data_AK_Dec2020", size="Abundance", sizebase=10,
                         shape="POCGroups20200115", base.spacing=0.04, plot.margin = 0.1, ladderize="right")

plotTreeFam + scale_color_manual(values = c("Inclusion"="gold2", "Infection_D1"="coral2",
                                            "Infection_D5"="chartreuse4",
                                            "Extubation"="cadetblue3"))

ggsave("PtreeFamInclInf_shapeNoInfPulmExtraPulmInf.pdf",
       width = 55, height = 30, dpi = 200, units = "cm", limitsize = FALSE)
```

### Taxa representation across total samples ----

```{r}
physeqGutHumRelFam <- tax_glom(physeqGutHumRel, "Family", NArm = TRUE) # Agglomeration

# Lachnospiraceae
physeqSubset <- subset_taxa(physeqGutHumRel, Family == "Lachnospiraceae")
physeqGlom <- tax_glom(physeqSubset, "Family", NArm = TRUE)
TaxName <- taxa_names(physeqGlom)
physeqGlomPresent <- prune_taxa(taxa_sums(physeqGlom) > 0, physeqGlom)

Lachnospiraceae <- get_sample(physeqGlomPresent, TaxName)
LachnospiraceaePresent <- Lachnospiraceae[Lachnospiraceae>0]
length(LachnospiraceaePresent)
```

### Composition analysis - Family-level - Per Phylum: Firmicutes -----

```{r}
psFirmiRel <- subset_taxa(physeqGutHumRel, Phylum == "Firmicutes", na.rm = TRUE)

psFirmiRelFam <- tax_glom(psFirmiRel, "Family", NArm  = TRUE) # Agglomeration
ntaxa(psFirmiRelFam)

psFirmiRelFamdf <- psmelt(psFirmiRelFam) # Obtain a data frame from phyloseq object
str(psFirmiRelFamdf)

psFirmiRelFamdf %>% 
  select(Phylum, Family, OTU, AKSampleID, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, OTU, AKSampleID)) %>% arrange(-(Abundance))

FirmiFamDescBoxPlot <- psFirmiRelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("FirmiFamily_TotalSamples_Gut_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

FirmiFamDescBoxPlotInclInf <- FirmiFamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("FirmiFamily_per_InclInfGroup_Gut_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
psFirmiRelFamdf$Filtered_Data_AK_Dec2020 <- factor(psFirmiRelFamdf$Filtered_Data_AK_Dec2020,
                                                   levels = c("Discharge","Extubation", "Infection D5", "Infection D1", "Inclusion"))
# Linear scale
psFirmiRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

# Log2 scale
psFirmiRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  scale_y_continuous(trans = "log2") +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("FirmiFamily_TotalSamples_InclInfGroup_log2_Lung_Boxplot.pdf",
       width = 15, height = 25, dpi = 200, units = "cm")


# Limit to the Top 10 abundant Firmicutes families
Top10Fam <-  sort(tapply(taxa_sums(psFirmiRelFam),
                         tax_table(psFirmiRelFam)[, "Family"], sum), decreasing = TRUE)[1:10]

# Top 10 ordered family names
Top10ListForGraph <- as.data.frame(Top10Fam) %>% rownames()
# Reverse order for vertical graph
Top10ListForGraphRev <- rev(Top10ListForGraph)

# Linear scale
psFirmiRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim(Top10ListForGraphRev) +
  ylim(c(0,1.0)) +
  coord_flip() +
  labs(x = "Firmicutes Phylum Families", y = "Relative abundance") +
  theme_classic()

ggsave("FirmiFamily_TotalSamples_InclInfGroup_Gut_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")

# Zoom on Firmicutes families 4 to 10
Top4to10Fam <-  sort(tapply(taxa_sums(psFirmiRelFam),
                            tax_table(psFirmiRelFam)[, "Family"], sum), decreasing = TRUE)[4:10]

# Top 4 to 10 ordered family names
Top4to10ListForGraph <- as.data.frame(Top4to10Fam) %>% rownames()
# Reverse order for vertical graph
Top4to10ListForGraphRev <- rev(Top4to10ListForGraph)

# Linear scale
psFirmiRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim (Top4to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.2)) +
  labs(x = "Firmicutes Phylum Families 4 to 10", y = "Relative abundance") +
  theme_classic()

ggsave("Firmi4to10Family_TotalSamples_InclInfGroup_Gut_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")
```

### Composition analysis - Family-level - Per Phylum: Proteobacteria -----

```{r}
psProteoRel <- subset_taxa(physeqGutHumRel, Phylum == "Proteobacteria", NArm = TRUE)

psProteoRelFam <- tax_glom(psProteoRel, "Family", NArm = TRUE) # Agglomeration
ntaxa(psProteoRelFam)

psProteoRelFamdf <- psmelt(psProteoRelFam) # Obtain a data frame from phyloseq object
str(psProteoRelFamdf)

psProteoRelFamdf %<>% 
  select(Phylum, Family, OTU, AKSampleID, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, OTU, AKSampleID)) %>% arrange(-(Abundance))

ProteoFamDescBoxPlot <- psProteoRelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("ProteoFamily_TotalSamples_Gut_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

ProteoFamDescBoxPlotInclInf <- ProteoFamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("ProteoFamily_per_InclInfGroup_Gut_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
psProteoRelFamdf$Filtered_Data_AK_Dec2020 <- factor(psProteoRelFamdf$Filtered_Data_AK_Dec2020,
                                                    levels = c("Discharge","Extubation", "Infection D5", "Infection D1", "Inclusion"))
# Linear scale - Order families by sum to show those with tipical pathogens
psProteoRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=sum), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("ProteoFamily_TotalSamples_InclInfGroup_Gut_linear_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

# Limit to the Top 10 abundant Proteobacteria families
Top10Fam <-  sort(tapply(taxa_sums(psProteoRelFam),
                         tax_table(psProteoRelFam)[, "Family"], sum), decreasing = TRUE)[1:10]

# Top 10 ordered family names
Top10ListForGraph <- as.data.frame(Top10Fam) %>% rownames()
# Reverse order for vertical graph
Top10ListForGraphRev <- rev(Top10ListForGraph)

# Linear scale
psProteoRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim(Top10ListForGraphRev) +
  ylim(c(0,1.0)) +
  coord_flip() +
  labs(x = "Proteobacteria Phylum Families", y = "Relative abundance") +
  theme_classic()

ggsave("ProteoFamily1to2_TotalSamples_InclInfGroup_Gut_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")

# Zoom on Proteobacteria families 3 to 10
Top3to10Fam <-  sort(tapply(taxa_sums(psProteoRelFam),
                            tax_table(psProteoRelFam)[, "Family"], sum), decreasing = TRUE)[3:10]

# Top 3 to 10 ordered family names
Top3to10ListForGraph <- as.data.frame(Top3to10Fam) %>% rownames()
# Reverse order for vertical graph
Top3to10ListForGraphRev <- rev(Top3to10ListForGraph)

# Linear scale
psProteoRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim (Top3to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.2)) +
  labs(x = "Proteobacteria Phylum Families 3 to 10", y = "Relative abundance") +
  theme_classic()

ggsave("ProteoFamily3to10_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")
```

### Composition analysis - Family-level - Per Phylum: Bacteroidetes -----

```{r}
psBactRel <- subset_taxa(physeqGutHumRel, Phylum == "Bacteroidetes", NArm = TRUE)

psBactRelFam <- tax_glom(psBactRel, "Family", NArm = TRUE) # Agglomeration
ntaxa(psBactRelFam)

psBactRelFamdf <- psmelt(psBactRelFam) # Obtain a data frame from phyloseq object
str(psBactRelFamdf)

psBactRelFamdf %<>% 
  select(Phylum, Family, OTU, AKSampleID, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, OTU, AKSampleID)) %>% arrange(-(Abundance))

BactFamDescBoxPlot <- psBactRelFamdf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Family", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("BactFamily_TotalSamples_Lung_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

BactFamDescBoxPlotInclInf <- BactFamDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("BactFamily_per_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Grouped boxplot
# Reverse Group levels for vertical plot
psBactRelFamdf$Filtered_Data_AK_Dec2020 <- factor(psBactRelFamdf$Filtered_Data_AK_Dec2020,
                                                  levels = c("Discharge","Extubation", "Infection D5", "Infection D1", "Inclusion"))
# Linear scale - Order families by sum to show those with tipical pathogens
psBactRelFamdf %>%
  ggplot(aes(x = reorder(Family, Abundance, FUN=sum), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Family", y = "Relative abundance") +
  theme_classic()

ggsave("BactFamily_TotalSamples_InclInfGroup_Gut_linear_Boxplot.pdf",
       width = 25, height = 25, dpi = 200, units = "cm")

# Limit to the Top 10 abundant Bacteroidetes families
Top10Fam <-  sort(tapply(taxa_sums(psBactRelFam),
                         tax_table(psBactRelFam)[, "Family"], sum), decreasing = TRUE)[1:10]

# Top 10 ordered family names
Top10ListForGraph <- as.data.frame(Top10Fam) %>% rownames()
# Reverse order for vertical graph
Top10ListForGraphRev <- rev(Top10ListForGraph)

# Linear scale
psBactRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim(Top10ListForGraphRev) +
  ylim(c(0,1.0)) +
  coord_flip() +
  labs(x = "Bacteroidetes Phylum Families", y = "Relative abundance") +
  theme_classic()

ggsave("BactFamily1to2_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")

# Zoom on Bacteroidetes families 4 to 10
Top4to10Fam <-  sort(tapply(taxa_sums(psBactRelFam),
                            tax_table(psBactRelFam)[, "Family"], sum), decreasing = TRUE)[4:10]

# Top 4 to 10 ordered family names
Top4to10ListForGraph <- as.data.frame(Top4to10Fam) %>% rownames()
# Reverse order for vertical graph
Top4to10ListForGraphRev <- rev(Top4to10ListForGraph)

# Linear scale
psBactRelFamdf %>%
  ggplot(aes(x = Family, y = Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  xlim (Top4to10ListForGraphRev) +
  coord_flip(ylim = c(0, 0.1)) +
  labs(x = "Bacteroidetes Phylum top 4-10 Families", y = "Relative abundance") +
  theme_classic()

ggsave("BactFamily3to10_TotalSamples_InclInfGroup_Lung_linear_Boxplot.pdf",
       width = 20, height = 12, dpi = 200, units = "cm")
```

### Composition analysis - Genus-level - Total samples ----

```{r}
physeqGutHumRelGen <- tax_glom(physeqGutHumRel, "Genus", NArm = TRUE) # Agglomeration
ntaxa(physeqGutHumRelGen)

physeqGutHumRelGendf <- psmelt(physeqGutHumRelGen) # Obtain a data frame from phyloseq object.
str(physeqGutHumRelGendf)

physeqGutHumRelGendf %<>% 
  select(OTU, Sample, Abundance, Phylum, Class, Order, Family, Genus, Filtered_Data_AK_Dec2020) %>%
  convert(fct(OTU, Sample, Phylum, Class, Order, Family, Genus)) %>% arrange(-(Abundance))

GenDescBoxPlot <- physeqGutHumRelGendf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum), alpha=0.7) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Genus_TotalSamples_Gut_Boxplot.pdf",
       width = 25, height = 45, dpi = 200, units = "cm")

GenDescBoxPlotInclInf <- GenDescBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("Genus_per_InclInfGroup_Gut_Boxplot.pdf",
       width = 25, height = 45, dpi = 200, units = "cm")

# plot_tree with shape corresponding to No Infect/Pulm Infect/Extra-pulm Infect
plotTreeGen <- plot_tree(physeqGutHumRelGen, nodelabf=nodeplotblank, label.tips="Genus",
                         color="Filtered_Data_AK_Dec2020", size="Abundance", sizebase=10,
                         shape="POCGroups20200115", base.spacing=0.2, plot.margin = 0.2,
                         text.size=5, ladderize="right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeGenInclInf_shapeNoInfPulmExtraPulmInf_Gut.pdf",
       width = 55, height = 80, dpi = 200, units = "cm", limitsize = FALSE)


plotTreeGenCirc <- plot_tree(physeqGutHumRelGen, label.tips="Genus", 
                             color="Filtered_Data_AK_Dec2020", sizebase=2, ladderize="left", 
                             plot.margin= 0) + coord_polar(theta="y")

ggsave("PtreeGenInclInf_Gut_Circ.pdf",
       width = 40, height = 40, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Species-level - Per Family: Streptococcaceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Streptococcaceae")

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeGen <- plot_tree(physeqSub, label.tips="Species",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         nodelabf = nodeplotblank,
                         plot.margin = 0.1,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeStreptoSpecInclInf_Gut.pdf",
       width = 35, height = 90, dpi = 200, units = "cm", limitsize = FALSE)

### Composition analysis - Genus-level - Per Family: Lachnospiraceae ----
physeqSub <- subset_taxa(physeqGutHumRelGen, Family == "Lachnospiraceae")

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeGen <- plot_tree(physeqSub, label.tips="Genus",
                         color = "Filtered_Data_AK_Dec2020",
                         size = "Abundance",
                         sizebase = 10,
                         plot.margin = 0.1,
                         nodelabf = nodeplotblank,
                         ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeLachnoGenInclInf_Gut.pdf",
       width = 37, height = 6, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Species-level - Per Family: Staphylococcaceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Staphylococcaceae")
tax_table(physeqSub)
physeqSub <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSub)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSub, label.tips="Species",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeStaphyloSpecInclInf_Lung.pdf",
       width = 40, height = 6, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Species-level - Per Family: Pasteurellaceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Pasteurellaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

physeqSubHaemo <- subset_taxa(physeqSubSpe, Genus == "Haemophilus")

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSubHaemo, label.tips="Species",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors)

ggsave("PtreeHaemophilusInclInf_Lung.pdf",
       width = 43, height = 7, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Genus-level - Per Family: Enterobacteraceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Enterobacteriaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSubSpe, label.tips="Species",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreeEnterobactSpecInclInf_Gut.pdf",
       width = 40, height = 6, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Genus-level - Per Family: Prevotelleceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Prevotellaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec <- plot_tree(physeqSubSpe, label.tips="Species",
                          color = "Filtered_Data_AK_Dec2020",
                          size = "Abundance",
                          sizebase = 10,
                          plot.margin = 0.1,
                          nodelabf = nodeplotblank,
                          ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePrevotSpecInclInf_Lung.pdf",
       width = 45, height = 65, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Genus-level - Per Family: Bacteroidaceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Bacteroidaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec 
plot_tree(physeqSubSpe, label.tips="Species",
          color = "Filtered_Data_AK_Dec2020",
          size = "Abundance",
          sizebase = 10,
          plot.margin = 0.1,
          nodelabf = nodeplotblank,
          ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePrevotSpecInclInf_Lung.pdf",
       width = 45, height = 65, dpi = 200, units = "cm", limitsize = FALSE)

```

### Composition analysis - Genus-level - Per Family: Peptoniphilaceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Peptoniphilaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec<-plot_tree(physeqSubSpe, label.tips="Species",
                        color = "Filtered_Data_AK_Dec2020",
                        size = "Abundance",
                        sizebase = 10,
                        plot.margin = 0.1,
                        nodelabf = nodeplotblank,
                        ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePrevotSpecInclInf_Lung.pdf",
       width = 45, height = 65, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Genus-level - Per Family: Veillonellaceae ----

```{r}
physeqSub <- subset_taxa(physeqGutHumRel, Family == "Veillonellaceae")
tax_table(physeqSub)
physeqSubSpe <- tax_glom(physeqSub, "Species", NArm = TRUE)
tax_table(physeqSubSpe)

# Plot tree and display prevalence and abundance according
# to the Inclusion-Infection group
plotTreeSpec<-plot_tree(physeqSubSpe, label.tips="Species",
                        color = "Filtered_Data_AK_Dec2020",
                        size = "Abundance",
                        sizebase = 10,
                        plot.margin = 0.1,
                        nodelabf = nodeplotblank,
                        ladderize = "right") +
  scale_color_manual(values = InclInfectColors5)

ggsave("PtreePrevotSpecInclInf_Gut.pdf",
       width = 45, height = 65, dpi = 200, units = "cm", limitsize = FALSE)
```

### Composition analysis - Patients with Pulm. infection AND (InfD1 + InfD5) available ----

```{r}
PatPulmInf_InfD1InfD5 <- c("1","2","4","8","9","11","12","26","28","30","33","36","37","40","41","42","43")

physPulmInf_InfD1InfD5 <- subset_samples(physeqGutHum, record_ID %in% PatPulmInf_InfD1InfD5)

# Focus on Proteobacteria genera
InfD1InfD5Prot <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Proteobacteria")

InfD1InfD5ProtGen <- tax_glom(InfD1InfD5Prot, "Genus", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5ProtPresent <- prune_taxa(taxa_sums(InfD1InfD5ProtGen) > 0, 
                                    InfD1InfD5ProtGen)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5ProtGenRel <- transform_sample_counts(InfD1InfD5ProtPresent, count_to_rel_abund)

InfD1InfD5ProtGenReldf <- psmelt(InfD1InfD5ProtGenRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5ProtGenReldf)

InfD1InfD5ProtGenReldf %<>% 
  select(Sample, Phylum, Genus, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Genus, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5ProtGenReldf)
tail(InfD1InfD5ProtGenReldf)

InfD1InfD5ProtGenBoxPlot <- InfD1InfD5ProtGenReldf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#1F78B4", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_ProteoGenera_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5ProtGenBoxPlot_perInfGroup <- InfD1InfD5ProtGenBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_ProteoGenera_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 15, dpi = 200, units = "cm")

# Focus on Firmicutes families
InfD1InfD5Firm <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Firmicutes")

InfD1InfD5FirmFam <- tax_glom(InfD1InfD5Firm, "Family", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5FirmFamPresent <- prune_taxa(taxa_sums(InfD1InfD5FirmFam) > 0, 
                                       InfD1InfD5FirmFam)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5FirmFamRel <- transform_sample_counts(InfD1InfD5FirmFamPresent, count_to_rel_abund)

InfD1InfD5FirmFamReldf <- psmelt(InfD1InfD5FirmFamRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5FirmFamReldf)

InfD1InfD5FirmFamReldf %<>% 
  select(Sample, Phylum, Family, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Family, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5FirmFamReldf)
tail(InfD1InfD5FirmFamReldf)

InfD1InfD5FirmFamBoxPlot <- InfD1InfD5FirmFamReldf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#33A02C", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Family" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_FirmFamily_PulmInf_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5FirmFamBoxPlot_perInfGroup <- InfD1InfD5FirmFamBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_FirmFamilies_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")

# Focus on Firmicutes genera
InfD1InfD5Firm <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Firmicutes")

InfD1InfD5FirmGen <- tax_glom(InfD1InfD5Firm, "Genus", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5FirmGenPresent <- prune_taxa(taxa_sums(InfD1InfD5FirmGen) > 0, 
                                       InfD1InfD5FirmGen)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5FirmGenRel <- transform_sample_counts(InfD1InfD5FirmGenPresent, count_to_rel_abund)

InfD1InfD5FirmGenReldf <- psmelt(InfD1InfD5FirmGenRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5FirmGenReldf)

InfD1InfD5FirmGenReldf %<>% 
  select(Sample, Phylum, Family, Genus, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Family, Genus, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5FirmGenReldf)
tail(InfD1InfD5FirmGenReldf)

InfD1InfD5FirmGenBoxPlot <- InfD1InfD5FirmGenReldf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#33A02C", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_FirmGenera_PulmInf_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5FirmGenBoxPlot_perInfGroup <- InfD1InfD5FirmGenBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_FirmFamilies_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")



# Focus on Bacteroidetes families
InfD1InfD5Bact <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Bacteroidetes")

InfD1InfD5BactFam <- tax_glom(InfD1InfD5Bact, "Family", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5BactFamPresent <- prune_taxa(taxa_sums(InfD1InfD5BactFam) > 0, 
                                       InfD1InfD5BactFam)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5BactFamRel <- transform_sample_counts(InfD1InfD5BactFamPresent, count_to_rel_abund)

InfD1InfD5BactFamReldf <- psmelt(InfD1InfD5BactFamRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5BactFamReldf)

InfD1InfD5BactFamReldf %<>% 
  select(Sample, Phylum, Family, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Family, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5BactFamReldf)
tail(InfD1InfD5BactFamReldf)

InfD1InfD5BactFamBoxPlot <- InfD1InfD5BactFamReldf %>% 
  ggplot(aes(x = reorder(Family, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#FF7F00", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Family" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_BactFamily_PulmInf_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5BactFamBoxPlot_perInfGroup <- InfD1InfD5BactFamBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_BactFamilies_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 7, dpi = 200, units = "cm")

# Focus on Bacteroidetes genera
InfD1InfD5Bact <- subset_taxa(physPulmInf_InfD1InfD5, Phylum == "Bacteroidetes")

InfD1InfD5BactGen <- tax_glom(InfD1InfD5Bact, "Genus", NArm = TRUE)

# Keep only represented taxa
InfD1InfD5BactPresent <- prune_taxa(taxa_sums(InfD1InfD5BactGen) > 0, 
                                    InfD1InfD5BactGen)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

InfD1InfD5BactGenRel <- transform_sample_counts(InfD1InfD5BactPresent, count_to_rel_abund)

InfD1InfD5BactGenReldf <- psmelt(InfD1InfD5BactGenRel) # Obtain a data frame from phyloseq object.
str(InfD1InfD5BactGenReldf)

InfD1InfD5BactGenReldf %<>% 
  select(Sample, Phylum, Genus, OTU, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Sample, Phylum, Genus, OTU)) %>% arrange(-(Abundance)) %>% filter(Abundance>0)

head(InfD1InfD5BactGenReldf)
tail(InfD1InfD5BactGenReldf)

InfD1InfD5BactGenBoxPlot <- InfD1InfD5BactGenReldf %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(fill = "#FF7F00", alpha = 0.7) +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("PulmInf_BactGenera_InclInfGroups_Lung_Boxplot.pdf",
       width = 15, height = 8, dpi = 200, units = "cm")

InfD1InfD5BactGenBoxPlot_perInfGroup <- InfD1InfD5BactGenBoxPlot + facet_wrap(~ Filtered_Data_AK_Dec2020, nrow = 1)

ggsave("PulmInf_BactGenera_InclInfGroup_Lung_Boxplot.pdf",
       width = 25, height = 15, dpi = 200, units = "cm")


### Inclusion group: Phylum-level analysis ----
Inclusion <- subset_samples(physeqGutHumRel, Filtered_Data_AK_Dec2020=="Inclusion")
InclusionPhy <- tax_glom(Inclusion, "Phylum", NArm = TRUE) # Agglomeration

InclusionPhydf <- psmelt(InclusionPhy) # Obtain a data frame from phyloseq object.
str(InclusionPhydf)

InclusionPhydf %<>% 
  select(OTU, Sample, Abundance, Phylum) %>% 
  convert(fct(OTU, Sample, Phylum)) %>% arrange(-(Abundance))

InclPhyDescBoxPlot <- InclusionPhydf %>% 
  ggplot(aes(x = reorder(Phylum, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(aes(fill=Phylum, alpha=0.9)) +
  scale_fill_manual(values = PhylumPalGut13) +
  ylim(c(0, 1.0)) +
  labs(x = "Phyla", y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("Inclusion_Phyla_Gut_Boxplot.pdf",
       width = 25, height = 12, dpi = 200, units = "cm")
```

### Phyloseq heatmaps ----

```{r}
# For plot_heatmap analysis, the sample and taxa names in secondary data, 
# OTU table and taxa table must be changed.
# Include a column in sec data with PatientID_Sample_InclInfGr for 
# labeling the x axis of heatmap
physeqGutHumDat_df <- data.frame(sample_data(physeqGutHum))
head(physeqGutHumDat_df)
rownames(physeqGutHumDat_df)

ForUnXlabeldf <- physeqGutHumDat_df %>% 
  unite(UnXlabel, SampleID, AKSampleID, record_ID, Filtered_Data_AK_Dec2020, sep = "_")

physeqGutHumDat_df$UnXlabel <- ForUnXlabeldf$UnXlabel
head(physeqGutHumDat_df)

# For labeling the y axis of heatmap: 
# Attach taxa names to OTU names in taxa table
physeqGutHumTax_df <- data.frame(tax_table(physeqGutHum))
head(physeqGutHumTax_df)
rownames(physeqGutHumTax_df)
physeqGutHumTax2_df <- physeqGutHumTax_df %>% rownames_to_column("SpID")
head(physeqGutHumTax2_df)

ForUnYlabeldf <- physeqGutHumTax2_df %>% 
  unite(UnYlabel, SpID, Phylum, Family, Genus, Species, sep = "_") %>% 
  arrange(UnYlabel)
head(ForUnYlabeldf)

physeqGutHumTax2Sor_df <- physeqGutHumTax2_df %>% arrange(SpID)
head(physeqGutHumTax2Sor_df)

physeqGutHumTax2Sor_df$UnYlabel <- ForUnYlabeldf$UnYlabel
head(physeqGutHumTax2Sor_df)

rownames(physeqGutHumTax2Sor_df) <- physeqGutHumTax2Sor_df$UnYlabel
head(physeqGutHumTax2Sor_df)

# Change rownames in OTU table to match taxa table
physeqGutHumOTU_df <- data.frame(otu_table(physeqGutHum))
physeqGutHumOTU2Sor_df <- physeqGutHumOTU_df %>% 
  rownames_to_column("SpID") %>% 
  arrange(SpID)
head(physeqGutHumOTU2Sor_df)
dim(physeqGutHumOTU2Sor_df)
dim(physeqGutHumTax2Sor_df)

rownames(physeqGutHumOTU2Sor_df) <- rownames(physeqGutHumTax2Sor_df)

# Remove SpID column to obtain further below a matrix
physeqGutHumOTU2Sor_df %<>% select(-1)
str(physeqGutHumOTU2Sor_df)

# Prepare phyloseq object for heatmap plotting
physeqGutHumUnDat <- sample_data(physeqGutHumDat_df)
class(physeqGutHumUnDat)
str(physeqGutHumUnDat)
sample_names(physeqGutHumUnDat)

physeqGutHumOTU2 <- otu_table(physeqGutHumOTU2Sor_df, taxa_are_rows = TRUE)
class(physeqGutHumOTU2)
sample_names(physeqGutHumOTU2)

physeqGutHumTax2 <- tax_table(as.matrix(physeqGutHumTax2Sor_df))
class(physeqGutHumTax2)

physeqGutHumUn <- phyloseq(physeqGutHumUnDat, physeqGutHumOTU2, physeqGutHumTax2)


# Plot for figure - 2021.01.21
PHeat1_20210121 <- plot_heatmap(physeqGutHumUn, "NMDS", "bray", sample.label="UnXlabel", 
                                taxa.label="UnYlabel", max.label = 1989)

ggsave(filename = "PHeat1_20210121_high_tax_names.pdf",
       width = 40, height = 450, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

ggsave(filename = "PHeat1_20210121.pdf",
       width = 25, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

PHeat1_20210121[["plot_env"]][["sample.order"]]
PHeat1_20210121[["plot_env"]][["taxa.order"]]
tail(PHeat1_20210121[["plot_env"]][["taxa.order"]])

PHeat2 <- plot_heatmap(physeqGutHumUn, "NMDS", "jaccard", sample.label="UnXlabel", 
                       species.label="Genus", max.label = 1989)

ggsave(filename = "PHeat1_20210121.pdf",
       width = 25, height = 600, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

PHeat1facetFiltered_Data_AK_Dec2020 <- PHeat1_20210121 + facet_grid(~Filtered_Data_AK_Dec2020,
                                                                    scales = "free_x")

plot(PHeat1facetFiltered_Data_AK_Dec2020)

ggsave(filename = "PHeat1facetFiltered_Data_AK_Dec2020.pdf",
       width = 25, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

PHeat1facetPatient <- PHeat1_20210121 + facet_grid(~record_ID, scales = "free_x")

plot(PHeat1facetPatient)

ggsave(filename = "PHeat1facetPatient.pdf",
       width = 50, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

```

### Venn diagram - Total samples at Inclusion versus Infection_D5 ----

```{r}
# Get taxa names at Inclusion
sample_data(physeqGutHum)
physeqGutHumIncl <- subset_samples(physeqGutHum, Filtered_Data_AK_Dec2020 %in% 
                                     "Inclusion")
ntaxa(physeqGutHumIncl)

# Keep only represented taxa
physeqGutHumInclPresent <- prune_taxa(taxa_sums(physeqGutHumIncl) > 0, 
                                      physeqGutHumIncl)

ntaxa(physeqGutHumInclPresent)

physeqGutHumInclPresent_df <- data.frame(tax_table(physeqGutHumInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(physeqGutHumInclPresent_df) # enter in area 1 in venn diagram
head(physeqGutHumInclPresent_df)

# Get taxa names at Infection_D5
physeqGutHumInfD5 <- subset_samples(physeqGutHum, Filtered_Data_AK_Dec2020 %in% 
                                      "Infection D5")
ntaxa(physeqGutHumInfD5)

# Keep only represented taxa
physeqGutHumInfD5Present <- prune_taxa(taxa_sums(physeqGutHumInfD5) > 0, 
                                       physeqGutHumInfD5)

ntaxa(physeqGutHumInfD5Present)

physeqGutHumInfD5Present_df <- data.frame(tax_table(physeqGutHumInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(physeqGutHumInfD5Present_df) # enter in area 2 in venn diagram
head(physeqGutHumInfD5Present_df)

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physeqGutHumInclPresent_df, 
                                 physeqGutHumInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
head(taxInclusionOnly_df)

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physeqGutHumInfD5Present_df, 
                             physeqGutHumInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
head(taxInfD5Only_df)

# Taxa present at Inclusion and Infection_D5
taxIncl_AND_InfD5 <- intersect(physeqGutHumInclPresent_df$ASV, physeqGutHumInfD5Present_df$ASV)
length(taxIncl_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create Venn diagram
VennDiagram::draw.pairwise.venn(area1 = 1012,area2 = 544, cross.area = 247,
                                fill = c("Inclusion"="gold2", "Infection D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Inclusion versus Infection_D1 ----

```{r}
Pat3 <- subset_samples(physeqLung9Un, record_ID %in% "3")

# Get taxa names at Inclusion
sample_data(Pat3)
Pat3Incl <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                             "Inclusion")
ntaxa(Pat3Incl)

# Keep only represented taxa
Pat3InclPresent <- prune_taxa(taxa_sums(Pat3Incl) > 0, 
                              Pat3Incl)
# n taxa at Inclusion + Inclusion AND Infection_D1
ntaxa(Pat3InclPresent)

Pat3InclPresent_df <- data.frame(tax_table(Pat3InclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InclPresent_df) # enter in area 1 in venn diagram
head(Pat3InclPresent_df)

# Get taxa names at Infection_D1
Pat3InfD1 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D1")
ntaxa(Pat3InfD1)

# Keep only represented taxa
Pat3InfD1Present <- prune_taxa(taxa_sums(Pat3InfD1) > 0, 
                               Pat3InfD1)
# n taxa at Infection_D1 + Inclusion AND Infection_D1
ntaxa(Pat3InfD1Present)

Pat3InfD1Present_df <- data.frame(tax_table(Pat3InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD1Present_df) # enter in area 2 in venn diagram
head(Pat3InfD1Present_df)

# Get taxa names present at Inclusion ONLY
Pat3taxInclusionOnly_df <- anti_join(Pat3InclPresent_df, 
                                     Pat3InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInclusionOnly_df)
head(Pat3taxInclusionOnly_df)

# Get taxa names present at Infection_D1 ONLY
Pat3taxInfD1Only_df <- anti_join(Pat3InfD1Present_df, 
                                 Pat3InclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD1Only_df)
head(Pat3taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
Pat3taxIncl_AND_InfD1 <- intersect(Pat3InclPresent_df$ASV, Pat3InfD1Present_df$ASV)
length(Pat3taxIncl_AND_InfD1) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 152,area2 = 144, cross.area = 94,
                                fill = c("Inclusion"="gold2", "Infection_D1"="coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D1"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD1_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Infection_D1 vs Infection_D5 ----

```{r}
# Get taxa names at Infection_D1
sample_data(Pat3)
Pat3InfD1 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D1")
ntaxa(Pat3InfD1)

# Keep only represented taxa
Pat3InfD1Present <- prune_taxa(taxa_sums(Pat3InfD1) > 0, 
                               Pat3InfD1)
# n taxa at Infection_D1 + Infection_D1 AND Infection_D5
ntaxa(Pat3InfD1Present)

Pat3InfD1Present_df <- data.frame(tax_table(Pat3InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD1Present_df) # enter in area 1 in venn diagram
head(Pat3InfD1Present_df)

# Get taxa names at Infection_D5
Pat3InfD5 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D5")
ntaxa(Pat3InfD5)

# Keep only represented taxa
Pat3InfD5Present <- prune_taxa(taxa_sums(Pat3InfD5) > 0, 
                               Pat3InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Infection_D1
ntaxa(Pat3InfD5Present)

Pat3InfD5Present_df <- data.frame(tax_table(Pat3InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD5Present_df) # enter in area 2 in venn diagram
head(Pat3InfD5Present_df)

# Get taxa names present at Infection_D1 ONLY
Pat3taxInfD1Only_df <- anti_join(Pat3InfD1Present_df, 
                                 Pat3InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD1Only_df)
head(Pat3taxInfD1Only_df)

# Get taxa names present at Infection_D5 ONLY
Pat3taxInfD5Only_df <- anti_join(Pat3InfD5Present_df, 
                                 Pat3InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD5Only_df)
head(Pat3taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
Pat3taxInfD1_AND_InfD5 <- intersect(Pat3InfD1Present_df$ASV, Pat3InfD5Present_df$ASV)
length(Pat3taxInfD1_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 144,area2 = 46, cross.area = 25,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Infection_D5 vs Extubation ----

```{r}
# Get taxa names at Infection_D5
sample_data(Pat3)
Pat3InfD5 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D5")
ntaxa(Pat3InfD5)

# Keep only represented taxa
Pat3InfD5Present <- prune_taxa(taxa_sums(Pat3InfD5) > 0, 
                               Pat3InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Extubation
ntaxa(Pat3InfD5Present)

Pat3InfD5Present_df <- data.frame(tax_table(Pat3InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD5Present_df) # enter in area 1 in venn diagram
head(Pat3InfD5Present_df)

# Get taxa names at Extubation
Pat3Extub <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Extubation")
ntaxa(Pat3Extub)

# Keep only represented taxa
Pat3ExtubPresent <- prune_taxa(taxa_sums(Pat3Extub) > 0, 
                               Pat3Extub)
# n taxa at Extubation + Extubation AND Infection_D5
ntaxa(Pat3ExtubPresent)

Pat3ExtubPresent_df <- data.frame(tax_table(Pat3ExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3ExtubPresent_df) # enter in area 2 in venn diagram
head(Pat3ExtubPresent_df)

# Get taxa names present at Infection_D5 ONLY
Pat3taxInfD5Only_df <- anti_join(Pat3InfD5Present_df, 
                                 Pat3ExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD5Only_df)
head(Pat3taxInfD5Only_df)

# Get taxa names present at Extubation ONLY
Pat3taxExtubOnly_df <- anti_join(Pat3ExtubPresent_df, 
                                 Pat3InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxExtubOnly_df)
head(Pat3taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
Pat3taxInfD5_AND_Extub <- intersect(Pat3InfD5Present_df$ASV, Pat3ExtubPresent_df$ASV)
length(Pat3taxInfD5_AND_Extub) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 46,area2 = 45, cross.area = 18,
                                fill = c("Infection_D5"="chartreuse4", "Extubation"="cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D5", "Extubation"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD5Extub_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Extrapulmonary infection Patient 3 - Inclusion vs Infection_D5 ----

```{r}
# Get taxa names at Inclusion
sample_data(Pat3)
Pat3Incl <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                             "Inclusion")
ntaxa(Pat3Incl)

# Keep only represented taxa
Pat3InclPresent <- prune_taxa(taxa_sums(Pat3Incl) > 0, 
                              Pat3Incl)
# n taxa at Inclusion + Inclusion AND Infection_D5
ntaxa(Pat3InclPresent)

Pat3InclPresent_df <- data.frame(tax_table(Pat3InclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InclPresent_df) # enter in area 1 in venn diagram
head(Pat3InclPresent_df)

# Get taxa names at Infection_D5
Pat3InfD5 <- subset_samples(Pat3, Filtered_Data_AK_Dec2020 %in% 
                              "Infection_D5")
ntaxa(Pat3InfD5)

# Keep only represented taxa
Pat3InfD5Present <- prune_taxa(taxa_sums(Pat3InfD5) > 0, 
                               Pat3InfD5)
# n taxa at Infection_D5 + Inclusion AND Infection_D5
ntaxa(Pat3InfD5Present)

Pat3InfD5Present_df <- data.frame(tax_table(Pat3InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat3InfD5Present_df) # enter in area 2 in venn diagram
head(Pat3InfD5Present_df)

# Get taxa names present at Inclusion ONLY
Pat3taxInclusionOnly_df <- anti_join(Pat3InclPresent_df, 
                                     Pat3InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInclusionOnly_df)
head(Pat3taxInclusionOnly_df)

# Get taxa names present at Infection_D5 ONLY
Pat3taxInfD5Only_df <- anti_join(Pat3InfD5Present_df, 
                                 Pat3InclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat3taxInfD5Only_df)
head(Pat3taxInfD5Only_df)

# Taxa present at Inclusion and Infection_D5
Pat3taxIncl_AND_InfD5 <- intersect(Pat3InclPresent_df$ASV, Pat3InfD5Present_df$ASV)
length(Pat3taxIncl_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 152,area2 = 46, cross.area = 20,
                                fill = c("Inclusion"="gold2", "Infection D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Pulmonary infection Patient 28 - Inclusion versus Infection_D1 ----

```{r}
Pat28 <- subset_samples(physeqLung9Un, record_ID %in% "28")

# Get taxa names at Inclusion
sample_data(Pat28)
Pat28Incl <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                              "Inclusion")
ntaxa(Pat28Incl)

# Keep only represented taxa
Pat28InclPresent <- prune_taxa(taxa_sums(Pat28Incl) > 0, 
                               Pat28Incl)
# n taxa at Inclusion + Inclusion AND Infection_D1
ntaxa(Pat28InclPresent)

Pat28InclPresent_df <- data.frame(tax_table(Pat28InclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InclPresent_df) # enter in area 1 in venn diagram
head(Pat28InclPresent_df)

# Get taxa names at Infection_D1
Pat28InfD1 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D1")
ntaxa(Pat28InfD1)

# Keep only represented taxa
Pat28InfD1Present <- prune_taxa(taxa_sums(Pat28InfD1) > 0, 
                                Pat28InfD1)
# n taxa at Infection_D1 + Inclusion AND Infection_D1
ntaxa(Pat28InfD1Present)

Pat28InfD1Present_df <- data.frame(tax_table(Pat28InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD1Present_df) # enter in area 2 in venn diagram
head(Pat28InfD1Present_df)

# Get taxa names present at Inclusion ONLY
Pat28taxInclusionOnly_df <- anti_join(Pat28InclPresent_df, 
                                      Pat28InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInclusionOnly_df)
head(Pat28taxInclusionOnly_df)

# Get taxa names present at Infection_D1 ONLY
Pat28taxInfD1Only_df <- anti_join(Pat28InfD1Present_df, 
                                  Pat28InclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD1Only_df)
head(Pat28taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
Pat28taxIncl_AND_InfD1 <- intersect(Pat28InclPresent_df$ASV, Pat28InfD1Present_df$ASV)
length(Pat28taxIncl_AND_InfD1) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 104,area2 = 39, cross.area = 13,
                                fill = c("Inclusion"="gold2", "Infection_D1"="coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection D1"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInclInfD1_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Pulmonary infection Patient 28 - Infection_D1 versus Infection_D5 ----

```{r}
# Get taxa names at Infection_D1
sample_data(Pat28)
Pat28InfD1 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D1")
ntaxa(Pat28InfD1)

# Keep only represented taxa
Pat28InfD1Present <- prune_taxa(taxa_sums(Pat28InfD1) > 0, 
                                Pat28InfD1)
# n taxa at Infection_D1 + Infection_D1 AND Infection_D5
ntaxa(Pat28InfD1Present)

Pat28InfD1Present_df <- data.frame(tax_table(Pat28InfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD1Present_df) # enter in area 1 in venn diagram
head(Pat28InfD1Present_df)

# Get taxa names at Infection_D5
Pat28InfD5 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D5")
ntaxa(Pat28InfD5)

# Keep only represented taxa
Pat28InfD5Present <- prune_taxa(taxa_sums(Pat28InfD5) > 0, 
                                Pat28InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Infection_D1
ntaxa(Pat28InfD5Present)

Pat28InfD5Present_df <- data.frame(tax_table(Pat28InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD5Present_df) # enter in area 2 in venn diagram
head(Pat28InfD5Present_df)

# Get taxa names present at Infection_D1 ONLY
Pat28taxInfD1Only_df <- anti_join(Pat28InfD1Present_df, 
                                  Pat28InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD1Only_df)
head(Pat28taxInfD1Only_df)

# Get taxa names present at Infection_D5 ONLY
Pat28taxInfD5Only_df <- anti_join(Pat28InfD5Present_df, 
                                  Pat28InfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD5Only_df)
head(Pat28taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
Pat28taxInfD1_AND_InfD5 <- intersect(Pat28InfD1Present_df$ASV, Pat28InfD5Present_df$ASV)
length(Pat28taxInfD1_AND_InfD5) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 39,area2 = 38, cross.area = 8,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Pulmonary infection Patient 28 - Infection_D5 versus Extubation ----

```{r}
# Get taxa names at Infection_D5
sample_data(Pat28)
Pat28InfD5 <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Infection_D5")
ntaxa(Pat28InfD5)

# Keep only represented taxa
Pat28InfD5Present <- prune_taxa(taxa_sums(Pat28InfD5) > 0, 
                                Pat28InfD5)
# n taxa at Infection_D5 + Infection_D5 AND Extubation
ntaxa(Pat28InfD5Present)

Pat28InfD5Present_df <- data.frame(tax_table(Pat28InfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28InfD5Present_df) # enter in area 1 in venn diagram
head(Pat28InfD5Present_df)

# Get taxa names at Extubation
Pat28Extub <- subset_samples(Pat28, Filtered_Data_AK_Dec2020 %in% 
                               "Extubation")
ntaxa(Pat28Extub)

# Keep only represented taxa
Pat28ExtubPresent <- prune_taxa(taxa_sums(Pat28Extub) > 0, 
                                Pat28Extub)
# n taxa at Extubation + Extubation AND Infection_D5
ntaxa(Pat28ExtubPresent)

Pat28ExtubPresent_df <- data.frame(tax_table(Pat28ExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
dim(Pat28ExtubPresent_df) # enter in area 2 in venn diagram
head(Pat28ExtubPresent_df)

# Get taxa names present at Infection_D5 ONLY
Pat28taxInfD5Only_df <- anti_join(Pat28InfD5Present_df, 
                                  Pat28ExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxInfD5Only_df)
head(Pat28taxInfD5Only_df)

# Get taxa names present at Extubation ONLY
Pat28taxExtubOnly_df <- anti_join(Pat28ExtubPresent_df, 
                                  Pat28InfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(Pat28taxExtubOnly_df)
head(Pat28taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
Pat28taxInfD5_AND_Extub <- intersect(Pat28InfD5Present_df$ASV, Pat28ExtubPresent_df$ASV)
length(Pat28taxInfD5_AND_Extub) # enter in cross.area in venn diagram

# Move to new plotting page
grid.newpage()

# Create pairwise venn diagram
VennDiagram::draw.pairwise.venn(area1 = 38,area2 = 24, cross.area = 2,
                                fill = c("Infection_D5"="chartreuse4", "Extubation"="cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D5", "Extubation"),
                                # Numbers
                                cex = 2, fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with No infection AND (Inclusion + Extub) available ----

```{r}
# Compare Inclusion with Extubation
PatNoInf_InclExtub <- c("7","19","20","24","39")

physNoInf_InclExtub <- subset_samples(physeqGutHum, record_ID %in% PatNoInf_InclExtub)

# Keep only represented taxa
physNoInf_InclExtubPresent <- prune_taxa(taxa_sums(physNoInf_InclExtub) > 0, 
                                         physNoInf_InclExtub)

TotalTaxaPresent <- ntaxa(physNoInf_InclExtubPresent)

# Subset taxa present either at Inclusion or (Inclusion AND Extubation)
physNoInfIncl <- subset_samples(physNoInf_InclExtubPresent,
                                Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physNoInfInclPresent <- prune_taxa(taxa_sums(physNoInfIncl) > 0, 
                                   physNoInfIncl)

physNoInfInclPresent_df <- data.frame(tax_table(physNoInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physNoInfInclPresent_df)[1] # for area 1 in venn diagram
head(physNoInfInclPresent_df)

# Subset taxa present either at Extubation or (Extubation AND Inclusion)
physNoInfExtub <- subset_samples(physNoInf_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                   "Extubation")

# Keep only represented taxa
physNoInfExtubPresent <- prune_taxa(taxa_sums(physNoInfExtub) > 0, 
                                    physNoInfExtub)

physNoInfExtubPresent_df <- data.frame(tax_table(physNoInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physNoInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physNoInfInclPresent_df, 
                                 physNoInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physNoInfExtubPresent_df, 
                             physNoInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physNoInfInclPresent_df$ASV, physNoInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                fill = c("gold2", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInclExtubAbs_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Inclusion + Extub) available ----

```{r}
# Compare Inclusion with Extubation
PatPulmInf_InclExtub <- c("16","28")

physPulmInf_InclExtub <- subset_samples(physeqGutHum, record_ID %in% PatPulmInf_InclExtub)

# Keep only represented taxa
physPulmInf_InclExtubPresent <- prune_taxa(taxa_sums(physPulmInf_InclExtub) > 0, 
                                           physPulmInf_InclExtub)

TotalTaxaPresent <- ntaxa(physPulmInf_InclExtubPresent)

# Subset taxa present either in Inclusion or Inclusion AND Extubation
physPulmInfIncl <- subset_samples(physPulmInf_InclExtubPresent,
                                  Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physPulmInfInclPresent <- prune_taxa(taxa_sums(physPulmInfIncl) > 0, 
                                     physPulmInfIncl)

physPulmInfInclPresent_df <- data.frame(tax_table(physPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physPulmInfInclPresent_df)

# Subset taxa present either in Extubation or Extubation AND Inclusion
physPulmInfExtub <- subset_samples(physPulmInf_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                     "Extubation")

# Keep only represented taxa
physPulmInfExtubPresent <- prune_taxa(taxa_sums(physPulmInfExtub) > 0, 
                                      physPulmInfExtub)

physPulmInfExtubPresent_df <- data.frame(tax_table(physPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physPulmInfInclPresent_df, 
                                 physPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physPulmInfExtubPresent_df, 
                             physPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physPulmInfInclPresent_df$ASV, physPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                fill = c("gold2", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInclExtubAbs_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Inclusion + Extub) available ----

```{r}
# Compare Inclusion with Extubation
PatExtraPulmInf_InclExtub <- c("3", "16", "18", "32")

physExtraPulmInf_InclExtub <- subset_samples(physeqGutHum, record_ID %in% PatExtraPulmInf_InclExtub)

# Keep only represented taxa
physExtraPulmInf_InclExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInf_InclExtub) > 0, 
                                                physExtraPulmInf_InclExtub)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InclExtubPresent)

# Subset taxa present either in Inclusion or Inclusion AND Extubation
physExtraPulmInfIncl <- subset_samples(physExtraPulmInf_InclExtubPresent,
                                       Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physExtraPulmInfInclPresent <- prune_taxa(taxa_sums(physExtraPulmInfIncl) > 0, 
                                          physExtraPulmInfIncl)

physExtraPulmInfInclPresent_df <- data.frame(tax_table(physExtraPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInclPresent_df)

# Subset taxa present either in Extubation or Extubation AND Inclusion
physExtraPulmInfExtub <- subset_samples(physExtraPulmInf_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                          "Extubation")

# Keep only represented taxa
physExtraPulmInfExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInfExtub) > 0, 
                                           physExtraPulmInfExtub)

physExtraPulmInfExtubPresent_df <- data.frame(tax_table(physExtraPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physExtraPulmInfInclPresent_df, 
                                 physExtraPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclusionOnly_df)
taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physExtraPulmInfExtubPresent_df, 
                             physExtraPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physExtraPulmInfInclPresent_df$ASV, physExtraPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                fill = c("gold2", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInclExtubAbs_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Infection_D1 + Infection_D5) available ----

```{r}
PatPulmInf_InfD1InfD5 <- c("1","2", "4", "8","9","27", "28", "30", "36", "37", "41")

physPulmInf_InfD1InfD5 <- subset_samples(physeqGutHum, record_ID %in% PatPulmInf_InfD1InfD5)

# Keep only represented taxa
physPulmInf_InfD1InfD5Present <- prune_taxa(taxa_sums(physPulmInf_InfD1InfD5) > 0, 
                                            physPulmInf_InfD1InfD5)

TotalTaxaPresent <- ntaxa(physPulmInf_InfD1InfD5Present)

# Subset taxa present either at Infection_D1 or (Infection D1 AND Infection D5)
physPulmInfInfD1 <- subset_samples(physPulmInf_InfD1InfD5Present,
                                   Filtered_Data_AK_Dec2020 %in% "Infection D1")

# Keep only represented taxa
physPulmInfInfD1Present <- prune_taxa(taxa_sums(physPulmInfInfD1) > 0, 
                                      physPulmInfInfD1)

physPulmInfInfD1Present_df <- data.frame(tax_table(physPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInfD1Present_df)[1] # for area 1 in venn diagram
head(physPulmInfInfD1Present_df)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Infection_D1)
physPulmInfInfD5 <- subset_samples(physPulmInf_InfD1InfD5Present, Filtered_Data_AK_Dec2020 %in% 
                                     "Infection D5")

# Keep only represented taxa
physPulmInfInfD5Present <- prune_taxa(taxa_sums(physPulmInfInfD5) > 0, 
                                      physPulmInfInfD5)

physPulmInfInfD5Present_df <- data.frame(tax_table(physPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfInfD5Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physPulmInfInfD1Present_df, 
                             physPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
taxInfD1Only_df[1:10,]

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physPulmInfInfD5Present_df, 
                             physPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
head(taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
taxInfD1_AND_InfD5 <- intersect(physPulmInfInfD1Present_df$ASV, physPulmInfInfD5Present_df$ASV)
Venn_cross_area <- length(taxInfD1_AND_InfD5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Infection_D1 + Infection_D5) available ----

```{r}
PatExtraPulmInf_InfD1InfD5 <- c("3", "6", "13", "16", "18", "23", "25","29", "32", "34")

physExtraPulmInf_InfD1InfD5 <- subset_samples(physeqGutHum, record_ID %in% PatExtraPulmInf_InfD1InfD5)

# Keep only represented taxa
physExtraPulmInf_InfD1InfD5Present <- prune_taxa(taxa_sums(physExtraPulmInf_InfD1InfD5) > 0, 
                                                 physExtraPulmInf_InfD1InfD5)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InfD1InfD5Present)

# Subset taxa present either at Infection_D1 or (Infection D1 AND Infection D5)
physExtraPulmInfInfD1 <- subset_samples(physExtraPulmInf_InfD1InfD5Present,
                                        Filtered_Data_AK_Dec2020 %in% "Infection D1")

# Keep only represented taxa
physExtraPulmInfInfD1Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD1) > 0, 
                                           physExtraPulmInfInfD1)

physExtraPulmInfInfD1Present_df <- data.frame(tax_table(physExtraPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInfD1Present_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInfD1Present_df)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Infection_D1)
physExtraPulmInfInfD5 <- subset_samples(physExtraPulmInf_InfD1InfD5Present, Filtered_Data_AK_Dec2020 %in% 
                                          "Infection D5")

# Keep only represented taxa
physExtraPulmInfInfD5Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD5) > 0, 
                                           physExtraPulmInfInfD5)

physExtraPulmInfInfD5Present_df <- data.frame(tax_table(physExtraPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfInfD5Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physExtraPulmInfInfD1Present_df, 
                             physExtraPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
taxInfD1Only_df[1:10,]

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physExtraPulmInfInfD5Present_df, 
                             physExtraPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
head(taxInfD5Only_df)

# Taxa present at Infection_D1 and Infection_D5
taxInfD1_AND_InfD5 <- intersect(physExtraPulmInfInfD1Present_df$ASV, physExtraPulmInfInfD5Present_df$ASV)
Venn_cross_area <- length(taxInfD1_AND_InfD5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Infection_D1"="coral2", "Infection_D5"="chartreuse4"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection D1", "Infection_D5"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Inclusion + Infection_D1) available ----

```{r}
PatPulmInf_InclInfD1 <- c("1", "9","26", "28", "33", "40", "41")

physPulmInf_InclInfD1 <- subset_samples(physeqGutHum, record_ID %in% PatPulmInf_InclInfD1)

# Keep only represented taxa
physPulmInf_InclInfD1Present <- prune_taxa(taxa_sums(physPulmInf_InclInfD1) > 0, 
                                           physPulmInf_InclInfD1)

TotalTaxaPresent <- ntaxa(physPulmInf_InclInfD1Present)

# Subset taxa present either at Inclusion or (Inclusion AND Infection_D1)
physPulmInfIncl <- subset_samples(physPulmInf_InclInfD1Present,
                                  Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physPulmInfInclPresent <- prune_taxa(taxa_sums(physPulmInfIncl) > 0, 
                                     physPulmInfIncl)

physPulmInfInclPresent_df <- data.frame(tax_table(physPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physPulmInfInclPresent_df)

# Subset taxa present either at Infection_D1 or (Infection_D1 AND Inclusion)
physPulmInfInfD1 <- subset_samples(physPulmInf_InclInfD1Present, Filtered_Data_AK_Dec2020 %in% 
                                     "Infection D1")

# Keep only represented taxa
physPulmInfInfD1Present <- prune_taxa(taxa_sums(physPulmInfInfD1) > 0, 
                                      physPulmInfInfD1)

physPulmInfInfD1Present_df <- data.frame(tax_table(physPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfInfD1Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclOnly_df <- anti_join(physPulmInfInclPresent_df, 
                            physPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclOnly_df)
taxInclOnly_df[1:10,]

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physPulmInfInfD1Present_df, 
                             physPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
head(taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
taxIncl_AND_InfD1 <- intersect(physPulmInfInclPresent_df$ASV, physPulmInfInfD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_InfD1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Inclusion" = "gold2", "Infection_D1" = "coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection_D1"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Inclusion + Infection_D1) available ----

```{r}
PatExtraPulmInf_InclInfD1 <- c("3", "16", "18", "32")

physExtraPulmInf_InclInfD1 <- subset_samples(physeqGutHum, record_ID %in% PatExtraPulmInf_InclInfD1)

# Keep only represented taxa
physExtraPulmInf_InclInfD1Present <- prune_taxa(taxa_sums(physExtraPulmInf_InclInfD1) > 0, 
                                                physExtraPulmInf_InclInfD1)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InclInfD1Present)

# Subset taxa present either at Inclusion or (Inclusion AND Infection_D1)
physExtraPulmInfIncl <- subset_samples(physExtraPulmInf_InclInfD1Present,
                                       Filtered_Data_AK_Dec2020 %in% "Inclusion")

# Keep only represented taxa
physExtraPulmInfInclPresent <- prune_taxa(taxa_sums(physExtraPulmInfIncl) > 0, 
                                          physExtraPulmInfIncl)

physExtraPulmInfInclPresent_df <- data.frame(tax_table(physExtraPulmInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInclPresent_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInclPresent_df)

# Subset taxa present either at Infection_D1 or (Infection_D1 AND Inclusion)
physExtraPulmInfInfD1 <- subset_samples(physExtraPulmInf_InclInfD1Present, Filtered_Data_AK_Dec2020 %in% 
                                          "Infection D1")

# Keep only represented taxa
physExtraPulmInfInfD1Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD1) > 0, 
                                           physExtraPulmInfInfD1)

physExtraPulmInfInfD1Present_df <- data.frame(tax_table(physExtraPulmInfInfD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfInfD1Present_df)[1] # for area 2 in venn diagram

# Get taxa names present at Inclusion ONLY
taxInclOnly_df <- anti_join(physExtraPulmInfInclPresent_df, 
                            physExtraPulmInfInfD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInclOnly_df)
taxInclOnly_df[1:10,]

# Get taxa names present at Infection_D1 ONLY
taxInfD1Only_df <- anti_join(physExtraPulmInfInfD1Present_df, 
                             physExtraPulmInfInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD1Only_df)
head(taxInfD1Only_df)

# Taxa present at Inclusion and Infection_D1
taxIncl_AND_InfD1 <- intersect(physExtraPulmInfInclPresent_df$ASV, physExtraPulmInfInfD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_InfD1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Inclusion" = "gold2", "Infection_D1" = "coral2"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Inclusion", "Infection_D1"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD1InfD5_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Pulmonary infection AND (Infection_D5 + Extubation) available ----

```{r}
PatPulmInf_InfD5Extub <- c("2", "4", "8", "28", "30", "37")

physPulmInf_InfD5Extub <- subset_samples(physeqGutHum, record_ID %in% PatPulmInf_InfD5Extub)

# Keep only represented taxa
physPulmInf_InfD5ExtubPresent <- prune_taxa(taxa_sums(physPulmInf_InfD5Extub) > 0, 
                                            physPulmInf_InfD5Extub)

TotalTaxaPresent <- ntaxa(physPulmInf_InfD5ExtubPresent)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Extubation)
physPulmInfInfD5 <- subset_samples(physPulmInf_InfD5ExtubPresent,
                                   Filtered_Data_AK_Dec2020 %in% "Infection D5")

# Keep only represented taxa
physPulmInfInfD5Present <- prune_taxa(taxa_sums(physPulmInfInfD5) > 0, 
                                      physPulmInfInfD5)

physPulmInfInfD5Present_df <- data.frame(tax_table(physPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physPulmInfInfD5Present_df)[1] # for area 1 in venn diagram
head(physPulmInfInfD5Present_df)

# Subset taxa present either at Extubation or (Extubation AND Infection_D5)
physPulmInfExtub <- subset_samples(physPulmInf_InfD5ExtubPresent, Filtered_Data_AK_Dec2020 %in% 
                                     "Extubation")

# Keep only represented taxa
physPulmInfExtubPresent <- prune_taxa(taxa_sums(physPulmInfExtub) > 0, 
                                      physPulmInfExtub)

physPulmInfExtubPresent_df <- data.frame(tax_table(physPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physPulmInfInfD5Present_df, 
                             physPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
taxInfD5Only_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physPulmInfExtubPresent_df, 
                             physPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
taxInfD5_AND_Extub <- intersect(physPulmInfInfD5Present_df$ASV, physPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxInfD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("chartreuse4", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection_D5", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD5Extub_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Venn diagram - Patients with Extra-pulmonary infection AND (Infection_D5 + Extubation) available ----

```{r}
PatExtraPulmInf_InfD5Extub <- c("3", "13", "16", "18", "25", "32") 

physExtraPulmInf_InfD5Extub <- subset_samples(physeqGutHum, record_ID %in% PatExtraPulmInf_InfD5Extub)

# Keep only represented taxa
physExtraPulmInf_InfD5ExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInf_InfD5Extub) > 0, 
                                                 physExtraPulmInf_InfD5Extub)

TotalTaxaPresent <- ntaxa(physExtraPulmInf_InfD5ExtubPresent)

# Subset taxa present either at Infection_D5 or (Infection_D5 AND Extubation)
physExtraPulmInfInfD5 <- subset_samples(physExtraPulmInf_InfD5ExtubPresent,
                                        Filtered_Data_AK_Dec2020 %in% "Infection D5")

# Keep only represented taxa
physExtraPulmInfInfD5Present <- prune_taxa(taxa_sums(physExtraPulmInfInfD5) > 0, 
                                           physExtraPulmInfInfD5)

physExtraPulmInfInfD5Present_df <- data.frame(tax_table(physExtraPulmInfInfD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area1 <- dim(physExtraPulmInfInfD5Present_df)[1] # for area 1 in venn diagram
head(physExtraPulmInfInfD5Present_df)

# Subset taxa present either at Extubation or (Extubation AND Infection_D5)
physExtraPulmInfExtub <- subset_samples(physExtraPulmInf_InfD5ExtubPresent,
                                        Filtered_Data_AK_Dec2020 %in% "Extubation")

# Keep only represented taxa
physExtraPulmInfExtubPresent <- prune_taxa(taxa_sums(physExtraPulmInfExtub) > 0, 
                                           physExtraPulmInfExtub)

physExtraPulmInfExtubPresent_df <- data.frame(tax_table(physExtraPulmInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physExtraPulmInfExtubPresent_df)[1] # for area 2 in venn diagram

# Get taxa names present at Infection_D5 ONLY
taxInfD5Only_df <- anti_join(physExtraPulmInfInfD5Present_df, 
                             physExtraPulmInfExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxInfD5Only_df)
taxInfD5Only_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physExtraPulmInfExtubPresent_df, 
                             physExtraPulmInfInfD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Infection_D5 and Extubation
taxInfD5_AND_Extub <- intersect(physExtraPulmInfInfD5Present_df$ASV, physExtraPulmInfExtubPresent_df$ASV)
Venn_cross_area <- length(taxInfD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Infection_D5" = "chartreuse4", "Extubation" = "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Infection_D5", "Extubation"),
                                # Numbers
                                cex = 2, print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammInfD5Extub_Lung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

### Heatmap Patient3 with represented taxa only ----

```{r}
Pat3Present <- prune_taxa(taxa_sums(Pat3) > 0, Pat3)
Pat3Present_PHeat1_20210125 <- plot_heatmap(Pat3Present, "NMDS", "bray", sample.label="UnXlabel", 
                                            taxa.label="UnYlabel", max.label = 240)

ggsave(filename = "Pat3Present_PHeat1_20210125_high_tax_names.pdf",
       width = 15, height = 45, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

ggsave(filename = "Pat3_PHeat1_20210125.pdf",
       width = 12, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

Pat3Present_PHeat1_20210125[["plot_env"]][["sample.order"]]
Pat3Present_PHeat1_20210125[["plot_env"]][["taxa.order"]]
tail(Pat3Present_PHeat1_20210125[["plot_env"]][["taxa.order"]])
```

### Heatmap Patient28 with represented taxa only ----

```{r}
Pat28Present <- prune_taxa(taxa_sums(Pat28) > 0, Pat28)
Pat28Present_PHeat1_20210127 <- plot_heatmap(Pat28Present, "NMDS", "bray", sample.label="UnXlabel", 
                                             taxa.label="UnYlabel", max.label = 166)

ggsave(filename = "Pat28Present_PHeat1_20210127_high_tax_names.pdf",
       width = 15, height = 45, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

ggsave(filename = "Pat28_PHeat1_20210127.pdf",
       width = 15, height = 15, dpi = 200, units = "cm", device='pdf', limitsize = FALSE)

Pat28Present_PHeat1_20210127[["plot_env"]][["sample.order"]]
Pat28Present_PHeat1_20210127[["plot_env"]][["taxa.order"]]
tail(Pat28Present_PHeat1_20210127[["plot_env"]][["taxa.order"]])
```

### Alpha diversity - Richness ----

/#/`/`/`{r} psGutAlphaRich_df /<- psmelt(physeqGutHumRel)

names(psGutAlphaRich_df)

psGutAlphaRich_df/<-subset(psGutAlphaRich_df,!(record_ID=="38"))

psGutAlphaRich_df %/</>% select(OTU,record_ID, Sample, Abundance, Phylum, Filtered_Data_AK_Dec2020, Shannon, InvSimpson,Chao1, Observed,POCGroups20200115) %/>% convert(fct(OTU, Sample, Phylum, Filtered_Data_AK_Dec2020)) %/>% arrange(-(Abundance))

psGutAlphaRich_df %/>% ggplot(aes(x = POCGroups20200115, y = Observed)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Group of patients", y = "Observed") + theme_classic()

Shannon_gut/<-psGutAlphaRich_df %/>% ggplot(aes(x = POCGroups20200115, y = Shannon)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Patient Group", y = "Shannon") + theme_classic()

ggsave("Shannon_Gut_Groups.tiff", width = 14, height = 9, dpi = 300, units = "cm")

Chao1_gut/<-psGutAlphaRich_df %/>% ggplot(aes(x = POCGroups20200115, y = Chao1)) + geom_boxplot( aes(fill = Filtered_Data_AK_Dec2020), alpha = 0.4, position = position_dodge(0.9), outlier.colour = NULL, outlier.size = 1.5, outlier.stroke = 0, outlier.alpha = 0.5, ) + scale_fill_manual(values = InclInfectColors5) + labs(x = "Patient Group", y = "Chao1") + theme_classic()

ggsave("Chao1_Gut_Groups.tiff", width = 14, height = 9, dpi = 300, units = "cm")

ggarrange(Shannon_gut,Chao1_gut, labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )

ggsave("Alpha_diversity_combined_gut.tiff", width = 18, height = 9, dpi = 300, units = "cm")

### Distance - beta diversity ----
#```{r}
# Adapted from https://joey711.github.io/phyloseq/distance.html
# Visualize sample distribution according to various distance metrics
dist_methods <- unlist(phyloseq::distanceMethodList)
print(dist_methods)

dist_methods <- dist_methods[-which(dist_methods=="ANY")]

plist <- vector("list", length(dist_methods))
names(plist) <-  dist_methods
names(plist)

for( i in dist_methods ){
  # Calculate distance matrix (distance function requires phyloseq::)
  iDist <- phyloseq::distance(physeqGutHum, method=i)
  # Calculate ordination
  iMDS  <- ordinate(physeqGutHum, "NMDS", distance=bin)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(physeqGutHum, iMDS, color="Filtered_Data_AK_Dec2020")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

multidiv <-ggplot(df, aes(Axis.1, Axis.2, color=Filtered_Data_AK_Dec2020)) + 
  geom_point(size=3, alpha=0.5) + 
  facet_wrap(~distance, scales="free") +
  ggtitle("MDS on various distance metrics")

multidiv

ggsave(filename = "Beta_diversity_multidiv.pdf",
       width = 45, height = 45, dpi = 200, units = "cm", device='pdf')

print(plist[["binary"]])
print(plist[["unifrac"]])

distBin <- phyloseq::distance(physeqGutHum, method="binary")

physeqGutHum_ord <- ordinate(physeqGutHum, "NMDS", distBin)

plot_ord_MDS_bin <- plot_ordination(physeqGutHum, physeqGutHum_ord, type="samples",
                                    color="Filtered_Data_AK_Dec2020") +
  scale_color_manual(values = InclInfectColors5) +
  theme_classic()+
  stat_ellipse(type="t")


ggsave(filename = "plot_ord_MDS_bin_InclInf_Gut.pdf",
       width = 12, height = 6, dpi = 200, units = "cm", device='pdf')
# stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(POCGroups20200115))  # this is just another form of elipse


# Assign a different shape to No infection, Pneumonia and Extra-pulmonary infection 
plot_ord_MDS_bin2 <- plot_ordination(physeqGutHum, physeqGutHum_ord, type="samples",
                                     axes=c(1,2), color="Filtered_Data_AK_Dec2020", shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Binary",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_bin2

ggsave(filename = "plot_ord_MDS_bin2_InclInf_shape_Gut.pdf",
       width = 12, height = 6, dpi = 200, units = "cm", device='pdf')



#weighted unifrac distance

# First subset physeqGutHum object to exlcude sample 38 which was not included and has an NA for the POCGroups20200115 variable

physeqGutHum2<-subset_samples(physeqGutHum,!(is.na(POCGroups20200115)))


View(psmelt(physeqGutHum2))


#weighted unifrac distance
physeqGutHum_ord2 <- ordinate(physeqGutHum2, "NMDS", "wunifrac")

plot_ord_MDS_wunifrac <- plot_ordination(physeqGutHum2, physeqGutHum_ord2, type="samples", axes=c(1,2), color="Filtered_Data_AK_Dec2020", label="SampleID",shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Weighted Unifrac",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_wunifrac

# Repeat the ordination plot after excludin the outlier S66 as shown in the graph

physeqGutHum3<-subset_samples(physeqGutHum,!(SampleID=="S1804"| SampleID=="S1804_rep"|SampleID=="S901"|SampleID=="S66")) # exclude sample S66 which is an outlier

physeqGutHum3<-subset_samples(physeqGutHum3,!(is.na(POCGroups20200115)))

View(psmelt(physeqGutHum))


###  Check for beta-diversity of the 3 groups only at inclusion

physeqGutHum3incl<-subset_samples(physeqGutHum3,Filtered_Data_AK_Dec2020=="Inclusion")

View(psmelt(physeqGutHum3incl))

# check for Bray distance
physeqGutHum_ordincl <- ordinate(physeqGutHum3incl, "NMDS", "bray")

bray_gut_incl <- plot_ordination(physeqGutHum3incl, physeqGutHum_ordincl, type="samples", axes=c(1,2), color="POCGroups20200115") +
  scale_color_brewer(palette="Dark2") +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Gut samples",subtitle="Bray-Curtis distance only at Inclusion",colour="Group of patients")

bray_gut_incl_elps<-bray_gut_incl+stat_ellipse(type="t")  # ad ellipses
bray_gut_incl_elps

# Repeat the same graph for jaccard

physeqGutHum_ordincljaccard <- ordinate(physeqGutHum3incl, "NMDS", "jaccard")

jaccard_gut_incl <- plot_ordination(physeqGutHum3incl,physeqGutHum_ordincljaccard, type="samples", axes=c(1,2), color="POCGroups20200115") +
  scale_color_brewer(palette="Dark2") +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Gut samples",subtitle="Jaccard distance only at Inclusion",colour="Group of patients")

jaccard_gut_incl

jaccard_gut_incl_elps<-jaccard_gut_incl+stat_ellipse(type="t")  # ad ellipses

jaccard_gut_incl_elps



#Combine graphs
ggarrange(bray_gut_incl_elps,jaccard_gut_incl_elps,
          labels = c("A", "B"),common.legend = TRUE, legend = "bottom" )

ggsave("Beta_diversity_Gut_inclusion.tiff",
       width = 24, height = 9, dpi = 300, units = "cm")



#wunifrac distance
physeqGutHum_ord2 <- ordinate(physeqGutHum3, "NMDS", "wunifrac")

plot_ord_MDS_wunifrac <- plot_ordination(physeqGutHum3, physeqGutHum_ord2, type="samples", axes=c(1,2),color="Filtered_Data_AK_Dec2020",shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Weighted Unifrac",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_wunifrac

wunifracdist<-plot_ord_MDS_wunifrac+stat_ellipse(type="t")  # ad ellipses

wunifracdist<-wunifracdist+facet_wrap(~POCGroups20200115) # combined graphs
wunifracdist


#Bray distance

physeqGutHum_ord3 <- ordinate(physeqGutHum3, "NMDS", "bray")

plot_ord_MDS_bray <- plot_ordination(physeqGutHum3, physeqGutHum_ord3, type="samples", axes=c(1,2), color="Filtered_Data_AK_Dec2020",shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Bray-Curtis",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_bray

braydist<-plot_ord_MDS_bray+stat_ellipse(type="t")  # ad ellipses

braydist<-braydist+facet_wrap(~POCGroups20200115) # combined graphs
braydist



# Taxa ordination according to patient group based on Bray-Curtis distance
taxord<-plot_ordination(physeqGutHum3, physeqGutHum_ord3, type="split", 
                        color="Phylum",title="Taxa ordination (distance:Bray-Curtis)")

taxord+facet_wrap(~Phylum)



gg_ordiplot(rda(otu_table(physeqGutHum3)),groups=sample_data(physeqGutHum3)$POCGroups20200115, ellipse=F, spider=T)

taxord2<-plot_ordination(fil, physeqGutHum_ord3,color="POCGroups20200115")

fil<-tax_glom(physeqGutHum3,"Family")

taxord2

taxord2+facet_wrap(~Phylum)

ordiplot(taxord2,dis="si",type="n")

points (taxord2, shape = physeqGutHum3$POCGroups20200115, pch = physeqGutHum3$POCGroups20200115 )

centroid.pcoa

sample_data(physeqGutHum3)

View(psmelt(physeqGutHum3))

sample_data(physeqGutHum3)
rank_names(physeqGutHum3)


taxord+facet_wrap(~Phylum)








#Jaccard distance


physeqGutHum_ord4 <- ordinate(physeqGutHum3, "NMDS", "jaccard")

plot_ord_MDS_jaccard <- plot_ordination(physeqGutHum3, physeqGutHum_ord4, type="samples", axes=c(1,2), color="Filtered_Data_AK_Dec2020",shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Non-metric multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Jaccard",colour="Sample timepoint",shape="Study group")

plot_ord_MDS_jaccard

jaccarddist<-plot_ord_MDS_jaccard+stat_ellipse(type="t")  # ad ellipses

jaccarddist<-jaccarddist+facet_wrap(~POCGroups20200115) # combined graphs

jaccarddist




###################


plot_ord_MDS_wunifrac_onlytimepoints <- plot_ordination(physeqGutHum2, physeqGutHum_ord2, type="samples",
                                                        axes=c(1,2), color="Filtered_Data_AK_Dec2020") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Weighted Unifrac",colour="Sample timepoint")+
  stat_ellipse(type="t")

plot_ord_MDS_wunifrac_onlytimepoints 


plot_ord_MDS_wunifrac_groups <- plot_ordination(physeqGutHum2, physeqGutHum_ord2, type="samples",
                                                axes=c(1,2), color="POCGroups20200115") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw()+
  labs(title="Multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Weighted Unifrac",colour="Group of patients")+
  stat_ellipse(type="t")

plot_ord_MDS_wunifrac_groups 

ggarrange(plot_ord_MDS_wunifrac,ggarrange(plot_ord_MDS_wunifrac_onlytimepoints,plot_ord_MDS_wunifrac_groups,ncol=2),nrow=2)


#Jaccard distance
physeqGutHum_ord3 <- ordinate(physeqGutHum2, "MDS", "jaccard")
plot_ord_MDS_jaccard <- plot_ordination(physeqGutHum2, physeqGutHum_ord3, type="samples",
                                        axes=c(1,2), color="Filtered_Data_AK_Dec2020", shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Jaccard",colour="Sample timepoint",shape="Study group")


jaccard1<-plot_ord_MDS_jaccard+stat_ellipse(type="t")

jaccard1<-jaccard1+facet_wrap(~POCGroups20200115)


#Bray-Curtis distance
physeqGutHum_ord4 <- ordinate(physeqGutHum2, "MDS", "bray")
plot_ord_MDS_bray <- plot_ordination(physeqGutHum2, physeqGutHum_ord4, type="samples",
                                     axes=c(1,2), color="Filtered_Data_AK_Dec2020", shape="POCGroups20200115") +
  scale_color_manual(values = InclInfectColors5) +
  theme_bw()+
  labs(title="Multidimensional Scale Plot for Gut samples",subtitle="Distance metric=Bray-Curtis",colour="Sample timepoint",shape="Study group")

bray1<-plot_ord_MDS_bray+stat_ellipse(type="t")

bray1<-bray1+facet_wrap(~POCGroups20200115)

# Plot in a single graph Wunifrac, Jaccard and Bray distances
ggarrange(wunifracdist,jaccarddist,braydist,labels=c("A","B","C"),common.legend=TRUE,nrow=3,legend="right")


ggsave(filename = "Beta_diversity_gut_combined.pdf",
       width = 24, height = 20, dpi = 300, units = "cm", device='pdf')



sample_data(physeqGutHum2)$record_ID %>% unique() %>% length()  #number of patients included in the analysis (37 in this case)
sample_data(physeqGutHum2)$Sample %>% unique() %>% length()   # number of samples analysed (106 in this case)

### Differential abundance - DESeq2 - Total patients - Inclusion vs Infection_D5 ----

```{r}
# From https://joey711.github.io/phyloseq-extensions/DESeq2.html
# Subset samples of Inclusion and InfectD5
InclInfD5_ps <- subset_samples(physeqGutHumUn,
                               Filtered_Data_AK_Dec2020 %in% c("Inclusion", "Infection D5"))

InclInfD5_present <- prune_taxa(taxa_sums(InclInfD5_ps)>0, InclInfD5_ps)

diagdds <-  phyloseq_to_deseq2(InclInfD5_present, ~ Filtered_Data_AK_Dec2020)

# calculate geometric means prior to estimate size factors
gm_mean  <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(diagdds), 1, gm_mean)
diagdds  <-  estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds  <-  DESeq(diagdds, fitType="local")

res <-  results(diagdds, cooksCutoff = FALSE)
alpha  <-  0.05
sigtab  <-  res[which(res$padj < alpha), ]
sigtab  <-  cbind(as(sigtab, "data.frame"), as(tax_table(InclInfD5_present)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_classic())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, colour=Phylum)) +
  geom_point(size=4, alpha=0.7) +
  scale_colour_manual(values = PhylumPalGut13) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(filename = "DESeq2_AllpatientsInclInfD5_Gut.pdf",
       width = 15, height = 12, dpi = 200, units = "cm", device='pdf')
```

### Differential abundance - DESeq2 - Patients with Pulmonary infection AND (Infection_D1, Infection_D5 AND/OR Extubation) available ----

```{r}
# From https://joey711.github.io/phyloseq-extensions/DESeq2.html
# Subset samples of Inclusion and InfectD5
PatPulmInf_InfD1InfD5 <- c("1","2","4", "8","9","27", "28", "30", "36", "37", "41")

physPulmInf_InfD1InfD5 <- subset_samples(physeqGutHum, record_ID %in% PatPulmInf_InfD1InfD5)

# Keep only represented taxa
physPulmInf_InfD1InfD5Present <- prune_taxa(taxa_sums(physPulmInf_InfD1InfD5) > 0, 
                                            physPulmInf_InfD1InfD5)

diagdds <-  phyloseq_to_deseq2(physPulmInf_InfD1InfD5Present, ~ Filtered_Data_AK_Dec2020)

# calculate geometric means prior to estimate size factors
gm_mean  <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(diagdds), 1, gm_mean)
diagdds  <-  estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds  <-  DESeq(diagdds, fitType="local")

res <-  results(diagdds, cooksCutoff = FALSE)
alpha  <-  0.05
sigtab  <-  res[which(res$padj < alpha), ]
sigtab  <-  cbind(as(sigtab, "data.frame"), as(tax_table(physPulmInf_InfD1InfD5Present)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_classic())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, colour=Phylum)) +
  geom_point(size=2.5, alpha=0.7) +
  scale_colour_manual(values = PhylumPalGut13) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(filename = "DESeq2_PulmInf_InfD1InfD5_Gut.pdf",
       width = 15, height = 12, dpi = 200, units = "cm", device='pdf')
```

### Differential abundance - DESeq2 - Patients with Extrapulmonary infection AND (Infection_D1, Infection_D5 AND/OR Extubation) available ----

```{r}
# From https://joey711.github.io/phyloseq-extensions/DESeq2.html
# Subset samples of Inclusion and InfectD5
PatExtraPulmInf_InfD1InfD5 <- c("3","6","13","16","18","23","25","29","32","34")

physExtraPulmInf_InfD1InfD5 <- subset_samples(physeqGutHum, record_ID %in% PatExtraPulmInf_InfD1InfD5)

# Keep only represented taxa
physExtraPulmInf_InfD1InfD5Present <- prune_taxa(taxa_sums(physExtraPulmInf_InfD1InfD5) > 0, 
                                                 physExtraPulmInf_InfD1InfD5)

diagdds <-  phyloseq_to_deseq2(physExtraPulmInf_InfD1InfD5Present, ~ Filtered_Data_AK_Dec2020)

# calculate geometric means prior to estimate size factors
gm_mean  <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(diagdds), 1, gm_mean)
diagdds  <-  estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds  <-  DESeq(diagdds, fitType="local")

res <-  results(diagdds, cooksCutoff = FALSE)
alpha  <-  0.05
sigtab  <-  res[which(res$padj < alpha), ]
sigtab  <-  cbind(as(sigtab, "data.frame"), as(tax_table(physExtraPulmInf_InfD1InfD5Present)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_classic())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Family, y=log2FoldChange, colour=Phylum)) +
  geom_point(size=2.5, alpha=0.7) +
  scale_colour_manual(values = PhylumPalGut13) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(filename = "DESeq2_ExtraPulmInf_InfD1InfD5_Gut_Fam.pdf",
       width = 15, height = 12, dpi = 200, units = "cm", device='pdf')

```

### Differential abundance - DESeq2 - Patients with No infection AND (Inclusion + Extub) available ----

```{r}
# From https://joey711.github.io/phyloseq-extensions/DESeq2.html
# Subset samples
PatNoInf_InclExtub <- c("7","20","24","39")

physNoInf_InclExtub <- subset_samples(physeqGutHum, record_ID %in% PatNoInf_InclExtub)

# Keep only represented taxa
physNoInf_InclExtubPresent <- prune_taxa(taxa_sums(physNoInf_InclExtub) > 0, 
                                         physNoInf_InclExtub)

sample_data(physNoInf_InclExtubPresent)

diagdds <-  phyloseq_to_deseq2(physNoInf_InclExtubPresent, ~ Filtered_Data_AK_Dec2020)

# calculate geometric means prior to estimate size factors
gm_mean  <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <-  apply(counts(diagdds), 1, gm_mean)
diagdds  <-  estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds  <-  DESeq(diagdds, fitType="local")

res <-  results(diagdds, cooksCutoff = FALSE)
alpha  <-  0.05 # significance level
sigtab  <-  res[which(res$padj < alpha), ]
sigtab  <-  cbind(as(sigtab, "data.frame"), as(tax_table(physNoInf_InclExtubPresent)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_classic())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, colour=Phylum)) +
  geom_point(size=3, alpha=0.6) +
  scale_colour_manual(values = PhylumPal) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(filename = "DESeq2_NoInf_InclExtub_Lung.pdf",
       width = 15, height = 12, dpi = 200, units = "cm", device='pdf')
```

### Taxa unique to Patients with No infection AND (Inclusion + Extub) available ----

```{r}
# See Venn diagram - Patients with No infection AND (Inclusion + Extub) available
# Inclusion only
InclOnlyASVNames <- taxInclusionOnly_df$ASV

physInclOnly <- prune_taxa(InclOnlyASVNames, physNoInfIncl)
tax_table(physInclOnly)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physInclOnlyRel <- transform_sample_counts(physInclOnly, count_to_rel_abund)

# Agglomeration
physInclOnlyAgg <- tax_glom(physInclOnlyRel, "Genus", NArm = TRUE)

# Top 20 genera of taxa exclusively represented at Inclusion
TopTaxaNames <-  names(sort(taxa_sums(physInclOnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physInclOnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPalLung13) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusInclOnly_Lung_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Extubation only
ExtubOnlyASVNames <- taxExtubOnly_df$ASV

physExtubOnly <- prune_taxa(ExtubOnlyASVNames, physNoInfExtub)
tax_table(physExtubOnly)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physExtubOnlyRel <- transform_sample_counts(physExtubOnly, count_to_rel_abund)

# Agglomeration
physExtubOnlyAgg <- tax_glom(physExtubOnlyRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa exclusively represented at Extubation
TopTaxaNames <-  names(sort(taxa_sums(physExtubOnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physExtubOnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPalLung13) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusExtubOnly_Lung_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Taxa shared between Inclusion and Extubation
physShared <- prune_taxa(taxIncl_AND_Extub, physNoInf_InclExtubPresent)
tax_table(physShared)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physSharedRel <- transform_sample_counts(physShared, count_to_rel_abund)

# Agglomeration
physSharedAgg <- tax_glom(physSharedRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa shared between Inclusion and Extubation
TopTaxaNames <-  names(sort(taxa_sums(physSharedAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physSharedAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPalLung13) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusSharedInclExtub_Lung_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")
```

### Taxa unique to Patients with Pulm infection AND (Infection_D1 + Infection_D5) available ----

```{r}
# See Venn diagram - Pulm infection AND (Infection_D1 + Infection_D5) available
# Infection_D1 only
InfD1OnlyASVNames <- taxInfD1Only_df$ASV

physInfD1Only <- prune_taxa(InfD1OnlyASVNames, physPulmInfInfD1Present)
tax_table(physInfD1Only)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physInfD1OnlyRel <- transform_sample_counts(physInfD1Only, count_to_rel_abund)

# Agglomeration
physInfD1OnlyAgg <- tax_glom(physInfD1OnlyRel, "Genus", NArm = TRUE)

# Top 20 genera of taxa exclusively represented at Infection_D1
TopTaxaNames <-  names(sort(taxa_sums(physInfD1OnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physInfD1OnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPal) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusPulmInf_InfD1Only_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Infection_D5 only
InfectD5OnlyASVNames <- taxInfD5Only_df$ASV

physInfD5Only <- prune_taxa(InfectD5OnlyASVNames, physPulmInfInfD5Present)
tax_table(physInfD5Only)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physInfD5OnlyRel <- transform_sample_counts(physInfD5Only, count_to_rel_abund)

# Agglomeration
physInfD5OnlyAgg <- tax_glom(physInfD5OnlyRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa exclusively represented at Infection_D5
TopTaxaNames <-  names(sort(taxa_sums(physInfD5OnlyAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physInfD5OnlyAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPal) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusPulmInf_InfD5Only_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Taxa shared between Infection_D1 and Infection_D5
physShared <- prune_taxa(taxInfD1_AND_InfD5, physPulmInf_InfD1InfD5Present)
tax_table(physShared)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

physSharedRel <- transform_sample_counts(physShared, count_to_rel_abund)

# Agglomeration
physSharedAgg <- tax_glom(physSharedRel, "Genus", NArm = TRUE)

# Identify the Top 20 genera of taxa shared between Infection_D1 and Infection_D5
TopTaxaNames <-  names(sort(taxa_sums(physSharedAgg), decreasing = TRUE)[1:20])
physTopTaxa <-  prune_taxa(TopTaxaNames, physSharedAgg)

TopTaxa_df <- psmelt(physTopTaxa)
str(TopTaxa_df)

TopTaxa_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(TopTaxa_df)

TopTaxaBoxPlot <- TopTaxa_df %>% 
  ggplot(aes(x = reorder(Genus, Abundance, FUN=sum), y=Abundance)) +
  geom_point(aes(colour = Phylum), size = 2.5, alpha=0.5) +
  scale_colour_manual(values = PhylumPal) +
  xlab("") +
  ylim(c(0, 1.0)) +
  labs(x = "Genus" , y = "Relative abundance") +
  coord_flip() +
  theme_classic()

ggsave("GenusPulmInf_SharedInfD1InfD5_Dotplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")
```

############################################################################################ 

############################################################################################ 

######################## THE FOLLOWING SCRIPT IS NOT RUN BECAUSE...

########### DIRECT LUNG AND GUT COMPARISON DOEAS NOT SEEM TO MAKE SENSE

/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#THIS PART COULD BE DISCUSSED LATER /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/# /#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#/#

##### Lung-Gut comparisons - Venn diagram - All patients and samples ----

```{r}
# Subset Lung samples with corresponding Gut samples available
psLung_wGutAvai <- subset_samples(physeqLung9, GutSampleAvailable == "Yes")

# Keep only represented taxa
psLung_wGutAvai_Present <- prune_taxa(taxa_sums(psLung_wGutAvai) > 0, 
                                      psLung_wGutAvai)

# Subset Gut samples with corresponding Lung samples available
psGut_wLungAvai <- subset_samples(physeqGutHum, LungSampleAvailable == "Yes")

# Remove replicates
psGut_wLungAvai2 <-  subset_samples(psGut_wLungAvai, Replicate  == "No")

# Keep only represented taxa
psGut_wLungAvai2_Present <- prune_taxa(taxa_sums(psGut_wLungAvai2) > 0, 
                                       psGut_wLungAvai2)

# Species present in Gut and (Gut and Lung) samples
# in patients with (Gut and Lung) samples available
Gut_GutLung_Species <- get_taxa_unique(psGut_wLungAvai2_Present, "Species")
Venn_area1 <- length(Gut_GutLung_Species)

# Species present in Lung and (Lung and Gut) samples
# in patients with Gut and Lung samples available
Lung_LungGut_Species <- get_taxa_unique(psLung_wGutAvai_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples
# in patients with Gut and Lung samples available
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for all samples with shared taxa
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("raw", "percent"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_diagrammGutLung.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with No infection and Samples at Inclusion ----

```{r}
# Subset samples
psGutIncl <- subset_samples(psGut_wLungAvai2,
                            POCGroups20200115 == "Control" & Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psGutIncl)[, c(7,8)]

psLungIncl <- subset_samples(psLung_wGutAvai,
                             POCGroups20200115 == "Control" &  Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psLungIncl)[, c(4,15)]

# Keep only represented taxa
psGutIncl_Present <- prune_taxa(taxa_sums(psGutIncl) > 0, 
                                psGutIncl)
psLungIncl_Present <- prune_taxa(taxa_sums(psLungIncl) > 0, 
                                 psLungIncl)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutIncl_Present)
otu_table(psGutIncl_Present)[ , "S2104"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungIncl_Present)
otu_table(psLungIncl_Present)[ , "S3905"] %>% unique %>% na.exclude %>% length

# Species present at Inclusion in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutIncl_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Inclusion in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungIncl_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_NoInf_Inclusion.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with No infection and Samples at Extubation ----

```{r}
# Subset samples
psGutExtub <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Control" & Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psGutExtub)[, c(7,8)]

psLungExtub <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Control" &  Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psLungExtub)[, c(4,15)]

# Keep only represented taxa
psGutExtub_Present <- prune_taxa(taxa_sums(psGutExtub) > 0, 
                                 psGutExtub)
psLungExtub_Present <- prune_taxa(taxa_sums(psLungExtub) > 0, 
                                  psLungExtub)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutExtub_Present)
otu_table(psGutExtub_Present)[ , "S3906"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungExtub_Present)
otu_table(psLungExtub_Present)[ , "S3907"] %>% unique %>% na.exclude %>% length

# Species present at Extubation in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutExtub_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Extubation in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungExtub_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_NoInf_Extubation.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Inclusion ----

```{r}
# Subset samples
psGutIncl <- subset_samples(psGut_wLungAvai2,
                            POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psGutIncl)[, c(7,8)]

psLungIncl <- subset_samples(psLung_wGutAvai,
                             POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Inclusion")
sample_data(psLungIncl)[, c(4,15)]

# Keep only represented taxa
psGutIncl_Present <- prune_taxa(taxa_sums(psGutIncl) > 0, 
                                psGutIncl)
psLungIncl_Present <- prune_taxa(taxa_sums(psLungIncl) > 0, 
                                 psLungIncl)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutIncl_Present)
otu_table(psGutIncl_Present)[ , "S3304"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungIncl_Present)
otu_table(psLungIncl_Present)[ , "S4204"] %>% unique %>% na.exclude %>% length

# Species present at Inclusion in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutIncl_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Inclusion in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungIncl_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_Inclusion.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Infection D1 ----

```{r}
# Subset samples
psGutInfD1 <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Infection D1")
sample_data(psGutInfD1)[, c(7,8)]

psLungInfD1 <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Infection_D1")
sample_data(psLungInfD1)[, c(4,15)]

# Keep only represented taxa
psGutInfD1_Present <- prune_taxa(taxa_sums(psGutInfD1) > 0, 
                                 psGutInfD1)
psLungInfD1_Present <- prune_taxa(taxa_sums(psLungInfD1) > 0, 
                                  psLungInfD1)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutInfD1_Present)
otu_table(psGutInfD1_Present)[ , "S3604"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungInfD1_Present)
otu_table(psLungInfD1_Present)[ , "S4305"] %>% unique %>% na.exclude %>% length

# Species present at Infection D1 in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutInfD1_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Infection D1 in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungInfD1_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_InfD1.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Infection D5 ----

```{r}
# Subset samples
psGutInfD5 <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Infection D5")
sample_data(psGutInfD5)[, c(7,8)]

psLungInfD5 <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Infection_D5")
sample_data(psLungInfD5)[, c(4,15)]

# Keep only represented taxa
psGutInfD5_Present <- prune_taxa(taxa_sums(psGutInfD5) > 0, 
                                 psGutInfD5)
psLungInfD5_Present <- prune_taxa(taxa_sums(psLungInfD5) > 0, 
                                  psLungInfD5)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutInfD5_Present)
otu_table(psGutInfD5_Present)[ , "S3610"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungInfD5_Present)
otu_table(psLungInfD5_Present)[ , "S4112"] %>% unique %>% na.exclude %>% length

# Species present at Infection D5 in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutInfD5_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Infection D5 in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungInfD5_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_InfD5.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Venn - Patients with Pulm infection and Samples at Extubation ----

```{r}
# Subset samples
psGutExtub <- subset_samples(psGut_wLungAvai2,
                             POCGroups20200115 == "Pneumonia" & Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psGutExtub)[, c(7,8)]

psLungExtub <- subset_samples(psLung_wGutAvai,
                              POCGroups20200115 == "Pneumonia" &  Filtered_Data_AK_Dec2020 == "Extubation")
sample_data(psLungExtub)[, c(4,15)]

# Keep only represented taxa
psGutExtub_Present <- prune_taxa(taxa_sums(psGutExtub) > 0, 
                                 psGutExtub)
psLungExtub_Present <- prune_taxa(taxa_sums(psLungExtub) > 0, 
                                  psLungExtub)

# Obtain the numbers of taxa per Sample type and per Patient
GutSampleNames <- sample_names(psGutExtub_Present)
otu_table(psGutExtub_Present)[ , "S3012"] %>% unique %>% na.exclude %>% length

LungSampleNames <- sample_names(psLungExtub_Present)
otu_table(psLungExtub_Present)[ , "S3714"] %>% unique %>% na.exclude %>% length

# Species present at Extubation in Gut and (Gut and Lung) samples
Incl_Gut_GutLung_Species <- get_taxa_unique(psGutExtub_Present, "Species")
Venn_area1 <- length(Incl_Gut_GutLung_Species)

# Species present at Extubation in Lung and (Lung and Gut) samples
Lung_LungGut_Species <- get_taxa_unique(psLungExtub_Present, "Species")
Venn_area2 <- length(Lung_LungGut_Species)

# Species shared between Gut and Lung samples at Inclusion
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)
Venn_cross_area <- length(SharedSpecies)

# Total number of taxa
TotalUniqueTaxa <- Venn_area1 + Venn_area2 - Venn_cross_area

# For Venn diagram - Move to new plotting page
grid.newpage()

# Create Venn diagram for samples at Inclusion
# set print.mode to ("raw", "percent")
VennDiagram::draw.pairwise.venn(area1 = Venn_area1,area2 = Venn_area2, cross.area = Venn_cross_area,
                                fill = c("Gut"="burlywood", "Lung"="skyblue1"),
                                alpha = 0.25,
                                lty = "blank",
                                category = c("Gut", "Lung"),
                                cex = 1.5, # Size of numbers
                                print.mode = c("percent", "raw"), sigdigs = 2,
                                fontfamily = "sans",
                                filename = "Venn_GutLung_PulmInf_Extub.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches
```

##### Lung-Gut - Analysis limited to taxa shared between Gut and Lung ----

```{r}
# Work with Objects obtained in Lung-Gut comparisons - Venn diagram - ...
psLung_wGutAvai_Present
psGut_wLungAvai2_Present

# Obtain lists of taxa to subset
SharedSpecies <- intersect(Gut_GutLung_Species, Lung_LungGut_Species)

psLungShared <- subset_taxa(psLung_wGutAvai_Present, Species %in% SharedSpecies)
psGutShared <- subset_taxa(psGut_wLungAvai2_Present, Species %in% SharedSpecies)

# Transform to Relative abundance
count_to_rel_abund <- function(x) {return( x / sum(x))}

psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

# Agglomeration
psLungSharedRel_Gen <- tax_glom(psLungSharedRel, "Genus", NArm = TRUE)
psGutSharedRel_Gen <- tax_glom(psGutSharedRel, "Genus", NArm = TRUE)

# Top x genera with shared taxa in the Lung
TopGenNamesLung <-  names(sort(taxa_sums(psLungSharedRel_Gen), decreasing = TRUE)[1:20])
psLungTopGen <-  prune_taxa(TopGenNamesLung, psLungSharedRel_Gen)

# Prepare a df for subsequent plot
psLungTopGen_df <- psmelt(psLungTopGen)
str(psLungTopGen_df)

# Convert as required
psLungTopGen_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(psLungTopGen_df)

# Grouped boxplot
# Reverse Group levels for vertical plot
psLungTopGen_df$Filtered_Data_AK_Dec2020 <- factor(psLungTopGen_df$Filtered_Data_AK_Dec2020,
                                                   levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungTopGen_df %>%
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Genus", y = "Relative abundance") +
  theme_classic()

ggsave("TopGenLung_sharedGutLung_Boxplot.pdf",
       width = 25, height = 30, dpi = 200, units = "cm")

# Top x genera with shared taxa in the Gut
TopGenNamesGut <-  names(sort(taxa_sums(psGutSharedRel_Gen), decreasing = TRUE)[1:20])
psGutTopGen <-  prune_taxa(TopGenNamesGut, psGutSharedRel_Gen)

# Prepare a df for subsequent plot
psGutTopGen_df <- psmelt(psGutTopGen)
str(psGutTopGen_df)

# Convert as required
psGutTopGen_df %<>% 
  select(Phylum, Family, Genus, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus)) %>% arrange(-(Abundance))

head(psGutTopGen_df)

# Grouped boxplot
# Reverse Group levels for vertical plot
psGutTopGen_df$Filtered_Data_AK_Dec2020 <- factor(psGutTopGen_df$Filtered_Data_AK_Dec2020,
                                                  levels = c("Discharge", "Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutTopGen_df %>%
  ggplot(aes(x = reorder(Genus, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Genus", y = "Relative abundance") +
  theme_classic()

ggsave("TopGenGut_sharedGutLung_Boxplot.pdf",
       width = 25, height = 30, dpi = 200, units = "cm")

# Analysis of Prevotella species in Lung and Gut
psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

LungPrevotellaShared <- subset_taxa(psLungSharedRel, Genus == "Prevotella")
GutPrevotellaShared <- subset_taxa(psGutSharedRel, Genus == "Prevotella")

# Agglomeration
psLungSharedRel_Spe <- tax_glom(LungPrevotellaShared, "Species", NArm = TRUE)
psGutSharedRel_Spe <- tax_glom(GutPrevotellaShared, "Species", NArm = TRUE)

# Prepare a df for subsequent plot
psLungSharedRel_Spe_df <- psmelt(psLungSharedRel_Spe)
str(psLungSharedRel_Spe_df)
psGutSharedRel_Spe_df <- psmelt(psGutSharedRel_Spe)
str(psGutSharedRel_Spe_df)

# Convert as required
psLungSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psLungSharedRel_Spe_df)

psGutSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psGutSharedRel_Spe_df)

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                          levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("LungPrevotella_sharedGutLung_Boxplot.pdf",
       width = 20, height = 14, dpi = 200, units = "cm")

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                         levels = c("Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("GutPrevotella_sharedGutLung_Boxplot.pdf",
       width = 20, height = 14, dpi = 200, units = "cm")

# Analysis of Enterobacter species in Lung and Gut
psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

LungEnterobacterShared <- subset_taxa(psLungSharedRel, Genus == "Enterobacter")
GutEnterobacterShared <- subset_taxa(psGutSharedRel, Genus == "Enterobacter")

# Agglomeration
psLungSharedRel_Spe <- tax_glom(LungEnterobacterShared, "Species", NArm = TRUE)
psGutSharedRel_Spe <- tax_glom(GutEnterobacterShared, "Species", NArm = TRUE)

# Prepare a df for subsequent plot
psLungSharedRel_Spe_df <- psmelt(psLungSharedRel_Spe)
str(psLungSharedRel_Spe_df)
psGutSharedRel_Spe_df <- psmelt(psGutSharedRel_Spe)
str(psGutSharedRel_Spe_df)

# Convert as required
psLungSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psLungSharedRel_Spe_df)

psGutSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psGutSharedRel_Spe_df)

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                          levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("LungEnterobacter_sharedGutLung_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                         levels = c("Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("GutPrevotella_sharedGutLung_Boxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Analysis of Enterobacteriaceae species in Lung and Gut
psLungSharedRel <- transform_sample_counts(psLungShared, count_to_rel_abund)
psGutSharedRel <- transform_sample_counts(psGutShared, count_to_rel_abund)

LungEnterobacterShared <- subset_taxa(psLungSharedRel, Family == "Enterobacteriaceae")
GutEnterobacterShared <- subset_taxa(psGutSharedRel, Family == "Enterobacteriaceae")

# Agglomeration
psLungSharedRel_Spe <- tax_glom(LungEnterobacterShared, "Species", NArm = TRUE)
psGutSharedRel_Spe <- tax_glom(GutEnterobacterShared, "Species", NArm = TRUE)

# Prepare a df for subsequent plot
psLungSharedRel_Spe_df <- psmelt(psLungSharedRel_Spe)
str(psLungSharedRel_Spe_df)
psGutSharedRel_Spe_df <- psmelt(psGutSharedRel_Spe)
str(psGutSharedRel_Spe_df)

# Convert as required
psLungSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psLungSharedRel_Spe_df)

psGutSharedRel_Spe_df %<>% 
  select(Phylum, Family, Genus, Species, Abundance, Filtered_Data_AK_Dec2020) %>% 
  convert(fct(Phylum, Family, Genus, Species)) %>% arrange(-(Abundance))
head(psGutSharedRel_Spe_df)

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psLungSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                          levels = c("Extubation", "Infection_D5", "Infection_D1", "Inclusion"))

psLungSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors4) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("LungEnterobacteriaceae_sharedGutLung_Boxplot.pdf",
       width = 15, height = 15, dpi = 200, units = "cm")

# Grouped boxplot - Lung
# Reverse Group levels for vertical plot
psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020 <- factor(psGutSharedRel_Spe_df$Filtered_Data_AK_Dec2020,
                                                         levels = c("Extubation", "Infection D5", "Infection D1", "Inclusion"))

psGutSharedRel_Spe_df %>%
  ggplot(aes(x = reorder(Species, Abundance, FUN=median), y=Abundance)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 2.0,
    outlier.stroke = 0,
    outlier.alpha = 0.5,
  ) +
  scale_fill_manual(values = InclInfectColors5) +
  ylim(c(0, 1.0)) +
  coord_flip() +
  labs(x = "Species", y = "Relative abundance") +
  theme_classic()

ggsave("GutPrevotella_sharedGutLung_Boxplot.pdf",
       width = 15, height = 12, dpi = 200, units = "cm")
```
