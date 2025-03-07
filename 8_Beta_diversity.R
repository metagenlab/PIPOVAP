# title: "PIPOVAP_AK_04/05/22"
# date: "04.05.22"
# Venn diagrams for shared species among different timepoints for the 3 groups of patients
# Output: Venn diagram figures 
---
  
  # Install packages
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(BiocManager)
BiocManager::install("phyloseq")
#install.packages("phyloseq")
library("phyloseq")
BiocManager::install("GenomeInfoDb",force=TRUE)
library(GenomeInfoDb)
BiocManager::install("BiocParallel")
library(BiocParallel)
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
install.packages("VennDiagram")
library(VennDiagram)
install.packages("svglite")
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
install.packages("ampvis2")
library(ampvis2)
BiocManager::install("DESeq2")
install.packages("DESeq2")
library(DESeq2)
install.packages("readxl")
library("readxl")
BiocManager::install("miaViz")
library(mia)
library(miaViz)
library(tidyselect)
devtools::install_github("gmteunisse/fantaxtic")
# install.packages("fantaxtic")
library(fantaxtic)
remotes::install_github("gmteunisse/Fantaxtic")
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
install.packages("metacoder")
library(metacoder)
remotes::install_github("kstagaman/phyloseqCompanion")
install.packages("ggplotify")
library("ggplotify")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library("pairwiseAdonis")
#####################################################################################################################


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

View(sam_lung1)

# Load filtered Lung metadata 
metadata_lung <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/8.Beta diversity/metadata_lung.xlsx",sheet = "Feuil1")
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


## Lung data preprocessing ----

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

View(sam_lung2)

# Correct the order of levels for different variables if necessary
levels(factor(get_variable(physeqLung2, "Filtered_Data_AK_Dec2020")))

levels(factor(get_variable(physeqLung2, "Subgroup")))

levels(factor(get_variable(physeqLung2, "time")))

Correct_order_filtered_data <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation")
Correct_order_subgroup <- c("Inclusion","Infection_D5", "Discharge")
Correct_order_time <- c("Inclusion", "Infection_D5","Discharge")


# Edit names of the levels
Correct_labels <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation")

sample_data(physeqLung2)$Filtered_Data_AK_Dec2020 <- factor(sample_data(physeqLung2)$Filtered_Data_AK_Dec2020,
                                                            levels = Correct_order_filtered_data,
                                                            labels = Correct_labels)
levels(get_variable(physeqLung2, "Filtered_Data_AK_Dec2020"))

sample_data(physeqLung2)$Subgroup <- factor(sample_data(physeqLung2)$Subgroup)

levels(get_variable(physeqLung2,"Subgroup"))

sample_data(physeqLung2)$time <- factor(sample_data(physeqLung2)$time,
                                        levels=Correct_order_time)

levels(get_variable(physeqLung2,"time"))

## Keep only taxa with total abundance at least 2 in at least 1 samples
keep_function <- function(x) {x >= 1}
TaxaToKeep <- genefilter_sample(physeqLung2, keep_function, A = 1)

TaxaToKeep
TaxaToKeep[1:20]

physeqLung3 <- prune_taxa(TaxaToKeep, physeqLung2)

sort(taxa_sums(physeqLung3), decreasing = F)[1:50]


View(sample_data(physeqLung3))
#####################################################################################################################

# Load Gut Data with No merge, Rarefaction = 100000, No collapse
Gut_base_with_tree <- readRDS("//file3.intranet.chuv/data3/SHARE/PIPOVAP/220131_pipovap_gut/DADA2/4_physeq/RDP_ezbiocloud201805.202011/Bacteria_in_Kingdom/rarefaction_100000/no_collapse/base_with_tree_normNONE_abund0_prev0.rds")
physeqGut1 <- Gut_base_with_tree
OTU_Gut1 <- otu_table(physeqGut1)
sam_Gut1 <- sample_data(physeqGut1)
tax_Gut1 <- tax_table(physeqGut1)
tree_Gut1 <- phy_tree(physeqGut1)
refseq_Gut1 <- refseq(physeqGut1)

View(sam_Gut1)

head(rownames(OTU_Gut1))
head(colnames(OTU_Gut1))

head(rownames(sam_Gut1))
head(colnames(sam_Gut1))

head(rownames(tax_Gut1))
head(colnames(tax_Gut1))

head(sam_Gut1)
variable.names(sam_Gut1)


# Load complemented Gut secondary data
metadata_gut <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/8.Beta diversity/metadata_gut.xlsx",sheet = "working_sheet")
View(metadata_gut)

spec(metadata_gut)
head(metadata_gut)

colnames(metadata_gut)

# Vector of Columns to be converted into factors 
ColAsFactors <- colnames(metadata_gut)[c(1:20)]

metadata_gut <- metadata_gut %>% hablar::convert(fct(all_of(ColAsFactors)))

head(metadata_gut)


## Gut data preprocessing ----

# Use phyloseq::sample_data function in view of preparing a phyloseq object
sam_Gut2 <- sample_data(metadata_gut)
View(sam_Gut2)

# For consistency with the other elements of phyloseq object, use "sample_label"
# as sample_names

sample_names(sam_Gut2) <- as.factor(sam_Gut2$Sample)
head(sample_names(sam_Gut2))

# Create new phyloseq object with sam_Gut2
physeqGut2 <- phyloseq(OTU_Gut1, tax_Gut1, sam_Gut2, refseq_Gut1, tree_Gut1)
head(sample_names(physeqGut2))
sample_variables(physeqGut2)

# Correct the order of levels for different variables if necessary
levels(factor(get_variable(physeqGut2, "Filtered_Data_AK_Dec2020")))

levels(factor(get_variable(physeqGut2, "Subgroup")))

levels(factor(get_variable(physeqGut2, "time")))

Correct_order_filtered_data <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation","Discharge")

Correct_order_time <- c("Inclusion", "Infection_D5","Discharge")

# Edit names of the levels
Correct_labels <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation","Discharge")

sample_data(physeqGut2)$Filtered_Data_AK_Dec2020 <- factor(sample_data(physeqGut2)$Filtered_Data_AK_Dec2020,
                                                           levels = Correct_order_filtered_data,
                                                           labels = Correct_labels)
levels(get_variable(physeqGut2, "Filtered_Data_AK_Dec2020"))

sample_data(physeqGut2)$Subgroup <- factor(sample_data(physeqGut2)$Subgroup)

levels(get_variable(physeqGut2,"Subgroup"))

sample_data(physeqGut2)$time<-factor(sample_data(physeqGut2)$time,
                                    levels=Correct_order_time)

levels(get_variable(physeqGut2,"time"))


## Keep only taxa with total abundance at least 2 in at least 1 samples
keep_function <- function(x) {x >= 1}
TaxaToKeep <- genefilter_sample(physeqGut2, keep_function, A = 1)
TaxaToKeep[1:20]

physeqGut3 <- prune_taxa(TaxaToKeep, physeqGut2)

sort(taxa_sums(physeqGut3), decreasing = F)[1:50]

View(sample_data(physeqGut3))
#####################################################################################################################


### Set color palettes for further graphical representation

# Set manual scales of colors for the following groups: 1)Filtered_data_AK_Dec2020, 2)Subgroup, 3)Group

#Color palette for groups comparison
Infectcolorgroup<-c("Control"="#80b2d3","Other_infection"="#d38088","Pneumonia"="#d3cc80")  # color palette for Groups

# Color palette for Subgroups comparison
InclInfectColorsubgroup <- c("Control T1"="#80b2d3", "Control T3"="#80d3cb",
                             "OtherInfection T1"="#d38088", "OtherInfection T2"="#d380b2","OtherInfection T3"="#d3a280",
                             "Pneumonia T1"="#d3cc80","Pneumonia T2"="#b1d380","Pneumonia T3"="#c8e0a6")
                             

# Color palette for Timepoints for Lung
InclInfectColors4 <- c("Inclusion"="gold2", "Infection_D1"="coral2",  #color palette for Filtered_data_AK_Dec2020
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3")


# Color palette for timepoints for Gut
InclInfectColors5 <- c("Inclusion"="gold2", "Infection_D1"="coral2",
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3", "Discharge"="orchid1")




#####################################################################################################################


## Subset phyloseq objects to prepare phyloseq objects containing only "Inclusion" and "Extubation/Discharge" samples for Lung and Gut respectively

# First for Lung samples
ps1<- subset_samples(physeqLung3,Filtered_Data_AK_Dec2020=="Inclusion") # only inclusion timepoints
ps2<-subset_samples(physeqLung3,Filtered_Data_AK_Dec2020=="Extubation") # only extubation timepoints
psc<-subset_samples(physeqLung3,Group=="Control") #only Control group
psp<-subset_samples(physeqLung3,Group=="Pneumonia") #only Pneumonia group
pso<-subset_samples(physeqLung3,Group=="Other_infection") #only Other_infection group


# Then for Gut samples
pg1<- subset_samples(physeqGut3,Filtered_Data_AK_Dec2020=="Inclusion") # only inclusion timepoints
pg2<-subset_samples(physeqGut3,Time_point=="T3") # only Discharge timepoints
pgc<-subset_samples(physeqGut3,Group=="Control") #only Control group
pgp<-subset_samples(physeqGut3,Group=="Pneumonia") #only Pneumonia group
pgo<-subset_samples(physeqGut3,Group=="Other_infection") #only Other_infection group

#####################################################################################################################


getwd()
setwd("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/8.Beta diversity")



######   BETA Diversity measures    ########

# First calculate beta diversity for the lung samples

# We will fisrt produce a NMDS ordination plot only for inclusion timepoints and then including all other timepoints


##### BRAY DISTANCE   ###########3

set.seed(1)

# Ordinate
ps1_bray <- ordinate(
  physeq = ps1, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)

set.seed(1)
ps2_bray <- ordinate(
  physeq = ps2, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)

set.seed(1)
psc_bray <- ordinate(
  physeq = psc, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)

set.seed(1)
psp_bray <- ordinate(
  physeq = psp, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)


set.seed(1)
pso_bray <- ordinate(
  physeq = pso, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)



set.seed(1)
ps1_jaccard <- ordinate(
  physeq = ps1, 
  method = "NMDS", 
  distance = "jaccard", binary=TRUE) # using the Jaccard index

ps2_jaccard <- ordinate(
  physeq = ps2, 
  method = "NMDS", 
  distance = "jaccard", binary=TRUE) # using the Jaccard index



## Plot Ordination

## Only Inclusion timepoint

a1<-plot_ordination(
  physeq = ps1,
  ordination = ps1_bray,
  color = "Group",
  shape="Group",
  title = "Bray-curtis dissimilarity (inclusion)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

a1
  
## Calculate Centroids vegan for the plot
phybeta <- ps1
readcounts <- t(ps1@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

a1<-a1+geom_point(data = cent, size = 8, alpha=0.8)
a1

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps1, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps1)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps1)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)


#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
lung_incl_bray <- phyloseq::distance(ps1, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps1))

# Adonis test
adonis2(lung_incl_bray ~ Group, data = sampledf)

pairwise.adonis2(lung_incl_bray ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_incl_bray, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
a1<-a1+annotate("text",x=-1,y=c(1,0.9),label=c("betadisper ns","adonis ns"))
ggsave("bray_lung_incl.tiff",width = 6,height = 6,dpi=300)



## Only Extubation/Discharge timepoint

a2<-plot_ordination(
  physeq = ps2,
  ordination = ps2_bray,
  color = "Group",
  shape="Group",
  title = "Bray-curtis dissimilarity (extubation)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

a2

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- ps2
readcounts <- t(ps2@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

a2<-a2+geom_point(data = cent, size = 8, alpha=0.8)
a2

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps2, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps2)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps2)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
lung_extub_bray <- phyloseq::distance(ps2, method = "bray")
lung_extub_bray
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps2))

# Adonis test
adonis2(lung_extub_bray ~ Group, data = sampledf)

pairwise.adonis2(lung_extub_bray ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_extub_bray, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
a2<-a2+annotate("text",x=-0.7,y=c(0.6,0.5),label=c("betadisper ***","adonis ns"))
ggsave("bray_lung_extub.tiff",width = 6,height = 6,dpi=300)


## Only Controls, all timepoint

ac<-plot_ordination(
  physeq = psc,
  ordination = psc_bray,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "Bray-curtis dissimilarity (Control)")+  
  scale_color_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3"))+
  scale_fill_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3"))+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

ac

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- psc
readcounts <- t(psc@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

ac<-ac+geom_point(data = cent, size = 8, alpha=0.8)
ac

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(psc, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(psc)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(psc)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
lung_control_bray <- phyloseq::distance(psc, method = "bray")
lung_control_bray

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psc))

# Adonis test
adonis2(lung_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(lung_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_control_bray, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
ac<-ac+annotate("text",x=-1.1,y=c(0.6,0.5),label=c("betadisper ns","adonis ns"))
ggsave("bray_lung_control.tiff",width = 6,height = 6,dpi=300)





## Only Pneumonia, all timepoint

ap<-plot_ordination(
  physeq = psp,
  ordination = psp_bray,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "Bray-curtis dissimilarity (Pneumonia)")+  
  scale_color_manual(values = InclInfectColors4)+
  scale_fill_manual(values = InclInfectColors4)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()



## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- psp
readcounts <- t(psp@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

ap<-ap+geom_point(data = cent, size = 8, alpha=0.8)
ap

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(psp, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(psp)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(psp)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
lung_control_bray <- phyloseq::distance(psp, method = "bray")
lung_control_bray

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psp))

# Adonis test
adonis2(lung_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(lung_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_control_bray, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
ap<-ap+annotate("text",x=-5,y=c(3.6,3),label=c("betadisper ns","adonis ns"))
ggsave("bray_lung_pneumonia.tiff",width = 6,height = 6,dpi=300)




## Only Other indections, all timepoint

ao<-plot_ordination(
  physeq = pso,
  ordination = pso_bray,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "Bray-curtis dissimilarity (Other infection)")+  
  scale_color_manual(values = InclInfectColors4)+
  scale_fill_manual(values = InclInfectColors4)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

ao


## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pso
readcounts <- t(pso@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

ao<-ao+geom_point(data = cent, size = 8, alpha=0.8)
ao

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pso, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pso)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pso)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
lung_control_bray <- phyloseq::distance(pso, method = "bray")
lung_control_bray

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pso))

# Adonis test
adonis2(lung_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(lung_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_control_bray, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
ao+annotate("text",x=-1.1,y=c(0.9,0.8),label=c("betadisper *","adonis ns"))
ggsave("bray_lung_other_infection.tiff",width = 6,height = 6,dpi=300)


#############################################
#############################################
#############################################
#############################################
#############################################


# Now calculate beta diversity for the GUT samples

# We will fisrt produce a NMDS ordination plot only for inclusion timepoints and then including all other timepoints
set.seed(1)

# Ordinate
pg1_bray <- ordinate(
  physeq = pg1, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)

set.seed(1)
pg2_bray <- ordinate(
  physeq = pg2, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)

set.seed(1)
pgc_bray <- ordinate(
  physeq = pgc, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)

set.seed(1)
pgp_bray <- ordinate(
  physeq = pgp, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)


set.seed(1)
pgo_bray <- ordinate(
  physeq = pgo, 
  method = "NMDS", 
  distance = "bray"  # calculating distance with Bray-Curtis index
)



## Plot Ordination

## Only Inclusion timepoint

g1<-plot_ordination(
  physeq = pg1,
  ordination = pg1_bray,
  color = "Group",
  shape="Group",
  title = "Bray-curtis dissimilarity (inclusion)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

g1

View(meta(pg1))

## Calculate Centroids vegan for the plot
phybeta <- pg1
readcounts <- t(pg1@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

g1<-g1+geom_point(data = cent, size = 8, alpha=0.8)
g1

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pg1, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pg1)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pg1)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)


#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
gut_incl_bray <- phyloseq::distance(pg1, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pg1))

# Adonis test
adonis2(gut_incl_bray ~ Group, data = sampledf)

pairwise.adonis2(gut_incl_bray ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_incl_bray, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
g1<-g1+annotate("text",x=-0.5,y=c(0.6,0.5),label=c("betadisper ns","adonis ns"))
ggsave("bray_gut_incl.tiff",width = 6,height = 6,dpi=300)



## Only Extubation/Discharge timepoint

g2<-plot_ordination(
  physeq = pg2,
  ordination = pg2_bray,
  color = "Group",
  shape="Group",
  title = "Bray-curtis dissimilarity (discharge)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

g2

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pg2
readcounts <- t(pg2@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

g2<-g2+geom_point(data = cent, size = 8, alpha=0.8)
g2

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pg2, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pg2)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pg2)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
gut_discharge_bray <- phyloseq::distance(pg2, method = "bray")
gut_discharge_bray
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pg2))

# Adonis test
adonis2(gut_discharge_bray ~ Group, data = sampledf)

pairwise.adonis2(gut_discharge_bray ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_discharge_bray, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
g2<-g2+annotate("text",x=-0.55,y=c(0.65,0.55),label=c("betadisper ns","adonis ***"))
ggsave("bray_gut_discharge(c_vs_o_signif).tiff",width = 6,height = 6,dpi=300)




## Only Controls, all timepoint

gc<-plot_ordination(
  physeq = pgc,
  ordination = pgc_bray,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "Bray-curtis dissimilarity (Control)")+  
  scale_color_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3","Discharge"="orchid1"))+
  scale_fill_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3","Discharge"="orchid1"))+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

gc

InclInfectColors5 <- c("Inclusion"="gold2", "Infection_D1"="coral2",
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3", "Discharge"="orchid1")




## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pgc
readcounts <- t(pgc@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

gc<-gc+geom_point(data = cent, size = 8, alpha=0.8)
gc

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pgc, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pgc)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pgc)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
gut_control_bray <- phyloseq::distance(pgc, method = "bray")
gut_control_bray

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pgc))

# Adonis test
adonis2(gut_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(gut_control_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_control_bray, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
gc<-gc+annotate("text",x=-1.5,y=c(1.5,1.3),label=c("betadisper ns","adonis ns"))
ggsave("bray_gut_control.tiff",width = 6,height = 6,dpi=300)




## Only Pneumonia, all timepoint

gp<-plot_ordination(
  physeq = pgp,
  ordination = pgp_bray,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "Bray-curtis dissimilarity (Pneumonia)")+  
  scale_color_manual(values = InclInfectColors5)+
  scale_fill_manual(values = InclInfectColors5)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

gp

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pgp
readcounts <- t(pgp@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

gp<-gp+geom_point(data = cent, size = 8, alpha=0.8)
gp

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pgp, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pgp)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pgp)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
gut_pneumonia_bray <- phyloseq::distance(pgp, method = "bray")
gut_pneumonia_bray

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pgp))

# Adonis test
adonis2(gut_pneumonia_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(gut_pneumonia_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_pneumonia_bray, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
gp<-gp+annotate("text",x=-1.3,y=c(1.6,1.4),label=c("betadisper **","adonis *"))
ggsave("bray_gut_pneumonia.tiff",width = 6,height = 6,dpi=300)




## Only Other indections, all timepoint

go<-plot_ordination(
  physeq = pgo,
  ordination = pgo_bray,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "Bray-curtis dissimilarity (Other infection)")+  
  scale_color_manual(values = InclInfectColors5)+
  scale_fill_manual(values = InclInfectColors5)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

go


## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pgo
readcounts <- t(pgo@otu_table)
iDist <- vegdist(readcounts, method = "bray")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

go<-go+geom_point(data = cent, size = 8, alpha=0.8)
go

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pgo, method = "bray") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pgo)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pgo)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate bray curtis distance matrix
gut_other_bray <- phyloseq::distance(pgo, method = "bray")
gut_other_bray

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pgo))

# Adonis test
adonis2(gut_other_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(gut_other_bray ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_other_bray, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
go<-go+annotate("text",x=-2.1,y=c(1.2,1),label=c("betadisper ns","adonis ns"))
ggsave("bray_gut_other_infection.tiff",width = 6,height = 6,dpi=300)





##########################################################################################
##########################################################################################
########## REPEAT ALL OF THE ABOVE ANALYSES FOR JACCARD DISTANCE METRICS  ##################
##########################################################################################
##########################################################################################




######   BETA Diversity measures    ########

# First calculate beta diversity for the lung samples

# We will fisrt produce a NMDS ordination plot only for inclusion timepoints and then including all other timepoints


   ###  JACCARD DISTANCE  ###

set.seed(1)

# Ordinate
ps1_jaccard <- ordinate(
  physeq = ps1, 
  method = "NMDS", 
  distance = "jaccard", binary=TRUE)  # calculating distance with jaccard index


set.seed(1)
ps2_jaccard <- ordinate(
  physeq = ps2, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)

set.seed(1)
psc_jaccard <- ordinate(
  physeq = psc, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)

set.seed(1)
psp_jaccard <- ordinate(
  physeq = psp, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)


set.seed(1)
pso_jaccard <- ordinate(
  physeq = pso, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)





## Plot Ordination

## Only Inclusion timepoint

a1<-plot_ordination(
  physeq = ps1,
  ordination = ps1_jaccard,
  color = "Group",
  shape="Group",
  title = "jaccard dissimilarity (inclusion)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

a1

## Calculate Centroids vegan for the plot
phybeta <- ps1
readcounts <- t(ps1@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

a1<-a1+geom_point(data = cent, size = 8, alpha=0.8)
a1

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps1, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps1)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps1)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)


#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
lung_incl_jaccard <- phyloseq::distance(ps1, method = "jaccard")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps1))

# Adonis test
adonis2(lung_incl_jaccard ~ Group, data = sampledf)

pairwise.adonis2(lung_incl_jaccard ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_incl_jaccard, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
a1<-a1+annotate("text",x=-1,y=c(1,0.9),label=c("betadisper *","adonis ns"))
ggsave("jaccard_lung_incl.tiff",width = 6,height = 6,dpi=300)



## Only Extubation/Discharge timepoint

a2<-plot_ordination(
  physeq = ps2,
  ordination = ps2_jaccard,
  color = "Group",
  shape="Group",
  title = "jaccard dissimilarity (extubation)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

a2

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- ps2
readcounts <- t(ps2@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

a2<-a2+geom_point(data = cent, size = 8, alpha=0.8)
a2

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps2, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps2)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps2)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
lung_extub_jaccard <- phyloseq::distance(ps2, method = "jaccard")
lung_extub_jaccard
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps2))

# Adonis test
adonis2(lung_extub_jaccard ~ Group, data = sampledf)

pairwise.adonis2(lung_extub_jaccard ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_extub_jaccard, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
a2<-a2+annotate("text",x=-0.7,y=c(0.6,0.5),label=c("betadisper ***","adonis ns"))
ggsave("jaccard_lung_extub.tiff",width = 6,height = 6,dpi=300)


## Only Controls, all timepoint

ac<-plot_ordination(
  physeq = psc,
  ordination = psc_jaccard,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "jaccard dissimilarity (Control)")+  
  scale_color_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3"))+
  scale_fill_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3"))+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

ac

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- psc
readcounts <- t(psc@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

ac<-ac+geom_point(data = cent, size = 8, alpha=0.8)
ac

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(psc, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(psc)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(psc)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
lung_control_jaccard <- phyloseq::distance(psc, method = "jaccard")
lung_control_jaccard

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psc))

# Adonis test
adonis2(lung_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(lung_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_control_jaccard, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
ac<-ac+annotate("text",x=-0.8,y=c(0.8,0.7),label=c("betadisper ns","adonis ns"))
ggsave("jaccard_lung_control.tiff",width = 6,height = 6,dpi=300)




## Only Pneumonia, all timepoint

ap<-plot_ordination(
  physeq = psp,
  ordination = psp_jaccard,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "jaccard dissimilarity (Pneumonia)")+  
  scale_color_manual(values = InclInfectColors4)+
  scale_fill_manual(values = InclInfectColors4)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

ap

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- psp
readcounts <- t(psp@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

ap<-ap+geom_point(data = cent, size = 8, alpha=0.8)
ap

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(psp, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(psp)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(psp)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
lung_control_jaccard <- phyloseq::distance(psp, method = "jaccard")
lung_control_jaccard

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psp))

# Adonis test
adonis2(lung_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(lung_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_control_jaccard, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
ap<-ap+annotate("text",x=-0.55,y=c(0.5,0.4),label=c("betadisper **","adonis ns"))
ggsave("jaccard_lung_pneumonia.tiff",width = 6,height = 6,dpi=300)




## Only Other indections, all timepoint

ao<-plot_ordination(
  physeq = pso,
  ordination = pso_jaccard,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "jaccard dissimilarity (Other infection)")+  
  scale_color_manual(values = InclInfectColors4)+
  scale_fill_manual(values = InclInfectColors4)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

ao


## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pso
readcounts <- t(pso@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

ao<-ao+geom_point(data = cent, size = 8, alpha=0.8)
ao

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pso, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pso)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pso)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
lung_control_jaccard <- phyloseq::distance(pso, method = "jaccard")
lung_control_jaccard

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pso))

# Adonis test
adonis2(lung_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(lung_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(lung_control_jaccard, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
ao<-ao+annotate("text",x=-0.7,y=c(1,0.9),label=c("betadisper ***","adonis ns"))
ggsave("jaccard_lung_other_infection.tiff",width = 6,height = 6,dpi=300)


#############################################
#############################################
#############################################
#############################################
#############################################


# Now calculate beta diversity for the GUT samples

# We will fisrt produce a NMDS ordination plot only for inclusion timepoints and then including all other timepoints
set.seed(1)

# Ordinate
pg1_jaccard <- ordinate(
  physeq = pg1, 
  method = "NMDS", 
  distance = "jaccard", binary=TRUE  # calculating distance with jaccard index
)

set.seed(1)
pg2_jaccard <- ordinate(
  physeq = pg2, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)

set.seed(1)
pgc_jaccard <- ordinate(
  physeq = pgc, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)

set.seed(1)
pgp_jaccard <- ordinate(
  physeq = pgp, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)


set.seed(1)
pgo_jaccard <- ordinate(
  physeq = pgo, 
  method = "NMDS", 
  distance = "jaccard",binary=TRUE  # calculating distance with jaccard index
)



## Plot Ordination

## Only Inclusion timepoint

g1<-plot_ordination(
  physeq = pg1,
  ordination = pg1_jaccard,
  color = "Group",
  shape="Group",
  title = "jaccard dissimilarity (inclusion)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

g1

View(meta(pg1))

## Calculate Centroids vegan for the plot
phybeta <- pg1
readcounts <- t(pg1@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

g1<-g1+geom_point(data = cent, size = 8, alpha=0.8)
g1

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pg1, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pg1)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pg1)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)


#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
gut_incl_jaccard <- phyloseq::distance(pg1, method = "jaccard")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pg1))

# Adonis test
adonis2(gut_incl_jaccard ~ Group, data = sampledf)

pairwise.adonis2(gut_incl_jaccard ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_incl_jaccard, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
g1<-g1+annotate("text",x=-0.3,y=c(0.4,0.3),label=c("betadisper *","adonis ns"))
ggsave("jaccard_gut_incl.tiff",width = 6,height = 6,dpi=300)



## Only Extubation/Discharge timepoint

g2<-plot_ordination(
  physeq = pg2,
  ordination = pg2_jaccard,
  color = "Group",
  shape="Group",
  title = "jaccard dissimilarity (discharge)") + 
  scale_color_manual(values = Infectcolorgroup)+
  scale_fill_manual(values = Infectcolorgroup)+
  geom_point(aes(color = Group), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Group),alpha=0.15)+
  theme_bw()

g2


## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pg2
readcounts <- t(pg2@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Group, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Group=meta$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Group, data = scrs, FUN = median)

g2<-g2+geom_point(data = cent, size = 8, alpha=0.8)
g2

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pg2, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pg2)$Group)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pg2)$Group)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
gut_discharge_jaccard <- phyloseq::distance(pg2, method = "jaccard")
gut_discharge_jaccard
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pg2))

# Adonis test
adonis2(gut_discharge_jaccard ~ Group, data = sampledf)

pairwise.adonis2(gut_discharge_jaccard ~ Group, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_discharge_jaccard, sampledf$Group)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
g2<-g2+annotate("text",x=-0.55,y=c(0.5,0.4),label=c("betadisper ns","adonis **"))
g2
ggsave("jaccard_gut_discharge(c_vs_o_signif).tiff",width = 6,height = 6,dpi=300)




## Only Controls, all timepoint

gc<-plot_ordination(
  physeq = pgc,
  ordination = pgc_jaccard,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "jaccard dissimilarity (Control)")+  
  scale_color_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3","Discharge"="orchid1"))+
  scale_fill_manual(values = c("Inclusion"="gold2","Extubation"="cadetblue3","Discharge"="orchid1"))+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

gc

InclInfectColors5 <- c("Inclusion"="gold2", "Infection_D1"="coral2",
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3", "Discharge"="orchid1")




## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pgc
readcounts <- t(pgc@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

gc<-gc+geom_point(data = cent, size = 8, alpha=0.8)
gc

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pgc, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pgc)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pgc)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
gut_control_jaccard <- phyloseq::distance(pgc, method = "jaccard")
gut_control_jaccard

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pgc))

# Adonis test
adonis2(gut_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(gut_control_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_control_jaccard, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
gc<-gc+annotate("text",x=-2.2,y=c(1.5,1.3),label=c("betadisper *","adonis ns"))
ggsave("jaccard_gut_control.tiff",width = 6,height = 6,dpi=300)




## Only Pneumonia, all timepoint

gp<-plot_ordination(
  physeq = pgp,
  ordination = pgp_jaccard,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "jaccard dissimilarity (Pneumonia)")+  
  scale_color_manual(values = InclInfectColors5)+
  scale_fill_manual(values = InclInfectColors5)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

gp

## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pgp
readcounts <- t(pgp@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

gp<-gp+geom_point(data = cent, size = 8, alpha=0.8)
gp

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pgp, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pgp)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pgp)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
gut_pneumonia_jaccard <- phyloseq::distance(pgp, method = "jaccard")
gut_pneumonia_jaccard

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pgp))

# Adonis test
adonis2(gut_pneumonia_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(gut_pneumonia_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_pneumonia_jaccard, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
gp<-gp+annotate("text",x=-1.1,y=c(1.2,1),label=c("betadisper ***","adonis *"))
ggsave("jaccard_gut_pneumonia.tiff",width = 6,height = 6,dpi=300)




## Only Other indections, all timepoint

go<-plot_ordination(
  physeq = pgo,
  ordination = pgo_jaccard,
  color = "Filtered_Data_AK_Dec2020",
  shape="Filtered_Data_AK_Dec2020",
  title = "jaccard dissimilarity (Other infection)")+  
  scale_color_manual(values = InclInfectColors5)+
  scale_fill_manual(values = InclInfectColors5)+
  geom_point(aes(color = Filtered_Data_AK_Dec2020), alpha = 0.7, size = 3) +
  stat_ellipse(geom="polygon",aes(fill=Filtered_Data_AK_Dec2020),alpha=0.15)+
  theme_bw()

go


## Calculate Centroids vegan for the plot
set.seed(1)
phybeta <- pgo
readcounts <- t(pgo@otu_table)
iDist <- vegdist(readcounts, method = "jaccard")

# 2 subgroups
### Loop over distances, generate beta-diversity matrix and ordinate

nmds <- metaMDS(iDist)

# Caption
meta<-meta(phybeta)
readcounts <- t(phybeta@otu_table)
perma <- adonis2(readcounts~Filtered_Data_AK_Dec2020, data=meta , permutations = 999)
perma

scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs),Filtered_Data_AK_Dec2020=meta$Filtered_Data_AK_Dec2020)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Filtered_Data_AK_Dec2020, data = scrs, FUN = median)

go<-go+geom_point(data = cent, size = 8, alpha=0.8)
go

## Follows script to generate boxplots that show the distance from centroids

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(pgo, method = "jaccard") 

#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(pgo)$Filtered_Data_AK_Dec2020)

dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(pgo)$Filtered_Data_AK_Dec2020)
plot(dispr)
boxplot(dispr,main="",xlab="")
permutest(dispr)

#Calculate statistics

# Inclusion
set.seed(1)

# Calculate jaccard curtis distance matrix
gut_other_jaccard <- phyloseq::distance(pgo, method = "jaccard")
gut_other_jaccard

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(pgo))

# Adonis test
adonis2(gut_other_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)

pairwise.adonis2(gut_other_jaccard ~ Filtered_Data_AK_Dec2020, data = sampledf)


# Beta dispersion
beta<-betadisper(gut_other_jaccard, sampledf$Filtered_Data_AK_Dec2020)
beta
permutest(beta)
TukeyHSD(beta)


# Annotate graph with statistical results (adonis and betadisper)
go<-go+annotate("text",x=-2,y=c(1,0.8),label=c("betadisper **","adonis ns"))
ggsave("jaccard_gut_other_infection.tiff",width = 6,height = 6,dpi=300)


##########################################################################################
##########################################################################################
####################################   END OF SCRIPT #####################################
##########################################################################################
##########################################################################################
