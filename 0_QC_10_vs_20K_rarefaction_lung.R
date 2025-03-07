# title: "PIPOVAP_AK_17/11/2021 based on initial script by EB-PNE 15/02/2021"
# date: "16.02.2022"
# Graphs for quality comparison of lung 10K and 20K rarefaction level data
# Output: Comparison tables for alpha and beta diversity
---

# Install packages
  
install.packages("phyloseq")

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

library("phyloseq")
install.packages("tidyverse")
library(tidyverse)
library(readr)
library(ggplot2)
library(hablar)
library(lubridate)
library(magrittr)
library(RColorBrewer)
library(backports)
library(ggpubr)
library(BiocManager)
library(scales)
library(openxlsx)
library(remotes)
library(rmarkdown)
library(VennDiagram)
library(svglite)
install.packages("nloptr")
library(nloptr)
install.packages("GUniFrac")
library(GUniFrac)
install.packages("GUniFrac")
library(reshape2)

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

library(tidyselect)


## Import data
setwd("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/Quality control for lung samples")
getwd()

list.files()


## Load data with rarefaction at 10.000
# Those data come from the "220203_pipovap_lung" file

lung_10k = readRDS("220203_lung_10k.rds") 
head(lung_10k)
View(phyloseq_to_df(lung_10k))

physeqLung_10k <- lung_10k
OTU_lung_10k <- otu_table(physeqLung_10k)
sam_lung_10k <- sample_data(physeqLung_10k)
tax_lung_10k <- tax_table(physeqLung_10k)
tree_lung_10k <- phy_tree(physeqLung_10k)
refseq_lung_10k <- refseq(physeqLung_10k)


# Load previous data with rarefaction level set at 20.000

# Load Lung Data with No merge, Rarefaction = 20000, No collapse
Lung_base_with_tree_VS_20201209 <- readRDS("L:/SHARE/PIPOVAP/metagenomics/20201209_PIPOVAP_lungs/DADA2/4_physeq/RDP_ezbiocloud201805.202012.nomerge/Bacteria_in_Kingdom/rarefaction_20000/no_collapse/base_with_tree.rds")
physeqLung_20k <- Lung_base_with_tree_VS_20201209
OTU_lung_20k <- otu_table(physeqLung1)
sam_lung_20k <- sample_data(physeqLung1)
tax_lung_20k <- tax_table(physeqLung1)
tree_lung_20k <- phy_tree(physeqLung1)
refseq_lung_20k <- refseq(physeqLung1)

# Load filtered Lung secondary data from AK - Dec2020 
SecDataLungFil <- read_delim("Secondary_data_Filtered_AK_20210122_previous_version.csv", 
                             ";", escape_double = FALSE, trim_ws = TRUE)
spec(SecDataLungFil)
class(SecDataLungFil)
head(SecDataLungFil)
colnames(SecDataLungFil)

View(SecDataLungFil)

View(sam_lung_20k)

# Vector of Columns to be converted into factors 
ColAsFactors <- colnames(SecDataLungFil)[c(1:19, 56:62)]

SecDataLungFilConv <- SecDataLungFil %>% hablar::convert(fct(all_of(ColAsFactors)))



## Lung data pre-processing

# Use phyloseq::sample_data function in view of preparing a phyloseq object
sam_lung2 <- sample_data(SecDataLungFilConv)

sam_lung2

# For consistency with the other elements of phyloseq object, use "sample_label"
# as sample_names

sample_names(sam_lung2) <- as.factor(sam_lung2$sample_label)
head(sample_names(sam_lung2))

# Create new phyloseq object with sam_lung2 for 20k
physeqLung20k <- phyloseq(OTU_lung_20k, tax_lung_20k, sam_lung2, refseq_lung_20k, tree_lung_20k)
head(sample_names(physeqLung20k))
sample_variables(physeqLung20k)

# Correct the order of levels for different variables if necessary
levels(get_variable(physeqLung20k, "Filtered_Data_AK_Dec2020"))

Correct_order_Incl_Inf <- c("Inclusion", "Infect1D1", "Infect1D5", "Extub", "Excluded")

# Edit names of the levels
Correct_labels <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation",
                    "Excluded")

sample_data(physeqLung20k)$Filtered_Data_AK_Dec2020 <- factor(sample_data(physeqLung20k)$Filtered_Data_AK_Dec2020,
                                                            levels = Correct_order_Incl_Inf,
                                                            labels = Correct_labels)
levels(get_variable(physeqLung20k, "Filtered_Data_AK_Dec2020"))

# Set a manual scale of colors for Inclusion-Infection-Extubation groups
InclInfectColors4 <- c("Inclusion"="gold2", "Infection_D1"="coral2",
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3")

# Repeat as above for 10k metadata

# Create new phyloseq object with sam_lung2 for 10k
physeqLung10k <- phyloseq(OTU_lung_10k, tax_lung_10k, sam_lung2, refseq_lung_10k, tree_lung_10k)
head(sample_names(physeqLung10k))
sample_variables(physeqLung10k)


#Alpha diversity comparison  between 10k and 20k rarefaction levels
#calculate alpha diversity for 10k level
psLungAlphaRich_df_10k <- psmelt(physeqLung_10k)
View(psLungAlphaRich_df_10k)

psLungAlphaRich_df_10k %<>% 
  select(OTU, Sample, Abundance, Phylum,
         Shannon, InvSimpson,Chao1, Observed) %>% 
  convert(fct(OTU, Sample, Phylum)) %>% arrange(-(Abundance))

chao1_10k<-psLungAlphaRich_df_10k %>%
  ggplot(aes(x ="" , y = Chao1)) +
  geom_boxplot()+
  labs(x = "rarefaction 10k", y = "Chao1")
  
chao1_10k

shannon_10k<-psLungAlphaRich_df_10k %>%
  ggplot(aes(x ="" , y = Shannon)) +
  geom_boxplot()+
  labs(x = "rarefaction 10k", y = "Shannon")

shannon_10k


# Repeat the same procedure for 20k rarefaction level

psLungAlphaRich_df_20k <- psmelt(physeqLung_20k)
View(psLungAlphaRich_df_20k)

psLungAlphaRich_df_20k %<>% 
  select(OTU, Sample, Abundance, Phylum,
         Shannon, InvSimpson,Chao1, Observed) %>% 
  convert(fct(OTU, Sample, Phylum)) %>% arrange(-(Abundance))

chao1_20k<-psLungAlphaRich_df_20k %>%
  ggplot(aes(x ="" , y = Chao1)) +
  geom_boxplot()+
  labs(x = "rarefaction 20k", y = "Chao1")

chao1_20k

shannon_20k<-psLungAlphaRich_df_20k %>%
  ggplot(aes(x ="" , y = Shannon)) +
  geom_boxplot()+
  labs(x = "rarefaction 20k", y = "Shannon")

shannon_20k

  
#Combine plots
ggarrange(chao1_10k,chao1_20k,ncol=2,nrow=1,
                            common.legend = T,legend="right")

ggsave("alpha_div_comparion_10k_vs_20k_chao1.tiff")

ggarrange(shannon_10k,shannon_20k,ncol=2,nrow=1,
          common.legend = T,legend="right")
ggsave("alpha_div_comparion_10k_vs_20k_shannon.tiff")



## Beta diversity comparison between 10k and 20k
View(psmelt(physeqLung10k))
View(psmelt(physeqLung20k))

levels(get_variable(physeqLung10k, "sample_group"))
levels(get_variable(physeqLung20k, "sample_group"))


readcounts <- t(physeqLung10k@otu_table)
View(readcounts)

iDist <- vegdist(readcounts, method = "jaccard", binary = T)
view(iDist)

paired_dist <- melt(as.matrix(iDist), varnames = c("row", "col"))

View(paired_dist)

meta <- meta(physeqLung10k)
View(meta)
metakey <-data.frame(meta[c("Group", "record_ID")])
metakey


row<-rownames(meta)
row
col<-rownames(meta)
col
rowpatient<-meta$record_ID
colpatient<-meta$record_ID

metakeyrow <-data.frame(row, rowpatient)
metakeycol <-data.frame(col, colpatient)

metakeyrow
metakeycol

paired_dist_r <-merge(paired_dist, metakeyrow, by="row", suffixes = "row")
paired_dist_rc <-merge(paired_dist_r, metakeycol, by="col", suffixes = "col")

View(paired_dist_r)
View(paired_dist_rc)

ggplot(paired_dist_rc,aes(x=colpatient,y=value))+
  geom_jitter()


#ad empty column to data.frame "paired_dist_rc"
paired_dist_rc_10<-paired_dist_rc
paired_dist_rc_10$code<-"code"
paired_dist_rc_10

#Paste col and row columns to create a unique code for every ASV comparison
paired_dist_rc_10$code<-paste(paired_dist_rc_10$col,paired_dist_rc_10$row)


paired_dist_rc_10k<-paired_dist_rc


## Repeat the procedure above to create a distance matrix for the 20k phyloseq object
readcounts <- t(physeqLung20k@otu_table)

readcounts

iDist <- vegdist(readcounts, method = "jaccard", binary = T)


paired_dist <- melt(as.matrix(iDist), varnames = c("row", "col"))
meta <- meta(physeqLung20k)

metakey <-data.frame(meta[c("Group", "record_ID")])

row<-rownames(meta)
col<-rownames(meta)
rowpatient<-meta$record_ID
colpatient<-meta$record_ID

metakeyrow <-data.frame(row, rowpatient)
metakeycol <-data.frame(col, colpatient)

paired_dist_r <-merge(paired_dist, metakeyrow, by="row", suffixes = "row")
paired_dist_rc <-merge(paired_dist_r, metakeycol, by="col", suffixes = "col")



ggplot(paired_dist_rc,aes(x=colpatient,y=value))+
  geom_jitter()

View(paired_dist_rc)

paired_dist_rc_20<-paired_dist_rc


#ad empty column to data.frame "paired_dist_rc"
paired_dist_rc_20$code<-"code20"

#Paste col and row columns to create a unique code for every ASV comparison
paired_dist_rc_20$code<-paste(paired_dist_rc_20$col,paired_dist_rc_20$row)



dist20<-select(paired_dist_rc_20,code,value,colpatient)
dist20<-dplyr::rename(dist20,code=code20,value20=value)


dist10<-select(paired_dist_rc_10,code,value,colpatient)

#Join the two tables together
dist<-left_join(dist10,dist20,by="code")

View(dist)

#Create a new column and calculate the difference of beta diversity among two datasets
View(dist)
dist$divdiff<-"diff"
dist$divdiff<-(dist$value-dist$value20)

#Plot the difference of beta diversity
a<-ggplot(dist,aes(x=colpatient.x,y=divdiff))+
  geom_jitter()+
  labs(title="Lung Beta diversity difference between 10k and 20k rarefied datasets",x="record_ID",y="Beta diversity difference (absolute number)")

ggsave("beta_div_diff_lung.tiff",width = 35, height = 25, dpi = 300, units = "cm")

b<-a+geom_smooth()
b
ggsave("beta_div_diff_lung_with_smooth_line.tiff",width = 35, height = 25, dpi = 300, units = "cm")


dist$colpatient.x<-as.numeric(dist$colpatient.x)

##Plot beta diversity of two datasets in the same graph
# A little more pre-processing. Our two datasets have to have the same column names in order to combine them and plot them later on

dist10
dist20<-dplyr::rename(dist20,value=value20)

#Ad columns to specify whic datasets correspond to each rarefaction level
dist10$rarefaction_level<-"10k"
dist20$rarefaction_level<-"20k"

all_dist<-rbind(dist10,dist20)

#Plot beta-diversity of both datasets
ggplot(all_dist,aes(x=colpatient,y=value,color=rarefaction_level))+
  geom_jitter(alpha=0.3)+
  labs(title="Lung Beta diversity (10k and 20k rarefied datasets)",x="record_ID",y="Beta diversity value")


ggsave("beta_div_lung.tiff",width = 35, height = 25, dpi = 300, units = "cm")

################################

a<-ggplot(all_dist,aes(x=colpatient,y=value,col=rarefaction_level))+
  geom_jitter(alpha=0.05)
  
  a+stat_smooth(color = rarefaction_level, fill = "#FC4E07",
                method = "loess")
  
a


#############################################
# Beta- diversity ordination plots comparison

## Load data with rarefaction at 10.000
# Those data come from the "220203_pipovap_lung" file

lung_10k = readRDS("220203_lung_10k.rds") 

physeqLung_10k <- lung_10k
OTU_lung_10k <- otu_table(physeqLung_10k)
sam_lung_10k <- sample_data(physeqLung_10k)
tax_lung_10k <- tax_table(physeqLung_10k)
tree_lung_10k <- phy_tree(physeqLung_10k)
refseq_lung_10k <- refseq(physeqLung_10k)

sample_names(sam_lung_10k)<-paste(sample_names(sam_lung_10k),"_10k",sep="")
sample_names(sam_lung_10k)
variable.names(sam_lung_10k)
View(melt(sam_lung_10k))

sam.new <- data.frame(rarefaction = sample("10k", size = nsamples(physeqLung_10k), replace = TRUE))
# Mix up the sample names (just for demonstration purposes)
rownames(sam.new) <- sample(sample_names(physeqLung_10k))
# Turn into `sample_data` 
sam.new <- sample_data(sam.new)
head(sam.new)

#merge with initial phyloseq object
physeqLung_10k <- merge_phyloseq(physeqLung_10k, sam.new)
head(sample_data(physeqLung_10k))
variable.names(sample_data(physeqLung_10k))
sample_names(physeqLung_10k)
sample_names(physeqLung_10k)<-paste(sample_names(physeqLung_10k),"_10k",sep="")

physeqLung_10k 
OTU_lung_10k <- otu_table(physeqLung_10k)
sam_lung_10k <- sample_data(physeqLung_10k)
tax_lung_10k <- tax_table(physeqLung_10k)
tree_lung_10k <- phy_tree(physeqLung_10k)
refseq_lung_10k <- refseq(physeqLung_10k)

variable.names(sam_lung_10k)
sample_names(sam_lung_10k)



# Load previous data with rarefaction level set at 20.000

# Load Lung Data with No merge, Rarefaction = 20000, No collapse
Lung_base_with_tree_VS_20201209 <- readRDS("L:/SHARE/PIPOVAP/metagenomics/20201209_PIPOVAP_lungs/DADA2/4_physeq/RDP_ezbiocloud201805.202012.nomerge/Bacteria_in_Kingdom/rarefaction_20000/no_collapse/base_with_tree.rds")
physeqLung_20k <- Lung_base_with_tree_VS_20201209
OTU_lung_20k <- otu_table(physeqLung1)
sam_lung_20k <- sample_data(physeqLung1)
tax_lung_20k <- tax_table(physeqLung1)
tree_lung_20k <- phy_tree(physeqLung1)
refseq_lung_20k <- refseq(physeqLung1)

sample_names(physeqLung_20k)
variable.names(sam_lung_20k)

sam.new <- data.frame(rarefaction = sample("20k", size = nsamples(physeqLung_20k), replace = TRUE))
# Mix up the sample names (just for demonstration purposes)
rownames(sam.new) <- sample(sample_names(physeqLung_20k))
# Turn into `sample_data` 
sam.new <- sample_data(sam.new)
head(sam.new)

#merge with initial phyloseq object
physeqLung_20k <- merge_phyloseq(physeqLung_20k, sam.new)
head(sample_data(physeqLung_20k))
variable.names(sample_data(physeqLung_20k))
sample_names(physeqLung_20k)
sample_names(physeqLung_20k)<-paste(sample_names(physeqLung_20k),"_20k",sep="")

physeqLung_20k 
OTU_lung_20k <- otu_table(physeqLung_20k)
sam_lung_20k <- sample_data(physeqLung_20k)
tax_lung_20k <- tax_table(physeqLung_20k)
tree_lung_20k <- phy_tree(physeqLung_20k)
refseq_lung_20k <- refseq(physeqLung_20k)

variable.names(sam_lung_20k)
sample_names(sam_lung_20k)

OTU<-merge_phyloseq(OTU_lung_10k,OTU_lung_20k)
OTU
sam<-merge_phyloseq(sam_lung_10k,sam_lung_20k)
sam
tax<-merge_phyloseq(tax_lung_10k,tax_lung_20k)
tax
refseq<-merge_phyloseq(refseq_lung_10k,refseq_lung_20k)


phyloseqlung<-phyloseq(OTU,sam,tax,refseq)
phyloseqlung

View(psmelt(phyloseqlung))

variable.names(sample_data(phyloseqlung))

readcounts <- t(phyloseqlung@otu_table)

dist_bray<-vegdist(readcounts, method = "bray",na.rm=TRUE)

k<-rowSums(readcounts)>0
readcounts<-readcounts[k,]

dist_bray<-vegdist(readcounts, method = "bray",na.rm=TRUE)
dist_jac<-vegdist(readcounts, method = "jaccard",na.rm=TRUE)


dist <- phyloseq::distance(phyloseqlung, method="jaccard",na.rm=TRUE)
psmelt(dist)

phyloseqlung_ord_bray <- ordinate(phyloseqlung, "NMDS", dist_bray)
phyloseqlung_ord_jac <- ordinate(phyloseqlung, "NMDS", dist_jac)


bray<-plot_ordination(phyloseqlung, phyloseqlung_ord_bray, type="samples",color="sample_label",shape="rarefaction") +
  theme_bw()+
  geom_point(size=5)+
  labs(title="NMDS:Bray-Curtis distance")

ggsave("beta_div_10k_vs_20k_lung_bray.tiff",width = 35, height = 25, dpi = 300, units = "cm")

bray+theme(legend.position="")
ggsave("beta_div_10k_vs_20k_lung_bray_nolegend.tiff",width = 35, height = 25, dpi = 300, units = "cm")


jaccard<-plot_ordination(phyloseqlung, phyloseqlung_ord_jac, type="samples",color="sample_label",shape="rarefaction") +
  theme_bw()+
  geom_point(size=5)+
  labs(title="NMDS:Jaccard distance")

ggsave("beta_div_10k_vs_20k_lung_jaccard.tiff",width = 35, height = 25, dpi = 300, units = "cm")

jaccard+theme(legend.position = "")

ggsave("beta_div_10k_vs_20k_lung_jaccard_nolegend.tiff",width = 35, height = 25, dpi = 300, units = "cm")







## Beta diversity comparison between 10k and 20k
View(psmelt(phyloseqlung))


levels(get_variable(physeqLung10k, "sample_group"))
levels(get_variable(physeqLung20k, "sample_group"))


readcounts <- t(phyloseqlung@otu_table)
View(readcounts)

iDist <- vegdist(readcounts, method = "jaccard", binary = T)


paired_dist <- melt(as.matrix(iDist), varnames = c("row", "col"))

View(paired_dist)

meta <- meta(phyloseqlung)
View(meta)
metakey<-data.frame(meta[c("sample_label","rarefaction")])
metakey


row<-rownames(meta)
row
col<-rownames(meta)
col
rowsample<-meta$sample_label
colsample<-meta$sample_label
rowrar<-meta$rarefaction
colrar<-meta$rarefaction

metakeyrow <-data.frame(row, rowsample,rowrar)
metakeycol <-data.frame(col, colsample,colrar)

metakeyrow
metakeycol

paired_dist_r <-merge(paired_dist, metakeyrow, by="row", suffixes = "row")
paired_dist_rc <-merge(paired_dist_r, metakeycol, by="col", suffixes = "col")

View(paired_dist_r)
View(paired_dist_rc)



ggplot(paired_dist_rc,aes(x=colpatient,y=value))+
  geom_jitter()