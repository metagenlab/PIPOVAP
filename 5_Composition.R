# title: "PIPOVAP_AK_04.05.2022 Taxonomical Heat trees to compare groups composition"
# date: "04.05.2022"
# output: html_document
---
  
  # Libraries
  
 # if (!require("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")
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
install.packages("fantaxtic")
library(fantaxtic)
remotes::install_github("gmteunisse/Fantaxtic")
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
install.packages("metacoder")
library(metacoder)
remotes::install_github("kstagaman/phyloseqCompanion")

# Loading data ----


#####################################################################################################################


# Load Lung Data with No merge, Rarefaction = 10000, No collapse
physeqLung1<-readRDS("Users/yangjichoi/Documents/PhD/Projects/3_PIPOVAP/drive-download-20240307T145542Z-001/220203_lung_10k.rds")
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
metadata_lung <- read_excel("Users/yangjichoi/Documents/PhD/Projects/3_PIPOVAP/drive-download-20240307T145542Z-001/metadata_lung.xlsx")
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

Correct_order_filtered_data <- c("Inclusion", "Infect_D1", "Infect_D5", "Extubation")
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


#####################################################################################################################

# Load Gut Data with No merge, Rarefaction = 100000, No collapse
Gut_base_with_tree <- readRDS("Users/yangjichoi/Documents/PhD/Projects/3_PIPOVAP//drive-download-20240307T145542Z-001/base_with_tree_normNONE_abund0_prev0.rds")
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
metadata_gut <- read_excel("Users/yangjichoi/Documents/PhD/Projects/3_PIPOVAP//drive-download-20240307T145542Z-001/metadata_gut.xlsx")
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

View(psmelt(physeqGut2))

# Correct the order of levels for different variables if necessary
levels(get_variable(physeqGut2, "Filtered_Data_AK_Dec2020"))

levels(get_variable(physeqGut2, "Subgroup"))

Correct_order_filtered_data <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation","Discharge")
Correct_order_subgroup <- c("Inclusion","Infection_D5", "Discharge")

# Edit names of the levels
Correct_labels <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation","Discharge")

sample_data(physeqGut2)$Filtered_Data_AK_Dec2020 <- factor(sample_data(physeqGut2)$Filtered_Data_AK_Dec2020,
                                                            levels = Correct_order_filtered_data,
                                                            labels = Correct_labels)
levels(get_variable(physeqGut2, "Filtered_Data_AK_Dec2020"))

sample_data(physeqGut2)$Subgroup <- factor(sample_data(physeqGut2)$Subgroup,
                                            levels = Correct_order_subgroup)

levels(get_variable(physeqGut2,"Subgroup"))


#####################################################################################################################


### Set color palettes for further graphical representation

# Set manual scales of colors for the following groups: 1)Filtered_data_AK_Dec2020, 2)Subgroup, 3)Group

#Color palette for groups comparison
Infectcolorgroup<-c("Control"="#8DD3C7","Other_infection"="#FDB462","Pneumonia"="#80B1D3")  # color palette for Groups

# Color palette for Subgroups comparison
InclInfectColors3 <- c("Inclusion"="gold2", "Infection_D5"="chartreuse4",   # Color palette for Subgroup
                       "Discharge"="cadetblue3")


# Color palette for Timepoints for Lung
InclInfectColors4 <- c("Inclusion"="gold2", "Infection_D1"="coral2",  #color palette for Filtered_data_AK_Dec2020
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3")


# Color palette for tiepoints for Gut
InclInfectColors5 <- c("Inclusion"="gold2", "Infection_D1"="coral2",
                       "Infection_D5"="chartreuse4",
                       "Extubation"="cadetblue3", "Discharge"="orchid1")




#####################################################################################################################

## Subset phyloseq objects to prepare phyloseq objects containing only "Inclusion" and "Extubation/Discharge" samples for Lung and Gut respectively

# First for Lung samples
ps1<- subset_samples(physeqLung2,Filtered_Data_AK_Dec2020=="Inclusion") # only inclusion timepoints
ps2<-subset_samples(physeqLung2,Filtered_Data_AK_Dec2020=="Extubation") # only extubation timepoints

# Then for Gut samples
pg1<- subset_samples(physeqGut2,Filtered_Data_AK_Dec2020=="Inclusion") # only inclusion timepoints
pg2<-subset_samples(physeqGut2,Filtered_Data_AK_Dec2020=="Discharge") # only extubation timepoints


#####################################################################################################################

display.brewer.pal(12,"Set3")
brewer.pal(12,"Set3")

mypal=c(  "#D9D9D9", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#8DD3C7", "#BC80BD", "#CCEBC5", "#FFED6F")

mypal

## Create Barplots with the most abundant species over time (inclusion vs infection_d5 vs discharge)

## Create barplots of the 10 most abundant families. First for Lung samples

ps_tmp_fam<-microbiomeutilities::aggregate_top_taxa2(physeqLung2, top=10, "Genus")

ps_tmp_fam <- name_taxa(ps_tmp_fam, label = "Unkown", species= T, other_label = "Other")

ps_tmp_fam_lung<-subset_samples(ps_tmp_fam,Filtered_Data_AK_Dec2020=="Inclusion"|Filtered_Data_AK_Dec2020=="Extubation")

View(psmelt(ps_tmp_fam_lung))

ps_tmp_fam_lung<-ps_tmp_fam_lung%>%microbiome::transform(transform="compositional")

tax_table(ps_tmp_fam_lung)


############################################################################################################################
otu_df <- psmelt(ps_tmp_fam_lung)

# Normalize the data by sample depths
# Replace `size` with your abundance column if it's named differently
otu_df$normalized_abundance <- otu_df$Abundance / sample_sums(ps_tmp_fam_lung)[otu_df$Sample]

# Convert to long format for ggplot2
otu_df <- otu_df %>%
    group_by(Sample) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    ungroup()

# Create the stacked bar plot, faceting by both 'Group' and 'Inclusion/Extubation' status
p <- ggplot(otu_df, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values=mypal)+scale_color_manual(values=mypal) + # Adjust the palette as necessary
   # facet_grid(Filtered_Data_AK_Dec2020 ~ Group,  scales = "free") + # Assuming you have Inclusion_Extubation as a combined factor
 facet_grid2(vars(Filtered_Data_AK_Dec2020), vars(Group), scales = "free_x", independent = "x") +
    theme_bw() +
    theme(axis.text.y = element_text(size=10, face="bold"),axis.title.y = element_text(size=12,face="bold"),
           legend.title = element_text(size=12,face="bold"),
           legend.text=element_text(size=10,face="bold"),
           axis.text.x = element_blank(),axis.title.x=element_text(size=12,face="bold"),
           axis.ticks.x = element_blank())+
    theme(strip.text.x=element_text(size=13,face="bold"),strip.text.y=element_text(size=13,face="bold"))

    labs(x = "Sample", y = "Normalized Abundance", fill = "Genus")

ggsave(plot = p , filename =  "/Users/yangjichoi/Documents/PhD/Projects/3_PIPOVAP/Bar_lung_230123V2.tiff",width = 8, height = 8, dpi = 300)
    
############################################################################################################################

##Repeat the same procedure as before for the Gut samples. Create barplots of the 10 most abundant Genera
ps_tmp_gut<-subset_samples(physeqGut2,Subgroup=="Inclusion"|Subgroup=="Discharge")

View(psmelt(ps_tmp_gut))
View(psmelt(physeqGut2))

ps_tmp_gut<-microbiomeutilities::aggregate_top_taxa2(ps_tmp_gut, top=10, "Genus")

ps_tmp_gut <- name_taxa(ps_tmp_gut, label = "Unkown", species= T, other_label = "Other")

ps_tmp_gut<-ps_tmp_gut%>%microbiome::transform(transform="compositional")

tax_table(ps_tmp_gut)


############################################################################################################################
otu_df <- psmelt(ps_tmp_gut)

# Normalize the data by sample depths
# Replace `size` with your abundance column if it's named differently
otu_df$normalized_abundance <- otu_df$Abundance / sample_sums(ps_tmp_gut)[otu_df$Sample]

# Convert to long format for ggplot2
otu_df <- otu_df %>%
    group_by(Sample) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    ungroup()

# Create the stacked bar plot, faceting by both 'Group' and 'Inclusion/Extubation' status
g <- ggplot(otu_df, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values=mypal)+scale_color_manual(values=mypal) + # Adjust the palette as necessary
    # facet_grid(Filtered_Data_AK_Dec2020 ~ Group,  scales = "free") + # Assuming you have Inclusion_Extubation as a combined factor
    facet_grid2(vars(Subgroup), vars(Group), scales = "free_x", independent = "x") +
    theme_bw() +
    theme(axis.text.y = element_text(size=10, face="bold"),axis.title.y = element_text(size=12,face="bold"),
          legend.title = element_text(size=12,face="bold"),
          legend.text=element_text(size=10,face="bold"),
          axis.text.x = element_blank(),axis.title.x=element_text(size=12,face="bold"),
          axis.ticks.x = element_blank())+
    theme(strip.text.x=element_text(size=13,face="bold"),strip.text.y=element_text(size=13,face="bold")) + labs(x = "Sample", y = "Normalized Abundance", fill = "Genus")

ggsave(plot = g, filename  = "/Users/yangjichoi/Documents/PhD/Projects/3_PIPOVAP/Bar_gut_230123V2.tiff",width = 8, height = 8, dpi = 300)

############################################################################################################################

