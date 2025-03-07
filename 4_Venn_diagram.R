# title: "PIPOVAP_AK_04/05/22"
# date: "04.05.22"
# Venn diagrams for shared species among different timepoints for the 3 groups of patients
# Output: Venn diagram figures 
---

# Install packages
  
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
install.packages("ggplotify")
library("ggplotify")

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


## Keep only taxa with total abundance at least 2 in at least 1 samples
keep_function <- function(x) {x >= 1}
TaxaToKeep <- genefilter_sample(physeqLung2, keep_function, A = 1)

TaxaToKeep
TaxaToKeep[1:20]

physeqLung3 <- prune_taxa(TaxaToKeep, physeqLung2)

sort(taxa_sums(physeqLung3), decreasing = F)[1:50]


View(sample_data(physeqLung2))
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
metadata_gut <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/220131_pipovap_gut/metadata_gut.xlsx")
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

levels(get_variable(physeqLung2,"Subgroup"))


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


##Produce Venn diagrams

### Venn diagram for Lung microbiota ###

## Number 1 :Patients with No infection AND (Inclusion + Discharge) available ----

# Compare Inclusion with discharge
PatNoInf_InclExtub <- c("7","20","24","39")

physNoInf_InclExtub <- subset_samples(physeqLung3, record_ID %in% PatNoInf_InclExtub)

# Keep only represented taxa
physNoInf_InclExtubPresent <- prune_taxa(taxa_sums(physNoInf_InclExtub) > 0, 
                                         physNoInf_InclExtub)

TotalTaxaPresent <- ntaxa(physNoInf_InclExtubPresent)
TotalTaxaPresent

# Subset taxa present either at Inclusion
physNoInfIncl <- subset_samples(physNoInf_InclExtubPresent,
                                Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physNoInfIncl))

# Keep only represented taxa
physNoInfInclPresent <- prune_taxa(taxa_sums(physNoInfIncl) > 0, 
                                   physNoInfIncl)

physNoInfInclPresent_df <- data.frame(tax_table(physNoInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physNoInfInclPresent_df)

Venn_area1 <- dim(physNoInfInclPresent_df)[1] # for area 1 in venn diagram
head(physNoInfInclPresent_df)
Venn_area1

# Subset taxa present at Extubation 
physNoInfExtub <- subset_samples(physNoInf_InclExtubPresent, Filtered_Data_AK_Dec2020 %in% 
                                   "Extubation")

# Keep only represented taxa
physNoInfExtubPresent <- prune_taxa(taxa_sums(physNoInfExtub) > 0, 
                                    physNoInfExtub)

physNoInfExtubPresent_df <- data.frame(tax_table(physNoInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physNoInfExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physNoInfInclPresent_df, 
                                 physNoInfExtubPresent_df, by = "ASV") %>%arrange(ASV)

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
Venn_cross_area

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
c<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                scaled=TRUE,
                                fill = c("gold2", "cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                cat.cex = 2,        # Numbers
                                cex = 2,print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_test.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)
c<-as.ggplot(grid.arrange(gTree(children = c), # Add title & subtitle
            bottom = textGrob("Patients (n)=4, ASVs (n)=436",gp=gpar(fontsize=18,font=7))))
c

ggsave("c.tiff",width = 4,height = 4,dpi=300)

#####################################################################################################################


## Number 2 :Patients with Pneumonia AND (Inclusion + Discharge) available ----


# Compare Inclusion with discharge
Patpneu_InclExtub <- c("28")

physpneu_InclExtub <- subset_samples(physeqLung3, record_ID %in% Patpneu_InclExtub)

# Keep only represented taxa
physpneu_InclExtubPresent <- prune_taxa(taxa_sums(physpneu_InclExtub) > 0, 
                                         physpneu_InclExtub)

TotalTaxaPresent <- ntaxa(physpneu_InclExtubPresent)
TotalTaxaPresent

# Subset taxa present either at Inclusion or (Inclusion AND Extubation)
physpneuIncl <- subset_samples(physpneu_InclExtubPresent,
                                Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                   physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Extubation or (Extubation AND Inclusion)
physpneuExtub <- subset_samples(physpneu_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                   "Extubation")

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                    physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                             physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physpneuInclPresent_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p1<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                   cross.area = Venn_cross_area,
                                   scaled=TRUE,
                                   fill = c("gold2","cadetblue3"),
                                   alpha = 0.25,
                                   lty = "blank",
                                   cat.cex = 2,
                                   # Numbers
                                   cex = 2,print.mode = "percent", sigdigs = 2,
                                   fontfamily = "sans",
                                   filename = "venn_test.pdf",
                                   output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)
p1<-as.ggplot(grid.arrange(gTree(children = p1), # Add title & subtitle
             bottom = textGrob("Patients (n)=1, ASVs (n)=165",gp=gpar(fontsize=18,font=7))))
             
p1
ggsave("p1.tiff",width=4,height=4,dpi=300)
 


#####################################################################################################################



## Number 3 :Patients with Pneumonia AND Inclusion + D1 available ----


# Compare Inclusion with discharge
Patpneu_InclD1 <- c("1","9","28","33","40","41")

physpneu_InclD1 <- subset_samples(physeqLung3, record_ID %in% Patpneu_InclD1)

View(sample_data(physeqLung))

# Keep only represented taxa
physpneu_InclD1Present <- prune_taxa(taxa_sums(physpneu_InclD1) > 0, 
                                        physpneu_InclD1)

TotalTaxaPresent <- ntaxa(physpneu_InclD1Present)
TotalTaxaPresent

# Subset taxa present either at Inclusion 
physpneuIncl <- subset_samples(physpneu_InclD1Present,
                               Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                  physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Infection_D1
physpneuD1<-subset_samples(physpneu_InclD1Present, 
                             Filtered_Data_AK_Dec2020 %in% "Infection_D1")

View(sample_data(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                   physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD1Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                             physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_D1 <- intersect(physpneuInclPresent_df$ASV, physpneuD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_D1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p2<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                   cross.area = Venn_cross_area,
                                   scaled=TRUE,
                                   fill = c("gold2", "coral2"),
                                   alpha = 0.25,
                                   lty = "blank",
                                   cat.cex = 2,
                                   # Numbers
                                   cex = 2,print.mode = "percent", sigdigs = 2,
                                   fontfamily = "sans",
                                   filename = "venn_test.pdf",
                                   output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

p2<-as.ggplot(grid.arrange(gTree(children = p2), # Add title & subtitle
                           bottom = textGrob("Patients (n)=6, ASVs (n)=689",gp=gpar(fontsize=18,font=7))))


ggsave("p2.tiff",width=4,height = 4,dpi=300)


#####################################################################################################################




## Number 4 :Patients with Pneumonia AND D1+D5 available ----


# Compare D1 vs D5
Patpneu_D1D5 <- c("2","8","28","30","36","41","43")

physpneu_D1D5 <- subset_samples(physeqLung3, record_ID %in% Patpneu_D1D5)

View(sample_data(physeqLung))

# Keep only represented taxa
physpneu_D1D5Present <- prune_taxa(taxa_sums(physpneu_D1D5) > 0, 
                                     physpneu_D1D5)

TotalTaxaPresent <- ntaxa(physpneu_D1D5Present)
TotalTaxaPresent

# Subset taxa present at D1 
physpneuD1 <- subset_samples(physpneu_D1D5Present,
                               Filtered_Data_AK_Dec2020 %in% "Infection_D1")


View(psmelt(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                  physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD1Present_df)

Venn_area1 <- dim(physpneuD1Present_df)[1] # for area 1 in venn diagram
head(physpneuD1Present_df)
Venn_area1

# Subset taxa present either at Infection_D5
physpneuD5<-subset_samples(physpneu_D1D5Present, 
                           Filtered_Data_AK_Dec2020 %in% "Infection_D5")

View(sample_data(physpneuD5))

# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD5Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                                 physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD1Only_df)

taxD1Only_df[1:10,]

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxD1_AND_D5 <- intersect(physpneuD1Present_df$ASV, physpneuD5Present_df$ASV)
Venn_cross_area <- length(taxD1_AND_D5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p3<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("coral2", "chartreuse4"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

p3<-as.ggplot(grid.arrange(gTree(children = p3), # Add title & subtitle
                           bottom = textGrob("Patients (n)=7, ASVs (n)=658",gp=gpar(fontsize=18,font=7))))



ggsave("p3.tiff",width = 4,height = 4,dpi=300)


#####################################################################################################################

## Number 5 :Patients with Pneumonia AND D5+Extub available ----


# Compare D5 vs Extub
Patpneu_D5Extub <- c("2","8","28","30")

physpneu_D5Extub <- subset_samples(physeqLung3, record_ID %in% Patpneu_D5Extub)

View(sample_data(physpneu_D5Extub))

# Keep only represented taxa
physpneu_D5ExtubPresent <- prune_taxa(taxa_sums(physpneu_D5Extub) > 0, 
                                   physpneu_D5Extub)

TotalTaxaPresent <- ntaxa(physpneu_D5ExtubPresent)
TotalTaxaPresent

# Subset taxa present at D5 
physpneuD5 <- subset_samples(physpneu_D5ExtubPresent,
                             Filtered_Data_AK_Dec2020 %in% "Infection_D5")


View(psmelt(physpneuD5))
View(sample_data(physpneuD5))


# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD5Present_df)

Venn_area1 <- dim(physpneuD5Present_df)[1] # for area 1 in venn diagram
head(physpneuD5Present_df)
Venn_area1

# Subset taxa present at Extubation
physpneuExtub<-subset_samples(physpneu_D5ExtubPresent, 
                           Filtered_Data_AK_Dec2020 %in% "Extubation")

View(sample_data(physpneuExtub))

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD5Only_df)

taxD5Only_df[1:10,]

# Get taxa names present at Extubatiob ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                          physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxD5Only_df)

# Taxa present at D5 and Extubation
taxD5_AND_Extub <- intersect(physpneuD5Present_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p4<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                cross.area = Venn_cross_area,
                                scaled=TRUE,
                                inverted = TRUE,
                                fill = c("chartreuse4","cadetblue3"),
                                alpha = 0.25,
                                lty = "blank",
                                cat.cex = 2,
                                # Numbers
                                cex = 2,print.mode = "percent", sigdigs = 2,
                                fontfamily = "sans",
                                filename = "venn_test.pdf",
                                output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

p4<-as.ggplot(grid.arrange(gTree(children = p4), # Add title & subtitle
                           bottom = textGrob("Patients (n)=4, ASVs (n)=485",gp=gpar(fontsize=18,font=7))))

ggsave("p4.tiff",width = 4,height = 4,dpi=300)



####################################################################################################################################################
####################################################################################################################################################
##    CAVE : We are going to repeat the above mentioned analyses only this time we will change the patient group and included patients.
##    Since all other parameters will be the same, for the shake of simplification we will keep the same script without changing the code.
##    This means that the scrupt keeps in the environment only the last run of code and that we will have to reproduce all the script if we want to 
##    print previous figures or data!!!
####################################################################################################################################################
####################################################################################################################################################


## Number 6 :Patients with other infection AND (Inclusion + Discharge) available ----


# Compare Inclusion with discharge
Patpneu_InclExtub <- c("3","18")

physpneu_InclExtub <- subset_samples(physeqLung3, record_ID %in% Patpneu_InclExtub)

# Keep only represented taxa
physpneu_InclExtubPresent <- prune_taxa(taxa_sums(physpneu_InclExtub) > 0, 
                                        physpneu_InclExtub)

TotalTaxaPresent <- ntaxa(physpneu_InclExtubPresent)
TotalTaxaPresent

# Subset taxa present either at Inclusion or (Inclusion AND Extubation)
physpneuIncl <- subset_samples(physpneu_InclExtubPresent,
                               Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                  physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Extubation or (Extubation AND Inclusion)
physpneuExtub <- subset_samples(physpneu_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                  "Extubation")

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                   physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                             physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physpneuInclPresent_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o1<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("gold2","cadetblue3"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)
o1<-as.ggplot(grid.arrange(gTree(children = o1), # Add title & subtitle
                           bottom = textGrob("Patients (n)=2, ASVs (n)=328",gp=gpar(fontsize=18,font=7))))

o1
ggsave("o1.tiff",width=4,height=4,dpi=300)



#####################################################################################################################



## Number 7 :Patients with Other infection AND Inclusion + D1 available ----


# Compare Inclusion with discharge
Patpneu_InclD1 <- c("3","16","18","32")

physpneu_InclD1 <- subset_samples(physeqLung3, record_ID %in% Patpneu_InclD1)

View(sample_data(physeqLung))

# Keep only represented taxa
physpneu_InclD1Present <- prune_taxa(taxa_sums(physpneu_InclD1) > 0, 
                                     physpneu_InclD1)

TotalTaxaPresent <- ntaxa(physpneu_InclD1Present)
TotalTaxaPresent

# Subset taxa present either at Inclusion 
physpneuIncl <- subset_samples(physpneu_InclD1Present,
                               Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                  physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Infection_D1
physpneuD1<-subset_samples(physpneu_InclD1Present, 
                           Filtered_Data_AK_Dec2020 %in% "Infection_D1")

View(sample_data(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD1Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                          physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_D1 <- intersect(physpneuInclPresent_df$ASV, physpneuD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_D1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o2<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("gold2", "coral2"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

o2<-as.ggplot(grid.arrange(gTree(children = o2), # Add title & subtitle
                           bottom = textGrob("Patients (n)=4, ASVs (n)=551",gp=gpar(fontsize=18,font=7))))


ggsave("o2.tiff",width=4,height = 4,dpi=300)


#####################################################################################################################




## Number 8 :Patients with Other infection AND D1+D5 available ----


# Compare D1 vs D5
Patpneu_D1D5 <- c("3","6","16","23","29","32","34")

physpneu_D1D5 <- subset_samples(physeqLung3, record_ID %in% Patpneu_D1D5)

View(sample_data(physeqLung))

# Keep only represented taxa
physpneu_D1D5Present <- prune_taxa(taxa_sums(physpneu_D1D5) > 0, 
                                   physpneu_D1D5)

TotalTaxaPresent <- ntaxa(physpneu_D1D5Present)
TotalTaxaPresent

# Subset taxa present at D1 
physpneuD1 <- subset_samples(physpneu_D1D5Present,
                             Filtered_Data_AK_Dec2020 %in% "Infection_D1")


View(psmelt(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD1Present_df)

Venn_area1 <- dim(physpneuD1Present_df)[1] # for area 1 in venn diagram
head(physpneuD1Present_df)
Venn_area1

# Subset taxa present either at Infection_D5
physpneuD5<-subset_samples(physpneu_D1D5Present, 
                           Filtered_Data_AK_Dec2020 %in% "Infection_D5")

View(sample_data(physpneuD5))

# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD5Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                          physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD1Only_df)

taxD1Only_df[1:10,]

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxD1_AND_D5 <- intersect(physpneuD1Present_df$ASV, physpneuD5Present_df$ASV)
Venn_cross_area <- length(taxD1_AND_D5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o3<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("coral2", "chartreuse4"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

o3<-as.ggplot(grid.arrange(gTree(children = o3), # Add title & subtitle
                           bottom = textGrob("Patients (n)=7, ASVs (n)=780",gp=gpar(fontsize=18,font=7))))


ggsave("o3.tiff",width = 4,height = 4,dpi=300)


#####################################################################################################################

## Number 9 :Patients with Other infection AND D5+Extub available ----


# Compare D5 vs Extub
Patpneu_D5Extub <- c("3","25")

physpneu_D5Extub <- subset_samples(physeqLung3, record_ID %in% Patpneu_D5Extub)

View(sample_data(physpneu_D5Extub))

# Keep only represented taxa
physpneu_D5ExtubPresent <- prune_taxa(taxa_sums(physpneu_D5Extub) > 0, 
                                      physpneu_D5Extub)

TotalTaxaPresent <- ntaxa(physpneu_D5ExtubPresent)
TotalTaxaPresent

# Subset taxa present at D5 
physpneuD5 <- subset_samples(physpneu_D5ExtubPresent,
                             Filtered_Data_AK_Dec2020 %in% "Infection_D5")


View(psmelt(physpneuD5))
View(sample_data(physpneuD5))


# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD5Present_df)

Venn_area1 <- dim(physpneuD5Present_df)[1] # for area 1 in venn diagram
head(physpneuD5Present_df)
Venn_area1

# Subset taxa present at Extubation
physpneuExtub<-subset_samples(physpneu_D5ExtubPresent, 
                              Filtered_Data_AK_Dec2020 %in% "Extubation")

View(sample_data(physpneuExtub))

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                   physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD5Only_df)

taxD5Only_df[1:10,]

# Get taxa names present at Extubatiob ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                             physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxD5Only_df)

# Taxa present at D5 and Extubation
taxD5_AND_Extub <- intersect(physpneuD5Present_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o4<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    inverted = TRUE,
                                    fill = c("chartreuse4","cadetblue3"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

o4<-as.ggplot(grid.arrange(gTree(children = o4), # Add title & subtitle
                           bottom = textGrob("Patients (n)=2, ASVs (n)=296",gp=gpar(fontsize=18,font=7))))

ggsave("o4.tiff",width = 4,height = 4,dpi=300)





####################################################################################################################################################
####################################################################################################################################################
##    CAVE : As previously we are going to repeat all the above analyses for the Gut microbiota.
##    In order to save time we will only change the corresponding Phyloseq object to wotk on the Gut and the patients of interest. 
##    Nevertheless and in order to save time we will keep the script and the abbreviations used meaning that even if we analyse the gut microbiota the 
##    factor codes we look like as if we treat the lung microbiota. So be careful at this point!!
####################################################################################################################################################
####################################################################################################################################################


##Produce Venn diagrams

### Venn diagram for Gut microbiota ###

## Number 1 :Patients with No infection AND Inclusion + Extubation available (we take extubation and not discharge in order to be comparable with lung samples)

# Compare Inclusion with discharge
PatNoInf_InclExtub <- c("7","19","20","24","39")

physNoInf_InclExtub <- subset_samples(physeqGut3, record_ID %in% PatNoInf_InclExtub)

# Keep only represented taxa
physNoInf_InclExtubPresent <- prune_taxa(taxa_sums(physNoInf_InclExtub) > 0, 
                                         physNoInf_InclExtub)

TotalTaxaPresent <- ntaxa(physNoInf_InclExtubPresent)
TotalTaxaPresent

# Subset taxa present either at Inclusion
physNoInfIncl <- subset_samples(physNoInf_InclExtubPresent,
                                Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physNoInfIncl))

# Keep only represented taxa
physNoInfInclPresent <- prune_taxa(taxa_sums(physNoInfIncl) > 0, 
                                   physNoInfIncl)

physNoInfInclPresent_df <- data.frame(tax_table(physNoInfInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physNoInfInclPresent_df)

Venn_area1 <- dim(physNoInfInclPresent_df)[1] # for area 1 in venn diagram
head(physNoInfInclPresent_df)
Venn_area1

# Subset taxa present at Extubation 
physNoInfExtub <- subset_samples(physNoInf_InclExtubPresent, Filtered_Data_AK_Dec2020 %in% 
                                   "Extubation")

# Keep only represented taxa
physNoInfExtubPresent <- prune_taxa(taxa_sums(physNoInfExtub) > 0, 
                                    physNoInfExtub)

physNoInfExtubPresent_df <- data.frame(tax_table(physNoInfExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physNoInfExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physNoInfInclPresent_df, 
                                 physNoInfExtubPresent_df, by = "ASV") %>%arrange(ASV)

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
Venn_cross_area

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
c<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                   cross.area = Venn_cross_area,
                                   scaled=TRUE,
                                   fill = c("gold2", "cadetblue3"),
                                   alpha = 0.25,
                                   lty = "blank",
                                   cat.cex = 2,        # Numbers
                                   cex = 2,print.mode = "percent", sigdigs = 2,
                                   fontfamily = "sans",
                                   filename = "venn_test.pdf",
                                   output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)
c<-as.ggplot(grid.arrange(gTree(children = c), # Add title & subtitle
                          bottom = textGrob("Patients (n)=5, ASVs (n)=1227",gp=gpar(fontsize=18,font=7))))
c

ggsave("c.tiff",width = 4,height = 4,dpi=300)

#####################################################################################################################


## Number 2 :Patients with Pneumonia AND (Inclusion + Extubation) available ----


# Compare Inclusion with discharge
Patpneu_InclExtub <- c("28")

physpneu_InclExtub <- subset_samples(physeqGut3, record_ID %in% Patpneu_InclExtub)

# Keep only represented taxa
physpneu_InclExtubPresent <- prune_taxa(taxa_sums(physpneu_InclExtub) > 0, 
                                        physpneu_InclExtub)

TotalTaxaPresent <- ntaxa(physpneu_InclExtubPresent)
TotalTaxaPresent

# Subset taxa present either at Inclusion or (Inclusion AND Extubation)
physpneuIncl <- subset_samples(physpneu_InclExtubPresent,
                               Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                  physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Extubation or (Extubation AND Inclusion)
physpneuExtub <- subset_samples(physpneu_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                  "Extubation")

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                   physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                             physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physpneuInclPresent_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p1<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("gold2","cadetblue3"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)
p1<-as.ggplot(grid.arrange(gTree(children = p1), # Add title & subtitle
                           bottom = textGrob("Patients (n)=1, ASVs (n)=737",gp=gpar(fontsize=18,font=7))))

p1
ggsave("p1.tiff",width=4,height=4,dpi=300)



#####################################################################################################################



## Number 3 :Patients with Pneumonia AND Inclusion + D1 available ----


# Compare Inclusion with discharge
Patpneu_InclD1 <- c("1","9","28","33","40","41")

physpneu_InclD1 <- subset_samples(physeqGut3, record_ID %in% Patpneu_InclD1)

View(sample_data(physeqLung))

# Keep only represented taxa
physpneu_InclD1Present <- prune_taxa(taxa_sums(physpneu_InclD1) > 0, 
                                     physpneu_InclD1)

TotalTaxaPresent <- ntaxa(physpneu_InclD1Present)
TotalTaxaPresent

# Subset taxa present either at Inclusion 
physpneuIncl <- subset_samples(physpneu_InclD1Present,
                               Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                  physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Infection_D1
physpneuD1<-subset_samples(physpneu_InclD1Present, 
                           Filtered_Data_AK_Dec2020 %in% "Infection_D1")

View(sample_data(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD1Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                          physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_D1 <- intersect(physpneuInclPresent_df$ASV, physpneuD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_D1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p2<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("gold2", "coral2"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

p2<-as.ggplot(grid.arrange(gTree(children = p2), # Add title & subtitle
                           bottom = textGrob("Patients (n)=6, ASVs (n)=1737",gp=gpar(fontsize=18,font=7))))


ggsave("p2.tiff",width=4,height = 4,dpi=300)


#####################################################################################################################




## Number 4 :Patients with Pneumonia AND D1+D5 available ----


# Compare D1 vs D5
Patpneu_D1D5 <- c("1","2","4","8","9","26","28","30","36","37","41","43")

physpneu_D1D5 <- subset_samples(physeqGut3, record_ID %in% Patpneu_D1D5)

View(sample_data(physpneu_D1D5))

# Keep only represented taxa
physpneu_D1D5Present <- prune_taxa(taxa_sums(physpneu_D1D5) > 0, 
                                   physpneu_D1D5)

TotalTaxaPresent <- ntaxa(physpneu_D1D5Present)
TotalTaxaPresent

# Subset taxa present at D1 
physpneuD1 <- subset_samples(physpneu_D1D5Present,
                             Filtered_Data_AK_Dec2020 %in% "Infection_D1")


View(psmelt(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD1Present_df)

Venn_area1 <- dim(physpneuD1Present_df)[1] # for area 1 in venn diagram
head(physpneuD1Present_df)
Venn_area1

# Subset taxa present either at Infection_D5
physpneuD5<-subset_samples(physpneu_D1D5Present, 
                           Filtered_Data_AK_Dec2020 %in% "Infection_D5")

View(sample_data(physpneuD5))

# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD5Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                          physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD1Only_df)

taxD1Only_df[1:10,]

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxD1_AND_D5 <- intersect(physpneuD1Present_df$ASV, physpneuD5Present_df$ASV)
Venn_cross_area <- length(taxD1_AND_D5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p3<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("coral2", "chartreuse4"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

p3<-as.ggplot(grid.arrange(gTree(children = p3), # Add title & subtitle
                           bottom = textGrob("Patients (n)=12, ASVs (n)=2413",gp=gpar(fontsize=18,font=7))))



ggsave("p3.tiff",width = 4,height = 4,dpi=300)


#####################################################################################################################

## Number 5 :Patients with Pneumonia AND D5+Extub available ----


# Compare D5 vs Extub
Patpneu_D5Extub <- c("28","30","37")

physpneu_D5Extub <- subset_samples(physeqGut3, record_ID %in% Patpneu_D5Extub)

View(sample_data(physpneu_D5Extub))

# Keep only represented taxa
physpneu_D5ExtubPresent <- prune_taxa(taxa_sums(physpneu_D5Extub) > 0, 
                                      physpneu_D5Extub)

TotalTaxaPresent <- ntaxa(physpneu_D5ExtubPresent)
TotalTaxaPresent

# Subset taxa present at D5 
physpneuD5 <- subset_samples(physpneu_D5ExtubPresent,
                             Filtered_Data_AK_Dec2020 %in% "Infection_D5")


View(psmelt(physpneuD5))
View(sample_data(physpneuD5))


# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD5Present_df)

Venn_area1 <- dim(physpneuD5Present_df)[1] # for area 1 in venn diagram
head(physpneuD5Present_df)
Venn_area1

# Subset taxa present at Extubation
physpneuExtub<-subset_samples(physpneu_D5ExtubPresent, 
                              Filtered_Data_AK_Dec2020 %in% "Extubation")

View(sample_data(physpneuExtub))

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                   physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD5Only_df)

taxD5Only_df[1:10,]

# Get taxa names present at Extubatiob ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                             physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxD5Only_df)

# Taxa present at D5 and Extubation
taxD5_AND_Extub <- intersect(physpneuD5Present_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
p4<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    inverted = TRUE,
                                    fill = c("chartreuse4","cadetblue3"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

p4<-as.ggplot(grid.arrange(gTree(children = p4), # Add title & subtitle
                           bottom = textGrob("Patients (n)=3, ASVs (n)=1229",gp=gpar(fontsize=18,font=7))))

ggsave("p4.tiff",width = 4,height = 4,dpi=300)



####################################################################################################################################################
####################################################################################################################################################
##    CAVE : We are going to repeat the above mentioned analyses only this time we will change the patient group and included patients.
##    Since all other parameters will be the same, for the shake of simplification we will keep the same script without changing the code.
##    This means that the scrupt keeps in the environment only the last run of code and that we will have to reproduce all the script if we want to 
##    print previous figures or data!!!
####################################################################################################################################################
####################################################################################################################################################


## Number 6 :Patients with other infection AND (Inclusion + Discharge) available ----


# Compare Inclusion with discharge
Patpneu_InclExtub <- c("16")

physpneu_InclExtub <- subset_samples(physeqGut3, record_ID %in% Patpneu_InclExtub)

# Keep only represented taxa
physpneu_InclExtubPresent <- prune_taxa(taxa_sums(physpneu_InclExtub) > 0, 
                                        physpneu_InclExtub)

TotalTaxaPresent <- ntaxa(physpneu_InclExtubPresent)
TotalTaxaPresent

# Subset taxa present either at Inclusion or (Inclusion AND Extubation)
physpneuIncl <- subset_samples(physpneu_InclExtubPresent,
                               Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                  physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Extubation or (Extubation AND Inclusion)
physpneuExtub <- subset_samples(physpneu_InclExtub, Filtered_Data_AK_Dec2020 %in% 
                                  "Extubation")

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                   physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at Extubation ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                             physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxExtubOnly_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_Extub <- intersect(physpneuInclPresent_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxIncl_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o1<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    inverted = TRUE,
                                    fill = c("gold2","cadetblue3"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)
o1<-as.ggplot(grid.arrange(gTree(children = o1), # Add title & subtitle
                           bottom = textGrob("Patients (n)=1, ASVs (n)=536",gp=gpar(fontsize=18,font=7))))

o1
ggsave("o1.tiff",width=4,height=4,dpi=300)



#####################################################################################################################



## Number 7 :Patients with Other infection AND Inclusion + D1 available ----


# Compare Inclusion with discharge
Patpneu_InclD1 <- c("3","16","32")

physpneu_InclD1 <- subset_samples(physeqGut3, record_ID %in% Patpneu_InclD1)

View(sample_data(physeqLung))

# Keep only represented taxa
physpneu_InclD1Present <- prune_taxa(taxa_sums(physpneu_InclD1) > 0, 
                                     physpneu_InclD1)

TotalTaxaPresent <- ntaxa(physpneu_InclD1Present)
TotalTaxaPresent

# Subset taxa present either at Inclusion 
physpneuIncl <- subset_samples(physpneu_InclD1Present,
                               Filtered_Data_AK_Dec2020 %in% "Inclusion")


View(psmelt(physpneuIncl))

# Keep only represented taxa
physpneuInclPresent <- prune_taxa(taxa_sums(physpneuIncl) > 0, 
                                  physpneuIncl)

physpneuInclPresent_df <- data.frame(tax_table(physpneuInclPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuInclPresent_df)

Venn_area1 <- dim(physpneuInclPresent_df)[1] # for area 1 in venn diagram
head(physpneuInclPresent_df)
Venn_area1

# Subset taxa present either at Infection_D1
physpneuD1<-subset_samples(physpneu_InclD1Present, 
                           Filtered_Data_AK_Dec2020 %in% "Infection_D1")

View(sample_data(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD1Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at Inclusion ONLY
taxInclusionOnly_df <- anti_join(physpneuInclPresent_df, 
                                 physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxInclusionOnly_df)

taxInclusionOnly_df[1:10,]

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                          physpneuInclPresent_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxIncl_AND_D1 <- intersect(physpneuInclPresent_df$ASV, physpneuD1Present_df$ASV)
Venn_cross_area <- length(taxIncl_AND_D1) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o2<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    inverted=TRUE,
                                    fill = c("gold2", "coral2"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

o2<-as.ggplot(grid.arrange(gTree(children = o2), # Add title & subtitle
                           bottom = textGrob("Patients (n)=3, ASVs (n)=1259",gp=gpar(fontsize=18,font=7))))


ggsave("o2.tiff",width=4,height = 4,dpi=300)


#####################################################################################################################




## Number 8 :Patients with Other infection AND D1+D5 available ----


# Compare D1 vs D5
Patpneu_D1D5 <- c("3","13","16","18","23","27","29","32","34","38")

physpneu_D1D5 <- subset_samples(physeqGut3, record_ID %in% Patpneu_D1D5)

View(sample_data(physeqLung))

# Keep only represented taxa
physpneu_D1D5Present <- prune_taxa(taxa_sums(physpneu_D1D5) > 0, 
                                   physpneu_D1D5)

TotalTaxaPresent <- ntaxa(physpneu_D1D5Present)
TotalTaxaPresent

# Subset taxa present at D1 
physpneuD1 <- subset_samples(physpneu_D1D5Present,
                             Filtered_Data_AK_Dec2020 %in% "Infection_D1")


View(psmelt(physpneuD1))

# Keep only represented taxa
physpneuD1Present <- prune_taxa(taxa_sums(physpneuD1) > 0, 
                                physpneuD1)

physpneuD1Present_df <- data.frame(tax_table(physpneuD1Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD1Present_df)

Venn_area1 <- dim(physpneuD1Present_df)[1] # for area 1 in venn diagram
head(physpneuD1Present_df)
Venn_area1

# Subset taxa present either at Infection_D5
physpneuD5<-subset_samples(physpneu_D1D5Present, 
                           Filtered_Data_AK_Dec2020 %in% "Infection_D5")

View(sample_data(physpneuD5))

# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuD5Present_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D1 ONLY
taxD1Only_df <- anti_join(physpneuD1Present_df, 
                          physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD1Only_df)

taxD1Only_df[1:10,]

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuD1Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxD5Only_df)
head(taxD5Only_df)

# Taxa present at Inclusion and Extubation
taxD1_AND_D5 <- intersect(physpneuD1Present_df$ASV, physpneuD5Present_df$ASV)
Venn_cross_area <- length(taxD1_AND_D5) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o3<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    inverted = TRUE,
                                    fill = c("coral2", "chartreuse4"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

o3<-as.ggplot(grid.arrange(gTree(children = o3), # Add title & subtitle
                           bottom = textGrob("Patients (n)=10, ASVs (n)=2464",gp=gpar(fontsize=18,font=7))))


ggsave("o3.tiff",width = 4,height = 4,dpi=300)


#####################################################################################################################

## Number 9 :Patients with Other infection AND D5+Extub available ----


# Compare D5 vs Extub
Patpneu_D5Extub <- c("13","16","18","25")

physpneu_D5Extub <- subset_samples(physeqGut3, record_ID %in% Patpneu_D5Extub)

View(sample_data(physpneu_D5Extub))

# Keep only represented taxa
physpneu_D5ExtubPresent <- prune_taxa(taxa_sums(physpneu_D5Extub) > 0, 
                                      physpneu_D5Extub)

TotalTaxaPresent <- ntaxa(physpneu_D5ExtubPresent)
TotalTaxaPresent

# Subset taxa present at D5 
physpneuD5 <- subset_samples(physpneu_D5ExtubPresent,
                             Filtered_Data_AK_Dec2020 %in% "Infection_D5")


View(psmelt(physpneuD5))
View(sample_data(physpneuD5))


# Keep only represented taxa
physpneuD5Present <- prune_taxa(taxa_sums(physpneuD5) > 0, 
                                physpneuD5)

physpneuD5Present_df <- data.frame(tax_table(physpneuD5Present)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)

View(physpneuD5Present_df)

Venn_area1 <- dim(physpneuD5Present_df)[1] # for area 1 in venn diagram
head(physpneuD5Present_df)
Venn_area1

# Subset taxa present at Extubation
physpneuExtub<-subset_samples(physpneu_D5ExtubPresent, 
                              Filtered_Data_AK_Dec2020 %in% "Extubation")

View(sample_data(physpneuExtub))

# Keep only represented taxa
physpneuExtubPresent <- prune_taxa(taxa_sums(physpneuExtub) > 0, 
                                   physpneuExtub)

physpneuExtubPresent_df <- data.frame(tax_table(physpneuExtubPresent)) %>% 
  rownames_to_column("ASV") %>% 
  arrange(ASV)
Venn_area2 <- dim(physpneuExtubPresent_df)[1] # for area 2 in venn diagram
Venn_area2

# Get taxa names present at D5 ONLY
taxD5Only_df <- anti_join(physpneuD5Present_df, 
                          physpneuExtubPresent_df, by = "ASV") %>% 
  arrange(ASV)

dim(taxD5Only_df)

taxD5Only_df[1:10,]

# Get taxa names present at Extubatiob ONLY
taxExtubOnly_df <- anti_join(physpneuExtubPresent_df, 
                             physpneuD5Present_df, by = "ASV") %>% 
  arrange(ASV)
dim(taxExtubOnly_df)
head(taxD5Only_df)

# Taxa present at D5 and Extubation
taxD5_AND_Extub <- intersect(physpneuD5Present_df$ASV, physpneuExtubPresent_df$ASV)
Venn_cross_area <- length(taxD5_AND_Extub) # for calculation of cross.area in venn diagram

# For venn diagram
Venn_area1
Venn_area2
Venn_cross_area
TotalTaxaPresent

# Move to new plotting page
grid.newpage()

# Create Venn diagram (set print.mode to "raw" or "percent")
o4<-VennDiagram::draw.pairwise.venn(area1 = Venn_area1, area2 = Venn_area2,
                                    cross.area = Venn_cross_area,
                                    scaled=TRUE,
                                    fill = c("chartreuse4","cadetblue3"),
                                    alpha = 0.25,
                                    lty = "blank",
                                    cat.cex = 2,
                                    # Numbers
                                    cex = 2,print.mode = "percent", sigdigs = 2,
                                    fontfamily = "sans",
                                    filename = "venn_test.pdf",
                                    output=TRUE) # Export > Save as PDF > 4 x 5 inches


require(gridExtra)

o4<-as.ggplot(grid.arrange(gTree(children = o4), # Add title & subtitle
                           bottom = textGrob("Patients (n)=4, ASVs (n)=1103",gp=gpar(fontsize=18,font=7))))

ggsave("o4.tiff",width = 4,height = 4,dpi=300)





######   END OF CURRENT SCRIPT HERE    ########


######   END OF CURRENT SCRIPT HERE    ########


######   END OF CURRENT SCRIPT HERE    ########


