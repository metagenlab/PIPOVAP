# title: "PIPOVAP"
# date: "15.02.23"
# Multivariate ordination plot
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

install.packages("xts")
library(xts)
BiocManager::install("RCM")
library("RCM")
browseVignettes("RCM")
library(ggplot2)
library(microViz)
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

# So for the multivariate model we will use the "pg2" phyloseq object which contains only the "discharge" timepoint of the Gut phyloseq object
View(meta(pg2))
meta<-(meta(pg2))
View(meta)

View(metadata_for_rda)
names(metadata_for_rda)
metaRCM<-select(metadata_for_rda,Sample,Sex,Age,BMI,days_ICU,tot_num_atb,days_atb,type_of_infect,atb_yes_no)
View(metaRCM)

meta<-left_join(meta,metaRCM,by="Sample")
View(meta)
names(meta)
meta<-select(meta,Sample,Count,record_ID,Group,Subgroup,Time_point,Filtered_Data_AK_Dec2020,Shannon,Chao1,Age,Sex,BMI,days_ICU,days_atb,tot_num_atb,atb_yes_no)
View(meta)

sample_names(pg2)

sam <- sample_data(meta)
View(sam)


View(sample_data(pg2))
sample_data(pg2)<-sam
View(psmelt(pg2))

# For consistency with the other elements of phyloseq object, use "sample_label"
# as sample_names

sample_names(sam) <- as.factor(sam$Sample)
head(sample_names(sam))

class(sam$Sample)
sample_names(sam)

# Create new phyloseq object with "pg2" phyloseq object containing only discharge samples from Gut

physeqRCM<-pg2 
head(sample_names(physeqRCM))
sample_variables(physeqRCM)

View(psmelt(physeqRCM))
     
    
# Rename our phyloseq object for the shake of subsequent analyses and to facilitate next analyses
Zeller=phyRCM
class(Zeller)
View(psmelt(Zeller))

##########################################
##So now we will use the Zeller phyloseq object to create our multivariate ordination plot
##########################################

View(psmelt(Zeller))

class(Zeller)

Zeller <- tax_fix(Zeller)
Zeller <- phyloseq_validate(Zeller, remove_undetected = TRUE)
View(psmelt(Zeller))

Zeller %>% 
  tax_transform(trans = "clr", rank = "Genus")


Zeller %>% 
  tax_transform(trans = "identity", rank = "Genus") %>% 
  dist_calc("bray")


Zeller %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Group",  size = 3) +
  scale_colour_manual(values = Infectcolorgroup)

Zeller %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Group", plot_taxa = 1:5, size = 2.5) +
  scale_colour_manual(values = Infectcolorgroup)


Zeller %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot_iris(tax_level = "Genus", ord_plot = "above", anno_colour = "Group")

Zeller %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("jaccard",binary=TRUE) %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "Group", size = 2) +
  scale_colour_brewer(palette = "Dark2")

Zeller %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "Group", size = 2) +
  scale_colour_brewer(palette = "Dark2")

Zeller %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("bray") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "Group", size = 2) +
  scale_colour_brewer(palette = "Dark2")


Zeller %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("bray") %>% 
  ord_calc("PCoA") %>% 
  ord_get() %>% 
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 6))


ibd %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("bray") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(axes = c(1, 3), color = "ibd", shape = "DiseaseState", size = 2) +
  scale_colour_brewer(palette = "Dark2") 


Zeller %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("bray") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(axes = c(1, 3), color = "Group", size = 2) +
  scale_colour_brewer(palette = "Dark2") 



Zeller %>% 
  tax_transform("identity", rank = "Genus") %>% 
  dist_calc(dist = "aitchison") %>% 
  ord_calc("NMDS") %>% 
  ord_plot(color = "Group", size = 2) +
  scale_colour_brewer(palette = "Dark2", aesthetics = c("fill", "colour")) +
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = Group, y = Group), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Group, x = Group), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()



Zeller %>%
  ps_mutate(
    ATB = as.numeric(atb_yes_no == "1"),
    Female = as.numeric(Sex == "2"),
    Days=as.numeric(days_ICU),
    BMI=as.numeric(BMI)
  ) %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c("ATB", "Female","Days","BMI"),
    # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "Group", size = 2, alpha = 0.5, plot_taxa = 1:8
  )


renamer <- function(x) str_replace(x, pattern = "_", replacement = " ")

View(sam)

Zeller %>%
  ps_mutate(
    ATB = as.numeric(atb_yes_no == "1"),
    Female = as.numeric(Sex == "2"),
    Days=as.numeric(days_ICU),
    BMI=as.numeric(BMI)
  ) %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c("ATB", "Female","Days","BMI"),
     # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "Group", size = 3, alpha = 0.5,
    auto_caption = NA, # remove the helpful automatic caption
    plot_taxa = 0.5, taxon_renamer = renamer, # renamer is the function we made earlier
    tax_vec_length = 5, # this value is a scalar multiplier for the biplot score vectors
    tax_lab_length = 5, # this multiplier moves the labels, independently of the arrowheads
    tax_lab_style = tax_lab_style(size = 2, alpha = 0.5,type="text"), # create a list of options to tweak the taxa labels' default styleconstraint_vec_length = 3, # this adjusts the length of the constraint arrows, and the labels track these lengths by default
    constraint_vec_style = vec_constraint(size = 2, alpha = 0.5),# this styles the constraint arrows  
    constraint_lab_style = constraint_lab_style(size = 3) # this styles the constraint labels
  ) +
  # the functions below are from ggplot2:
  # You can pick a different colour scale, such as a color_brewer palette
  scale_colour_brewer(palette = "Set1")+
  # You can set any scale's values manually, such as the shapes used
  # this is how you add a title and subtitle
  ggtitle(
    label = "[Insert your exciting interpretations here?]",
    subtitle = "RDA with clr-transformed genera: constraints in red, taxa in black"
  )+
  # and this is how you make your own caption
  labs(caption = "23 samples, 392 genera. Type 2 scaling.") +
  # this is one way to set the aspect ratio of the plot
  coord_fixed(ratio = 1, clip = "off")


## Create numeric variables in our phyloseq object
ps <- Zeller %>%
  ps_mutate(
    female = if_else(Sex == "1", true = 1, false = 0),
    extract_C=if_else(Group=="Control",true=1,false=0),
    extract_P=if_else(Group=="Pneumonia",true=1,false=0),
    extract_O=if_else(Group=="Other_infection",true=1,false=0),
    days=as.numeric(days_ICU)
  ) 

View(psmelt(ps))


# PERMANOVA when testing with Jaccard
perm<-ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>%
  dist_calc("jaccard",binary=TRUE) %>%
  dist_permanova(
    variables = c("days","tot_num_atb","days_atb","BMI","Sex","Group"),
    n_perms = 999, # this is a low number of permutations for speed, you should set more e.g. 9999
    seed = 12345, complete_cases = TRUE, verbose = "max"
  )


perm_get(perm)

# PERMANOVA when testing with Bray
perm<-ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>%
  dist_calc("bray") %>%
  dist_permanova(
    variables = c("days","tot_num_atb","days_atb","BMI","Sex","Group"),
    n_perms = 999, # this is a low number of permutations for speed, you should set more e.g. 9999
    seed = 12345, complete_cases = TRUE, verbose = "max"
  )


perm_get(perm)



# constraints need to be on the same or similar scales for comparability
# so make binary variables and scale the weight variable
ps <- ps %>%
  ps_mutate(
    Days_in_ICU = c(scale(days, center = TRUE, scale = TRUE)),
    Days_of_ATB = c(scale(days_atb, center = TRUE, scale = TRUE)),
    BMI = c(scale(BMI, center = TRUE, scale = TRUE)),
    Female_sex = if_else(Sex == "2", 1, 0))

View(ps_melt(ps))


constr_ord<-ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>%
  tax_transform("clr") %>%
  ord_calc(method = "RDA", 
           constraints = c("Female_sex", "Days_in_ICU", "Days_of_ATB", "BMI"))

ord_plot(constr_ord,
         plot_taxa= 1:1,  
         colour = "Group",size=3,
         tax_vec_length = 4,tax_lab_length =4,tax_lab_style =tax_lab_style(size = 0, alpha = 0.4,type="text"),
         constraint_vec_length = 5.5,
         constraint_lab_length = 5.8,
         constraint_lab_style =constraint_lab_style(size = 3.5,max_angle = 90,perpendicular = TRUE,fontface = "bold.italic") ) +
  scale_color_brewer(palette = "Set2", name = "Infection Group")+
  stat_ellipse(aes(linetype=Group,colour=Group),size=0.5)+
  coord_fixed(ratio = 0.5, clip = "off")+
  theme_bw() +
  ggside::geom_xsideboxplot(aes(fill = Group, y = Group), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Group, x = Group), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()+
  scale_colour_manual(values=Infectcolorgroup, name = "Infection Group")+
  scale_fill_manual(values=Infectcolorgroup)  


ggsave("RDA analysis_Jaccard_Gut_Discharge.tiff", width =16,height = 10)



