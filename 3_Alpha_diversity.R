# title: "PIPOVAP_AK_05.05.2022 Alpha diversity composition over time among groups"
# date: "05.05.2022"
# output: ggplot object
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

## Import data

## Import Gut metadata file
meta_gut <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/3.Alpha diversity/metadata_gut.xlsx",sheet="working_sheet")
View(meta_gut)

## Import Lung metadata file
meta_lung <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/3.Alpha diversity/metadata_lung.xlsx",sheet = "Feuil1" )
View(meta_lung)

## Pre-processing gut metadata
colnames(meta_gut) # reval names of meta_gut columns

meta_gut[sapply(meta_gut, is.character)] <- lapply(meta_gut[sapply(meta_gut, is.character)], 
                                                   as.factor)   #transform character vectors to factors
class(meta_gut$Filtered_Data_AK_Dec2020)
class(meta_gut$Group)

levels(meta_gut$Filtered_Data_AK_Dec2020)
levels(meta_gut$Subgroup)

Correct_order_Incl_Inf_gut <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation", "Discharge")  
Correct_order_subgroup_gut<-c("Inclusion","Infection_D5","Discharge")

meta_gut$Filtered_Data_AK_Dec2020 <- factor(meta_gut$Filtered_Data_AK_Dec2020,
                                                            levels = Correct_order_Incl_Inf_gut)

meta_gut$Subgroup<-factor(meta_gut$Subgroup,levels=Correct_order_subgroup_gut)

levels(meta_gut$Filtered_Data_AK_Dec2020)
levels(meta_gut$Subgroup)


## Pre-processing lung metadata

colnames(meta_lung) # reveal names of meta_gut columns

meta_lung[sapply(meta_lung, is.character)] <- lapply(meta_lung[sapply(meta_lung, is.character)], 
                                                   as.factor)   #transform character vectors to factors

class(meta_lung$Filtered_Data_AK_Dec2020)
class(meta_lung$Group)

levels(meta_lung$Filtered_Data_AK_Dec2020)

Correct_order_Incl_Inf_lung <- c("Inclusion", "Infection_D1", "Infection_D5", "Extubation")  

meta_lung$Filtered_Data_AK_Dec2020 <- factor(meta_lung$Filtered_Data_AK_Dec2020,
                                            levels = Correct_order_Incl_Inf_lung)

levels(meta_lung$Filtered_Data_AK_Dec2020)

Correct_order_subgroup_lung<-c("Inclusion","Infection_D5","Discharge")


meta_lung$Filtered_Data_AK_Dec2020 <- factor(meta_lung$Filtered_Data_AK_Dec2020,
                                            levels = Correct_order_Incl_Inf_lung)

meta_lung$Subgroup<-factor(meta_lung$Subgroup,levels=Correct_order_subgroup_lung)

levels(meta_lung$Filtered_Data_AK_Dec2020)
levels(meta_lung$Subgroup)



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


## Plot Alpha diversity index as calculated directly by the output of the pipeline.

### Gut alpha-diversity

## Alpha diversity only on inclusion time_point

## Shannon
s<-meta_gut %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(label.x = 2, label.y = 1)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Shannon diversity index", 
       x = "Inclusion", y = "Shannon", fill = "Group")

ggsave("shannon_inclusion_gut.tiff",width = 25, height = 25, dpi = 300, units = "cm")

## same graph as before only showing labels instead of jitter_points
meta_gut %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(label.x = 2, label.y = 5)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Shannon diversity index", 
       x = "Inclusion", y = "Shannon", fill = "Group")+
  geom_text(position = position_jitter(seed=1),size=3,check_overlap = T)

ggsave("shannon_inclusion_gut_labels_insteadof_points.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Chao1
ts<-meta_gut %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(label.x = 2, label.y = 1)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Chao1 diversity index", 
       x = "Inclusion", y = "Chao1", fill = "Group")

ggsave("Chao1_inclusion_gut_.tiff",width = 25, height = 25, dpi = 300, units = "cm")

ggarrange(s,ts,nrow=1,common.legend=TRUE,legend = "right")
ggsave("gut_alpha_diver_inclusion.tiff",width = 8,height = 6,dpi=300)

getwd()
setwd("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/3.Alpha diversity")

# Control again alpha diversity on inclusion but test for pairwise statistical significance

my_comparisons <- list( c("Control", "Other_infection"), c("Other_infection", "Pneumonia"), c("Control", "Pneumonia"))

meta_gut %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 50,label.x=2)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Chao1 diversity index", 
       x = "Inclusion", y = "Chao1", fill = "Group")


ggsave("Chao1_inclusion_gut_pairwisecomparison.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Alpha-diversity in different timepoints:inclusion-->D1-->D5-->Extubation-->Discharge

## GUT alpha diversity measure ##

## Shannon
meta_gut %>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity over time",subtitle = "Shannon diversity index", 
       x = "Group", y = "Shannon", fill = "Timepoint")+
  scale_fill_manual(values = InclInfectColors5)

ggsave("shannon_overtime_gut.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Same graph as above but with sample names instead of jitter points
meta_gut %>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  theme_classic()+
  geom_jitter(position = position_jitter(seed = 1))+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity over time",subtitle = "Shannon diversity index", 
       x = "Group", y = "Shannon", fill = "Timepoint")+
  geom_text(position = position_jitter(seed=1),size=3,check_overlap = T)+
  scale_fill_manual(values = InclInfectColors5)

ggsave("shannon_overtime_gut_with_labels.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Perform individual plots for every group with pair-wise comparisons

over_time <- list( c("Inclusion", "Infection_D1"), c("Inclusion", "Infection_D5"), c("Inclusion", "Extubation"),c("Inclusion", "Discharge"))

over_time2 <- list( c("Inclusion", "Extubation"),c("Inclusion", "Discharge"))

c
meta_gut %>%subset(Group=="Control")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  scale_y_continuous(breaks=seq(from=1,to=6,by=1),limits =c(1,7))+
  theme_classic()+
  stat_compare_means(comparisons=over_time2)+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "",subtitle = "",
       x="Control",  y = "Shannon", fill = "")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)

c<-c+guides(fill=FALSE)

c
c<-c+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
      axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 
c

## Repeat for other infections and pneumonia patients

o<-meta_gut %>%subset(Group=="Other_infection")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
   theme_classic()+
  stat_compare_means(comparisons=over_time)+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "",subtitle = "",
       x="Other_infection",  y = "Shannon", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=1,to=6,by=1),limits =c(1,7))

o<-o+guides(fill=FALSE)
o
o<-o+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
        axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 

o<-o+theme(axis.title.y = element_blank())
o

p<-meta_gut %>%subset(Group=="Pneumonia")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time)+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "",subtitle = "",
       x="Pneumonia",  y = "Shannon", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=1,to=6,by=1),limits =c(1,7))
p

p<-p+guides(fill = FALSE) 
p
p<-p+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
        axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 
p

## Combine plots

shannon_combined<-ggarrange(c+guides(fill = FALSE),o+guides(fill = FALSE),p+guides(fill = FALSE),ncol=3,nrow=1,
          common.legend = T,legend="right")

shannon_combined

annotate_figure(shannon_combined,top=text_grob("Alpha diversity over time",face="bold",size=14))

ggsave("shannon_overtime_gut_with_pairwise_comparison.tiff",width = 10, height = 6, dpi = 300)


# Test alpha diversity by separating subgroups

plot<-meta_gut %>%
  ggplot(aes(x = Subgroup, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Subgroup),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity overtime",subtitle = "Stratification by Subgroup",
       x="",  y = "Shannon", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)


plot+facet_wrap(~Group)

ggsave("shannon_overtime_gut_subgroups.tiff",width = 25, height = 25, dpi = 300, units = "cm")


### Repeat the above analyses for the Gut for Chao1 index

cch<-meta_gut %>%subset(Group=="Control")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time2)+
  scale_fill_manual(values = InclInfectColors5)+
  labs(title = "",subtitle = "",
       x="Control",  y = "Chao1", fill = "timepint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=seq(from=100,to=600,by=100),limits =c(100,800))

cch
# CAVE: If we ad the following command it will produce comparison of all groups among them
# stat_compare_means(label.x = 2, label.y = 0.5)
  

cch<-cch+guides(fill=FALSE)

cch<-cch+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 

cch
## Repeat for other infections and pneumonia patients

och<-meta_gut %>%subset(Group=="Other_infection")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time)+
  labs(title = "",subtitle = "",
       x="Other_infection",  y = "Chao1", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=100,to=600,by=100),limits =c(100,800))

# stat_compare_means(label.x = 2, label.y = 105)

och

och<-och+guides(fill=FALSE)
och<-och+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
              axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
              axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black"))+
  theme(axis.title.y = element_blank())
och


pch<-meta_gut %>%subset(Group=="Pneumonia")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time)+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "",subtitle = "",
       x="Pneumonia",  y = "Chao1", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=100,to=600,by=100),limits =c(100,800))

pch<-pch+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black"))+
  theme(axis.title.y = element_blank())

pch

## Combine plots

chao1_combined<-ggarrange(cch,och,pch,ncol=3,nrow=1,
                            common.legend = T,legend="right")

chao1_combined

annotate_figure(chao1_combined,top=text_grob("Alpha diversity over time",face="bold",size=14))

ggsave("chao1_overtime_gut_with_pairwise_comparison.tiff",width = 10, height = 6, dpi = 300)


ggarrange(shannon_combined,chao1_combined,labels=c("C","D"),ncol=2,nrow=1,
          common.legend = T,legend="right")

ggsave("gut_diversity_overtime.tiff",width = 20, height = 6, dpi = 300)


# Test alpha diversity by separating subgroups

plot<-meta_gut %>%
  ggplot(aes(x = Subgroup, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Subgroup),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity overtime",subtitle = "Stratification by Subgroup",
       x="",  y = "Chao1", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)


plot+facet_wrap(~Group)

ggsave("chao1_overtime_gut_subgroups.tiff",width = 25, height = 25, dpi = 300, units = "cm")



#### Plot Lung alpha-diversity as calculated by the output of the pipeline

## LUNG alpha diversity measure ##

## Alpha diversity only on inclusion time_point

## Shannon
s<-meta_lung %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(label.x = 2, label.y = 1)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Shannon diversity index", 
       x = "Inclusion", y = "Shannon", fill = "Group")+
  scale_y_continuous(breaks=seq(from=1,to=6,by=1),limits =c(1,6))

ggsave("shannon_inclusion_lung.tiff",width = 25, height = 25, dpi = 300, units = "cm")

## same graph as before only showing labels instead of jitter_points
meta_lung %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(label.x = 2, label.y = 5)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Shannon diversity index", 
       x = "Inclusion", y = "Shannon", fill = "Group")+
  geom_text(position = position_jitter(seed=1),size=3,check_overlap = T)

ggsave("shannon_inclusion_lung_labels_insteadof_points.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Chao1
ts<-meta_lung %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(label.x = 2, label.y = 1)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Chao1 diversity index", 
       x = "Inclusion", y = "Chao1", fill = "Group")

ts

ggarrange(s,ts,nrow = 1,common.legend = TRUE,legend="right")
ggsave("lung_alpha_diver_incl.tiff",width=8,height = 6,dpi=300)

ggsave("Chao1_inclusion_lung_.tiff",width = 25, height = 25, dpi = 300, units = "cm")


# Control again alpha diversity on inclusion but test for pairwise statistical significance

my_comparisons <- list( c("Control", "Other_infection"), c("Other_infection", "Pneumonia"), c("Control", "Pneumonia"))

meta_lung %>%subset(Filtered_Data_AK_Dec2020=="Inclusion")%>%
  ggplot(aes(x = Group, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Group),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 50,label.x=2)+
  labs(title = "Alpha diversity on inclusion",subtitle = "Chao1 diversity index", 
       x = "Inclusion", y = "Chao1", fill = "Group")


ggsave("Chao1_inclusion_lung_pairwisecomparison.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Alpha-diversity in different timepoints:inclusion-->D1-->D5-->Extubation-->Discharge

## Shannon
meta_lung %>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity over time",subtitle = "Shannon diversity index", 
       x = "Group", y = "Shannon", fill = "Timepoint")+
  scale_fill_manual(values = InclInfectColors5)

ggsave("shannon_overtime_lung.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Same graph as above but with sample names instead of jitter points
meta_lung %>%
  ggplot(aes(x = Group, y = Shannon,label=Sample)) +
  geom_boxplot(aes(fill = Filtered_Data_AK_Dec2020),
               alpha = 0.4,
               position = position_dodge(0.9),
               outlier.colour = NULL,
               outlier.size = 1.5,
               outlier.stroke = 0,
               outlier.alpha = 0.5) +
  theme_classic()+
  geom_jitter(position = position_jitter(seed = 1))+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity over time",subtitle = "Shannon diversity index", 
       x = "Group", y = "Shannon", fill = "Timepoint")+
  geom_text(position = position_jitter(seed=1),size=3,check_overlap = T)+
  scale_fill_manual(values = InclInfectColors5)

ggsave("shannon_overtime_lung_with_labels.tiff",width = 25, height = 25, dpi = 300, units = "cm")


## Perform individual plots for every group with pair-wise comparisons

over_time <- list( c("Inclusion", "Infection_D1"), c("Inclusion", "Infection_D5"), c("Inclusion", "Extubation"))

over_time2 <- list( c("Inclusion", "Extubation"))


c<-meta_lung %>%subset(Group=="Control")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time2)+
  stat_compare_means(label.x = 1.5, label.y = 0.5)+
  labs(title = "",subtitle = "",
       x="Control",  y = "Shannon", fill = "")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=1,to=6,by=1),limits =c(1,6))

c

c<-c+guides(fill=FALSE)
c
c<-c+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 
c


## Repeat for other infections and pneumonia patients

o<-meta_lung %>%subset(Group=="Other_infection")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time)+
  stat_compare_means(label.x = 2, label.y = 0.3)+
  labs(title = "",subtitle = "",
       x="Other_infection",  y = "Shannon", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=1,to=6,by=1),limits =c(1,6))

o
o<-o+guides(fill=FALSE)

o<-o+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 
o
o<-o+theme(axis.title.y = element_blank())
o



p<-meta_lung %>%subset(Group=="Pneumonia")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time)+
  stat_compare_means(label.x = 2, label.y = 0.25)+
  labs(title = "",subtitle = "",
       x="Pneumonia",  y = "Shannon", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=1,to=6,by=1),limits =c(1,6))

p
p<-p+guides(fill=FALSE)

p<-p+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 

p<-p+theme(axis.title.y = element_blank())
p

## Combine plots

shannon_combined<-ggarrange(c,o,p,ncol=3,nrow=1,
                            common.legend = T,legend="right")
shannon_combined

annotate_figure(shannon_combined,top=text_grob("Alpha diversity over time",face="bold",size=14))

ggsave("shannon_overtime_lung_with_pairwise_comparison.tiff",width = 10, height = 6, dpi = 300)


# Test alpha diversity by separating subgroups

plot<-meta_lung %>%
  ggplot(aes(x = Subgroup, y = Shannon,label=Sample)) +
  geom_boxplot(
    aes(fill = Subgroup),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity overtime",subtitle = "Stratification by Subgroup",
       x="",  y = "Shannon", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)


plot+facet_wrap(~Group)

ggsave("shannon_overtime_lung_subgroups.tiff",width = 30, height = 25, dpi = 300, units = "cm")


### Repeat the above analyses for the lung for Chao1 index

cch<-meta_lung %>%subset(Group=="Control")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time2)+
  stat_compare_means(label.x = 1.5, label.y = 0.5)+
  scale_fill_manual(values = InclInfectColors5)+
  labs(title = "",subtitle = "",
       x="Control",  y = "Chao1", fill = "timepint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=seq(from=0,to=200,by=50),limits =c(1,300))

cch

cch<-cch+guides(fill=FALSE)
cch

cch<-cch+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 
cch


## Repeat for other infections and pneumonia patients

och<-meta_lung %>%subset(Group=="Other_infection")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time)+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "",subtitle = "",
       x="Other_infection",  y = "Chao1", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=0,to=200,by=50),limits =c(1,300))

och
och<-och+guides(fill=FALSE)

och<-och+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 
och
och<-och+theme(axis.title.y = element_blank())
och


pch<-meta_lung %>%subset(Group=="Pneumonia")%>%
  ggplot(aes(x = Filtered_Data_AK_Dec2020, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Filtered_Data_AK_Dec2020),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(comparisons=over_time)+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "",subtitle = "",
       x="Pneumonia",  y = "Chao1", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)+
  scale_y_continuous(breaks=seq(from=0,to=200,by=50),limits =c(1,300))

pch

pch<-pch+theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0),size=10, face="bold", family= "Times New Roman"),
           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12, face="bold", family= "Times New Roman"),
           axis.text.y = element_text(size=10, face="bold", family= "Times New Roman", colour = "black")) 
pch
pch<-pch+theme(axis.title.y = element_blank())
pch

## Combine plots

chao1_combined<-ggarrange(cch,och,pch,ncol=3,nrow=1,
                          common.legend = T,legend="right")

chao1_combined

annotate_figure(chao1_combined,top=text_grob("Alpha diversity over time",face="bold",size=14))


ggarrange(shannon_combined,chao1_combined,labels=c("A","B"),ncol=2,nrow=1,
          common.legend = T,legend="right")

ggsave("lung_alpha_diversity_overtime.tiff",width = 20, height = 6, dpi = 300)





# Test alpha diversity by separating subgroups

plot<-meta_lung %>%
  ggplot(aes(x = Subgroup, y = Chao1,label=Sample)) +
  geom_boxplot(
    aes(fill = Subgroup),
    alpha = 0.4,
    position = position_dodge(0.9),
    outlier.colour = NULL,
    outlier.size = 1.5,
    outlier.stroke = 0,
    outlier.alpha = 0.5) +
  geom_jitter(position = position_jitter(seed = 1))+
  theme_classic()+
  stat_compare_means(label.x = 2, label.y = 0.5)+
  labs(title = "Alpha diversity overtime",subtitle = "Stratification by Subgroup",
       x="",  y = "Chao1", fill = "Timepoint")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = InclInfectColors5)

plot
plot+facet_wrap(~Group)

ggsave("chao1_overtime_lung_subgroups.tiff",width = 35, height = 25, dpi = 300, units = "cm")



#########
#########
#########

