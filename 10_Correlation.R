install.packages("corrplot")
library(corrplot)
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
install.packages("Hmisc")
library(Hmisc)
install.packages("psych")
library(psych)
library (readxl)

getwd()
setwd("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/Correlation analysis")

correlation_matrix <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/Correlation analysis/correlation_matrix.xlsx",sheet="correlation_matrix")

View(correlation_matrix)

correlation_matrix<-correlation_matrix[,-(2:4)]
View(correlation_matrix)
correlation_matrix<-correlation_matrix[,-1]
View(correlation_matrix)

sapply(correlation_matrix,class)
correlation_matrix$"Lactates_infection D1"<-as.numeric(correlation_matrix$"Lactates_infection D1")

class(correlation_matrix)
colnames(correlation_matrix)

# cor(correlation_matrix, method = c("pearson"))

# mcor<-cor(correlation_matrix, use = "complete.obs")

# mcor<-cor(correlation_matrix, use = "pairwise.complete.obs",method="pearson")

cor_2<-rcorr(as.matrix(correlation_matrix))
mcor<-cor_2$r
p.mat<-cor_2$P

View(p.mat)

View(mcor)
class(mcor)
mcor2<-as.data.frame(mcor)
class(mcor2)  
View(mcor2)
str(mcor2)

mcor2<-mcor2[1:9,]
View(mcor2)
mcor2<-mcor2[,10:17]

class(mcor2)
mcor2<-as.matrix(mcor2)

View(p.mat)
class(p.mat)
p.mat<-as.data.frame(p.mat)
class(p.mat)  
View(p.mat)
str(p.mat)

p.mat<-p.mat[1:9,]
View(p.mat)
p.mat<-p.mat[,10:17]

class(p.mat)
p.mat<-as.matrix(p.mat)


corrplot(mcor2, method="circle",tl.col="black",tl.srt=45,tl.cex=1.2,
         addCoef.col="black")


# cor.mtest <- function(mat, ...) {
#  mat <- as.matrix(mat)
#  n <- ncol(mat)
#  p.mat<- matrix(NA, n, n)
#  diag(p.mat) <- 0
#  for (i in 1:(n - 1)) {
#    for (j in (i + 1):n) {
#      tmp <- cor.test(mat[, i], mat[, j], ...)
#      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#    }
#  }
#  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#  p.mat
# }

# Matrice de p-value de la corr?lation
# p.mat <- cor.mtest(correlation_matrix)

View(p.mat)

corrplot(mcor2,method="color",
         addCoef.col="black",
         tl.col="black",tl.srt=45,
         p.mat=p.mat,sig.level = 0.05,insig="blank")


tiff("corrplot_alpha_divers.tiff", units="in", width=5.7, height=6, res=300)

corrplot(mcor2,method="circle",
         tl.col="black",tl.srt=90,
         p.mat=p.mat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.5,
         insig = 'label_sig', pch.col = 'grey20')

dev.off()


tiff("corrplot_alpha_divers_stat.tiff", units="in", width=8, height=6, res=300)

corrplot(mcor2,method="circle",
         tl.col="black",tl.srt=90)

dev.off()


####################################################################################################

## Create a dataframe with relative abundance of taxa of interest to test correlation for

####################################################################################################

## First check for the Lung taxa

View(psmelt(physeqLung2))

count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqLungcor <- transform_sample_counts(physeqLung2, count_to_rel_abund)

View(psmelt(physeqLungcor))

phylo<-psmelt(physeqLungcor)
class(phylo)
View(phylo)

phylo_bergeyella<-filter(phylo,Genus=="Bergeyella"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_bergeyella)
phylo_bergeyella<-select(phylo_bergeyella,OTU,Sample,Abundance,Genus)
phylo_bergeyella<-rename(phylo_bergeyella,Bergeyella=Abundance)
phylo_bergeyella<-select(phylo_bergeyella,Sample,Bergeyella)
phylo_bergeyella<-aggregate(. ~ Sample, phylo_bergeyella, sum)
phylo_bergeyella<-rename(phylo_bergeyella,"Bergeyella"=Bergeyella)
View(phylo_bergeyella)

phylo_rumino<-filter(phylo,Genus=="PAC000661_g"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_rumino)
phylo_rumino<-select(phylo_rumino,Sample,Abundance)
phylo_rumino<-rename(phylo_rumino,"Ruminococcaceae_genus"=Abundance)
phylo_rumino<-aggregate(. ~ Sample, phylo_rumino, sum)
View(phylo_rumino)

phylo_rothia<-filter(phylo,Genus=="Rothia"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_rothia)
phylo_rothia<-select(phylo_rothia,Sample,Abundance)
phylo_rothia<-rename(phylo_rothia,"Rothia"=Abundance)
View(phylo_rothia)
phylo_rothia<-aggregate(. ~ Sample, phylo_rothia, sum)

phylo_coryne<-filter(phylo,Genus=="Corynebacterium"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_coryne)
phylo_coryne<-select(phylo_coryne,Sample,Abundance)
phylo_coryne<-rename(phylo_coryne,"Corynebacterium"=Abundance)
View(phylo_coryne)
phylo_coryne<-aggregate(. ~ Sample, phylo_coryne, sum)

phylo_treponema<-filter(phylo,Genus=="Treponema"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_treponema)
phylo_treponema<-select(phylo_treponema,Sample,Abundance)
phylo_treponema<-rename(phylo_treponema,"Treponema"=Abundance)
View(phylo_treponema)
phylo_treponema<-aggregate(. ~ Sample, phylo_treponema, sum)


phylo_actinobacteria<-filter(phylo,Phylum=="Actinobacteria"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_actinobacteria)
phylo_actinobacteria<-select(phylo_actinobacteria,Sample,Abundance)
phylo_actinobacteria<-rename(phylo_actinobacteria,"Actinobacteria"=Abundance)
phylo_actinobacteria<-aggregate(. ~ Sample, phylo_actinobacteria, sum)

phylo_firmicutes<-filter(phylo,Phylum=="Firmicutes"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_firmicutes)
phylo_firmicutes<-select(phylo_firmicutes,Sample,Abundance)
phylo_firmicutes<-rename(phylo_firmicutes,"Firmicutes"=Abundance)
View(phylo_firmicutes)
phylo_firmicutes<-aggregate(. ~ Sample, phylo_firmicutes, sum)

phylo_bacteroidetes<-filter(phylo,Phylum=="Bacteroidetes"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_bacteroidetes)
phylo_bacteroidetes<-select(phylo_bacteroidetes,Sample,Abundance)
phylo_bacteroidetes<-rename(phylo_bacteroidetes,"Bacteroidetes"=Abundance)
View(phylo_bacteroidetes)
phylo_bacteroidetes<-aggregate(. ~ Sample, phylo_bacteroidetes, sum)

phylo_spirochaetes<-filter(phylo,Phylum=="Spirochaetes"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_spirochaetes)
phylo_spirochaetes<-select(phylo_spirochaetes,Sample,Abundance)
phylo_spirochaetes<-rename(phylo_spirochaetes,"Spirochaetes"=Abundance)
View(phylo_spirochaetes)
phylo_spirochaetes<-aggregate(. ~ Sample, phylo_spirochaetes, sum)
phylo_spirochaetes<-rename(phylo_spirochaetes,"#SampleID"="Sample")


lung_taxa<-left_join(phylo_bergeyella,phylo_rumino,by="Sample")
lung_taxa<-left_join(lung_taxa,phylo_gemella,by="Sample")
lung_taxa<-left_join(lung_taxa,phylo_rothia,by="Sample")
lung_taxa<-left_join(lung_taxa,phylo_coryne,by="Sample")
lung_taxa<-left_join(lung_taxa,phylo_treponema,by="Sample")
lung_taxa<-left_join(lung_taxa,phylo_actinobacteria,by="Sample")
lung_taxa<-left_join(lung_taxa,phylo_firmicutes,by="Sample")
lung_taxa<-left_join(lung_taxa,phylo_spirochaetes,by="#SampleID")
View(lung_taxa)
lung_taxa<-left_join(lung_taxa,phylo_bacteroidetes,by="Sample")
View(lung_taxa)
lung_taxa<-rename(lung_taxa,"#SampleID"="Sample")

#Calculate Actinobacteria/Bacteroidetes ratio
lung_taxa$ratio=lung_taxa$Actinobacteria/lung_taxa$Bacteroidetes
lung_taxa<-rename(lung_taxa,"Actinobacteria/Bacteroidetes ratio (inclusion)"=ratio)
lung_taxa<-rename(lung_taxa,"Actinobacteria/Bacteroidetes ratio"="Actinobacteria/Bacteroidetes ratio (inclusion)")
View(lung_taxa)

#Calculate Firmicutes/Bacteroidetes ratio
lung_taxa$ratio=lung_taxa$Firmicutes/lung_taxa$Bacteroidetes
lung_taxa<-rename(lung_taxa,"Firmictutes/Bacteroidetes ratio"=ratio)
View(lung_taxa)

#Calculate Spirochaetes/Bacteroidetes ratio
lung_taxa$ratio=lung_taxa$Spirochaetes/lung_taxa$Bacteroidetes
lung_taxa<-rename(lung_taxa,"Spirochaetes/Bacteroidetes ratio"=ratio)
View(lung_taxa)

lung_taxa<-lung_taxa[,-12]


#Reload the correlation_matix because the #SampleID was previously deleted for the shake of previous analyses
correlation_matrix <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/Correlation analysis/correlation_matrix.xlsx",sheet="correlation_matrix")
View(correlation_matrix)

taxa_corr<-left_join(correlation_matrix,lung_taxa,by="#SampleID")
View(taxa_corr)


### Create correlation plot
taxa_corr<-taxa_corr[,-(2:4)]
View(taxa_corr)
taxa_corr<-taxa_corr[,-1]
View(taxa_corr)

# We can moreover delete from this dataframe data regarding alpha-diversity since the correlation plot was already produced previously about alpha-diversity.
View(taxa_corr)
taxa_corr<-taxa_corr[,-(10:17)]
View(taxa_corr)

sapply(taxa_corr,class)
taxa_corr$"Lactates_infection D1"<-as.numeric(taxa_corr$"Lactates_infection D1")

class(taxa_corr)

cor_2<-rcorr(as.matrix(taxa_corr))
cor_2
mcor<-cor_2$r
View(mcor)
p.mat<-cor_2$P


mcor2<-as.data.frame(mcor)
class(mcor2)  
View(mcor2)
str(mcor2)

# Visualize scatter plots
pairs(mcor)
chart.Correlation(mcor,histogram = TRUE,pch=19)

mcor2<-mcor2[1:9,]
View(mcor2)
mcor2<-mcor2[,10:20]
View(mcor2)

class(mcor2)
mcor2<-as.matrix(mcor2)

View(p.mat)
class(p.mat)
p.mat<-as.data.frame(p.mat)
class(p.mat)  
View(p.mat)
str(p.mat)

p.mat<-p.mat[1:9,]
View(p.mat)
p.mat<-p.mat[,10:20]

class(p.mat)
p.mat<-as.matrix(p.mat)
View(p.mat)

corrplot(mcor2, method="circle",tl.col="black",tl.srt=45,tl.cex=1.2,
         addCoef.col="black")

corrplot(mcor2,method="color",
         addCoef.col="black",
         tl.col="black",tl.srt=45,
         p.mat=p.mat,sig.level = 0.05,insig="blank")


# Define the pathway to save the image of the correlation plot

tiff("taxacorr_lung_inclusion.tiff", units="in", width=8, height=6, res=300)

corrplot(mcor2,method="circle",
         tl.col="black",tl.srt=90,
         p.mat=p.mat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.5,
         insig = 'label_sig', pch.col = 'grey20')

dev.off()


tiff("test2.tiff", units="in", width=8, height=6, res=300)

corrplot(mcor2,method="circle",
         tl.col="black",tl.srt=45)

dev.off()


## Now check for the gut microbiome taxa

View(psmelt(physeqGut2))

count_to_rel_abund <- function(x) {return( x / sum(x))}

physeqGutcor <- transform_sample_counts(physeqGut2, count_to_rel_abund)

View(psmelt(physeqGutcor))

phylo<-psmelt(physeqGutcor)
class(phylo)
View(phylo)

phylo_actinotignum<-filter(phylo,Genus=="Actinotignum"&Filtered_Data_AK_Dec2020=="Inclusion")
phylo_actinotignum<-select(phylo_actinotignum,Sample,Abundance)
phylo_actinotignum<-rename(phylo_actinotignum,"Actinotignum"=Abundance)
phylo_actinotignum<-aggregate(. ~ Sample, phylo_actinotignum, sum)
View(phylo_actinotignum)

phylo_christen<-filter(phylo,Genus=="PAC001360_g"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_christen)
phylo_christen<-select(phylo_christen,Sample,Abundance)
phylo_christen<-rename(phylo_christen,"Christensenellaceae_genus"=Abundance)
phylo_christen<-aggregate(. ~ Sample, phylo_christen, sum)
View(phylo_christen)

phylo_eubacterium_g21<-filter(phylo,Genus=="Eubacterium_g21"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_eubacterium_g21)
phylo_eubacterium_g21<-select(phylo_eubacterium_g21,Sample,Abundance)
phylo_eubacterium_g21<-rename(phylo_eubacterium_g21,"Eubacterium_g21"=Abundance)
phylo_eubacterium_g21<-aggregate(. ~ Sample, phylo_eubacterium_g21, sum)
View(phylo_eubacterium_g21)

phylo_lachno<-filter(phylo,Genus=="PAC001046_g"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_lachno)
phylo_lachno<-select(phylo_lachno,Sample,Abundance)
phylo_lachno<-rename(phylo_lachno,"Lachnospiraceae_genus"=Abundance)
phylo_lachno<-aggregate(. ~ Sample, phylo_lachno, sum)
View(phylo_lachno)

phylo_lachno2<-filter(phylo,Genus=="LLKB_g"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_lachno2)
phylo_lachno2<-select(phylo_lachno2,Sample,Abundance)
phylo_lachno2<-rename(phylo_lachno2,"Lachnospiraceae_genus(2)"=Abundance)
phylo_lachno2<-aggregate(. ~ Sample, phylo_lachno2, sum)
View(phylo_lachno2)

phylo_butyricimonas<-filter(phylo,Genus=="Butyricimonas"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_butyricimonas)
phylo_butyricimonas<-select(phylo_butyricimonas,Sample,Abundance)
phylo_butyricimonas<-rename(phylo_butyricimonas,"Butyricimonas"=Abundance)
phylo_butyricimonas<-aggregate(. ~ Sample, phylo_butyricimonas, sum)
View(phylo_butyricimonas)

phylo_atopobium<-filter(phylo,Genus=="Atopobium"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_atopobium)
phylo_atopobium<-select(phylo_atopobium,Sample,Abundance)
phylo_atopobium<-rename(phylo_atopobium,"Atopobium"=Abundance)
phylo_atopobium<-aggregate(. ~ Sample, phylo_atopobium, sum)
View(phylo_atopobium)

phylo_lachno3<-filter(phylo,Genus=="PAC001177_g"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_lachno3)
phylo_lachno3<-select(phylo_lachno3,Sample,Abundance)
phylo_lachno3<-rename(phylo_lachno3,"Lachnospiraceae_genus(3)"=Abundance)
phylo_lachno3<-aggregate(. ~ Sample, phylo_lachno3, sum)
View(phylo_lachno3)

phylo_coprobacter<-filter(phylo,Genus=="Coprobacter"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_coprobacter)
phylo_coprobacter<-select(phylo_coprobacter,Sample,Abundance)
phylo_coprobacter<-rename(phylo_coprobacter,"Coprobacter"=Abundance)
phylo_coprobacter<-aggregate(. ~ Sample, phylo_coprobacter, sum)
View(phylo_coprobacter)

phylo_caproiciproducens<-filter(phylo,Genus=="Caproiciproducens"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_caproiciproducens)
phylo_caproiciproducens<-select(phylo_caproiciproducens,Sample,Abundance)
phylo_caproiciproducens<-rename(phylo_caproiciproducens,"Caproiciproducens"=Abundance)
phylo_caproiciproducens<-aggregate(. ~ Sample, phylo_caproiciproducens, sum)
View(phylo_caproiciproducens)

phylo_mogibacterium<-filter(phylo,Genus=="Mogibacterium"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_mogibacterium)
phylo_mogibacterium<-select(phylo_mogibacterium,Sample,Abundance)
phylo_mogibacterium<-rename(phylo_mogibacterium,"Mogibacterium"=Abundance)
phylo_mogibacterium<-aggregate(. ~ Sample, phylo_mogibacterium, sum)
View(phylo_mogibacterium)

phylo_lachno4<-filter(phylo,Genus=="Eubacterium_g5"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_lachno4)
phylo_lachno4<-select(phylo_lachno4,Sample,Abundance)
phylo_lachno4<-rename(phylo_lachno4,"Eubacterium_g5"=Abundance)
phylo_lachno4<-aggregate(. ~ Sample, phylo_lachno4, sum)
View(phylo_lachno4)

phylo_ruminoco<-filter(phylo,Genus=="PAC001100_g"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_ruminoco)
phylo_ruminoco<-select(phylo_ruminoco,Sample,Abundance)
phylo_ruminoco<-rename(phylo_ruminoco,"Ruminococcaceae_genus"=Abundance)
phylo_ruminoco<-aggregate(. ~ Sample, phylo_ruminoco, sum)
View(phylo_ruminoco)


phylo_actinobacteria_gut<-filter(phylo,Phylum=="Actinobacteria"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_actinobacteria_gut)
phylo_actinobacteria_gut<-select(phylo_actinobacteria_gut,Sample,Abundance)
phylo_actinobacteria_gut<-rename(phylo_actinobacteria_gut,"Actinobacteria"=Abundance)
phylo_actinobacteria_gut<-aggregate(. ~ Sample, phylo_actinobacteria_gut, sum)

phylo_firmicutes_gut<-filter(phylo,Phylum=="Firmicutes"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_firmicutes_gut)
phylo_firmicutes_gut<-select(phylo_firmicutes_gut,Sample,Abundance)
phylo_firmicutes_gut<-rename(phylo_firmicutes_gut,"Firmicutes"=Abundance)
View(phylo_firmicutes_gut)
phylo_firmicutes_gut<-aggregate(. ~ Sample, phylo_firmicutes_gut, sum)

phylo_bacteroidetes_gut<-filter(phylo,Phylum=="Bacteroidetes"&Filtered_Data_AK_Dec2020=="Inclusion")
View(phylo_bacteroidetes_gut)
phylo_bacteroidetes_gut<-select(phylo_bacteroidetes_gut,Sample,Abundance)
phylo_bacteroidetes_gut<-rename(phylo_bacteroidetes_gut,"Bacteroidetes"=Abundance)
View(phylo_bacteroidetes_gut)
phylo_bacteroidetes_gut<-aggregate(. ~ Sample, phylo_bacteroidetes_gut, sum)



## Recapitulatif des phyla que nous allons tester
## : phylo_actinotignum; phylo_christen; phylo_eubacterium_g21;phylo_lachno;
# phylo_lachno2;phylo_butyricimonas; phylo_atopobium;
# phylo_lachno3; phylo_coprobacter; phylo_caproiciproducens
# phylo_mogibacterium; phylo_lachno4; phylo_ruminoco



gut_taxa<-left_join(phylo_actinotignum,phylo_christen,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_eubacterium_g21,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_lachno,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_lachno2,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_butyricimonas,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_atopobium,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_lachno3,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_coprobacter,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_caproiciproducens,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_mogibacterium,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_lachno4,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_ruminoco,by="Sample")

View(gut_taxa)

gut_taxa<-left_join(gut_taxa,phylo_actinobacteria_gut,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_firmicutes_gut,by="Sample")
gut_taxa<-left_join(gut_taxa,phylo_bacteroidetes_gut,by="Sample")

View(gut_taxa)
gut_taxa<-rename(gut_taxa,"#SampleID"="Sample")

#Calculate Firmicutes/Bacteroidetes ratio
gut_taxa$ratio=gut_taxa$Firmicutes/gut_taxa$Bacteroidetes
gut_taxa<-rename(gut_taxa,"Firmictutes/Bacteroidetes ratio"=ratio)
View(gut_taxa)

#Calculate Firmicutes/Actinobacteria ratio
gut_taxa$ratio=gut_taxa$Firmicutes/gut_taxa$Actinobacteria
gut_taxa<-rename(gut_taxa,"Firmictutes/Actinobacteria ratio"=ratio)
View(gut_taxa)


#Reload the correlation_matix because the #SampleID was p reviously deleted for the shake of previous analyses
correlation_matrix <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/Correlation analysis/correlation_matrix.xlsx",sheet="correlation_matrix")
View(correlation_matrix)

taxa_corr<-left_join(correlation_matrix,gut_taxa,by="#SampleID")
View(taxa_corr)


### Create correlation plot
taxa_corr<-taxa_corr[,-(2:4)]
View(taxa_corr)
taxa_corr<-taxa_corr[,-1]
View(taxa_corr)

# We can moreover delete from this dataframe data regarding alpha-diversity since the correlation plot was already produced previously about alpha-diversity.
View(taxa_corr)
taxa_corr<-taxa_corr[,-(10:17)]
View(taxa_corr)

sapply(taxa_corr,class)
taxa_corr$"Lactates_infection D1"<-as.numeric(taxa_corr$"Lactates_infection D1")

class(taxa_corr)

cor_2<-rcorr(as.matrix(taxa_corr))
mcor<-cor_2$r
p.mat<-cor_2$P

mcor2<-as.data.frame(mcor)
class(mcor2)  
View(mcor2)
str(mcor2)

# Visualize scatter plots
pairs(mcor)
chart.Correlation(mcor,histogram = TRUE,pch=19)

mcor2<-mcor2[1:9,]
View(mcor2)
mcor2<-mcor2[,10:27]
View(mcor2)

class(mcor2)
mcor2<-as.matrix(mcor2)

View(p.mat)
class(p.mat)
p.mat<-as.data.frame(p.mat)
class(p.mat)  
View(p.mat)
str(p.mat)

p.mat<-p.mat[1:9,]
View(p.mat)
p.mat<-p.mat[,10:27]

class(p.mat)
p.mat<-as.matrix(p.mat)
View(p.mat)

corrplot(mcor2, method="circle",tl.col="black",tl.srt=45,tl.cex=1.2,
         addCoef.col="black")

corrplot(mcor2,method="color",
         addCoef.col="black",
         tl.col="black",tl.srt=45,
         p.mat=p.mat,sig.level = 0.05,insig="blank")


# Define the pathway to save the image of the correlation plot

tiff("taxacorr_gut_inclusion.tiff", units="in", width=13.2, height=6, res=300)

cor_gut<-corrplot(mcor2,method="circle",
         tl.col="black",tl.srt=90,
         p.mat=p.mat, sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.5,
         insig = 'label_sig', pch.col = 'grey20')

dev.off()


tiff("test2.tiff", units="in", width=8, height=6, res=300)

corrplot(mcor2,method="circle",
         tl.col="black",tl.srt=45)

dev.off()
