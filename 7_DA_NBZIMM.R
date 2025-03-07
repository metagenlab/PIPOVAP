
## Load Libraries
```{r libraries}
library(rmarkdown)
library(knitr)
library(phyloseq)
library(vegan)
library(DT)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(pheatmap)
library(viridis)
library(gplots)
library(reshape)
library(readr)
library(microbiome)
library(ANCOMBC)
library(ggnewscale)
library(NBZIMM)
```

## Load and Process Phyloseq Objects
### Lung Phyloseq Object
```{r lung_phyloseq}
physeq_lung <- readRDS("lung10k.rds")
lungmeta <- read.csv("metadata_20210122_lung.csv")
lungmeta$timepoint <- lungmeta$InclusNoInfect_Infect1D1PossiblyInclus_Infect1D5_Extub_old
lungmeta <- lungmeta[c("SeqSampleID", "AKSampleID", "record_ID", "Group", "GutSampleAvailable", "Time_point", "Subgroup", "timepoint", "Subgroup4", "days_ICU", "days_atb")]

# Merge metadata
meta <- meta(physeq_lung)
metamerge <- merge(meta, lungmeta, by.x = "OldSampleName", by.y = "SeqSampleID", all.x = TRUE)
metamerge$GroupPneumonia <- ifelse(metamerge$Group == "Pneumonia", "Pneumonia", "non_Pneumonia")
rownames(metamerge) <- metamerge$sample_label
sample_data(physeq_lung) <- metamerge
metalung <- meta(physeq_lung)
```

### Gut Phyloseq Object
```{r gut_phyloseq}
physeq_gut <- readRDS("gut100k.rds")
physeq_gut <- subset_samples(physeq_gut, sample_type == 'rectal_swab')

gutmeta <- read.csv("metadata_20210122_gut.csv")
gutmeta$timepoint <- gutmeta$InclusNoInfect_Infect1D1PossiblyInclus_Infect1D5_Extub_old
gutmeta <- gutmeta[c("Sample", "ancien_smpl", "days_ICU", "days_atb")]

meta <- meta(physeq_gut)
metamerge <- merge(meta, gutmeta, by.x = "OldSampleName", by.y = "ancien_smpl", all.x = TRUE)
metamerge$GroupPneumonia <- ifelse(metamerge$Group == "Pneumonia", "Pneumonia", "non_Pneumonia")
rownames(metamerge) <- metamerge$Sample
sample_data(physeq_gut) <- metamerge
metagut <- meta(physeq_gut)
```

## Define Color Palettes
```{r palettes}
Infectcolorgroup <- c("Control" = "#8DD3C7", "Other_infection" = "#FDB462", "Pneumonia" = "#80B1D3")
InclInfectColors3 <- c("Inclusion" = "gold2", "Infection_D5" = "chartreuse4", "Discharge" = "cadetblue3")
InclInfectColors5 <- c("Inclusion" = "gold2", "Infection_D1" = "coral2", "Infection_D5" = "chartreuse4", "Extubation" = "cadetblue3", "Discharge" = "orchid1")
```

## Analyze Lung Phyloseq Object
```{r lung_analysis}
phyobject1 <- subset_samples(physeq_lung, timepoint == "Inclusion")
phyobject <- aggregate_taxa(phyobject1, level = "Genus")

myotu <- as.matrix(t(phyobject@otu_table))
mymeta <- meta(phyobject)
mymeta$days_ICU <- as.numeric(as.character(mymeta$days_ICU))
mymeta$days_atb <- as.numeric(as.character(mymeta$days_atb))

f_all_time <- mms(
    y = myotu,
    fixed = ~ GroupPneumonia + days_ICU,
    data = mymeta,
    random = ~ 1 | record_ID,
    zi_fixed = ~1,
    zi_random = NULL,
    method = "zinb",
    min.p = 0.3
)

out <- fixed(f_all_time)$dist
out <- out[out[, 2] != "(Intercept)", ]
res <- out[, 3:5]
Cova <- "GroupPneumoniaPneumonia"
res_cov <- res[grep(paste("--", Cova, sep = ""), rownames(res)), ]
rownames(res_cov) <- gsub(pattern = paste("--", Cova, sep = ""), replacement = "", x = rownames(res_cov))

res_cov_cutoff <- subset(res_cov, pvalue < 0.05 & abs(Estimate) > 0)
bactaxa_levels <- res_cov_cutoff

bactaxa_levels$p <- ""
bactaxa_levels$p[bactaxa_levels$pvalue <= 0.05 & bactaxa_levels$pvalue > 0.01] <- "*"
bactaxa_levels$p[bactaxa_levels$pvalue <= 0.01 & bactaxa_levels$pvalue > 0.001] <- "**"
bactaxa_levels$p[bactaxa_levels$pvalue <= 0.001] <- "***"

bactaxa_levels$Compare <- ifelse(bactaxa_levels$Estimate > 0, "Pneumonia", "Non-pneumonia")

piponb <- ggplot(bactaxa_levels, aes(x = Estimate, y = reorder(taxa, Estimate), fill = Compare)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("grey", "#80B1D3")) +
    theme_light() +
    theme(
        text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom"
    ) +
    ylab("") +
    xlab("Effect size") +
    ggtitle("Pneumonia vs Non-pneumonia groups at inclusion") +
    geom_text(aes(label = p, hjust = Estimate < 0), size = 3)

piponb
```
