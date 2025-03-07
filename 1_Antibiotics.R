
remotes::install_github("giocomai/ganttrify")

library("ganttrify")
install.packages("vistime")
library("vistime")
library("readxl")
library("openxlsx")
library("ggplot2")
library("tidyverse")
library("tidyselect")
library("tidyr")


Classeur2<-read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/Demographics and Data Filtering/wp_for_antibiotics_pneumonia.xlsx")
View(Classeur2)

spots <- read_excel("//file3.intranet.chuv/data3/SHARE/PIPOVAP/PIPOVAP script _171121_AK/Demographics and Data Filtering/spots_for_antibiotics_pneumonia.xlsx")
View(spots)


atb_pneumonia<-ganttrify(Classeur2,
          spots,
          project_start_date = "2021-01",
          font_family = "Roboto Condensed",exact_date = FALSE,
          month_date_label=FALSE,
          by_date = FALSE,
          axis_text_align = "left")+
  labs(title = "Antibiotic treatment received by patients with pneumonia",
       caption = "I=inclusion, D1=1st day of infection, D5=5th day of infection, I2=second infectious episode, Ex=extubation, R=release from ICU",
       fontface="bold")


atb_pneumonia

ggsave("atb_pneumonia.tiff",
width = 35, height = 25, dpi = 300, units = "cm")




##########################################################
#2nd part: test for another kind of graph for atb followup

tggplot(data = followup_atb, aes(x = record_id, y = days, fill = received_atb)) + 
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_brewer(palette="Set3")


followup_atb$group<-as.factor(followup_atb$group)
followup_atb$record_id<-as.factor(followup_atb$record_id)

followup_atb$received_atb<-as.factor(followup_atb$received_atb)

followup_atb$received_atb<-revalue(followup_atb$received_atb,c("1"="amoxicillin","5"="cefepime","12"="amikacin","13"="gentamicin","2"="amoxicillin/clavulanic acid","3"="ceftriaxone","9"="levofloxacin","11"="piperacillin/tazobactam","16"="vancomycin","20"="meropenem","23"="clarithromycin","25"="ertapenem"))

ggplot(filter(followup_atb,followup_atb$group=="pneumonia"), aes(x = record_id, y = days, fill = received_atb)) + 
  geom_bar(stat="identity")+
  scale_fill_brewer(colScale)+
  geom_text(aes(label = days), size = 3, hjust = 0.1, vjust = 0.3, position = "stack")+
  coord_flip()+
  labs(title = "Antibiotic treatment received by pneumonia patients", 
       x = "Patient identification number", y = "Days of treatment", fill = "Administered antibiotic")

table1::label(followup_atb$days) <- "Days of treatment"
table1::label(followup_atb$record_id) <- "Patient"

atb

atb_pneumonia<-atb+facet_wrap(~ infectious_episode, ncol=2)


atb2<-ggplot(filter(followup_atb,followup_atb$group=="other_infection"), aes(x = record_id, y = days, fill = received_atb)) + 
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Set3")+
  geom_text(aes(label = days), size = 3, hjust = 0.1, vjust = 0.3, position = "stack")+
  coord_flip()+
  labs(title = "Antibiotic treatment received by patients with other infections", 
       x = "Patient identification number", y = "Days of treatment", fill = "Administered antibiotic")


atb2

atb_other<-atb2+facet_wrap(~ infectious_episode, ncol=2)

atb_other

ggarrange(atb_pneumonia,atb_other,ncol=1,common.legend = TRUE)


class(followup_atb$received_atb)
levels(followup_atb$received_atb)

myColors <- brewer.pal(12,"Set3")
names(myColors) <- levels(followup_atb$received_atb)
names(myColors)
colScale <- scale_colour_manual(name = "received_atb",values = myColors)
