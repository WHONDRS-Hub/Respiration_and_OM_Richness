#Loading packages 
require(reshape2); require(ggplot2); require(easycsv)
library(TPD); library(ggplot2); library(reshape2);library(vegan)
library(ggthemes);library(ggpubr);library(dplyr)
library(stringr);library(e1071) 
library(lattice)
library(cowplot)
library (ggpmisc)

#devtools::install_github("thomasp85/patchwork")

library(ggplot2)
library(patchwork)

#options(digits=10)
rm(list=ls());graphics.off()

#####################################################################################
# Set directories
setwd("//PNL/Projects/SBR_SFA/RC4/Respiration_ms/Input/")
output.dir = ("//PNL/Projects/SBR_SFA/RC4/Respiration_ms/Output/")


# read in the data

respiration.data = read.csv("WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv")
respiration.mass = read.csv("WHONDRS_S19S_Respiration_IncubationData.csv")
field.data = read.csv("Processed_ICR_Field_Hawkes_Peaks.csv")
npoc.data = read.csv("WHONDRS_S19S_Sediment_NPOC.csv")
npoc.to.filter = read.csv("NPOC_ordered_to_filter.csv")
pts.to.remove = 30 # this is an output from 07_Figure_S1 code
respiration.data = merge(respiration.data,respiration.mass, by = "Sample_ID")

respiration.data$rate_mass = respiration.data$rate_mg_per_L_per_h/respiration.data$Sed_mass_g

# Sample ID formatting

field.data$Sample_ID = field.data$Sample_ID_ICR
field.data$Sample_ID = gsub("_p15","",field.data$Sample_ID)
field.data$Sample_ID = gsub("_p2","",field.data$Sample_ID)
field.data$Sample_ID = gsub("_p1","",field.data$Sample_ID)

###############################################################
# Removing the samples in the flagged NPOC list to filter
npoc.to.filter$Sample_ID = paste0(substr(npoc.to.filter$Sample,1,9),"_Sed_Field_ICR",substr(npoc.to.filter$Sample,10,11))

filter.2.threshold = npoc.to.filter[c(1:pts.to.remove),] # List of samples to be removed. Remove the first 30 pts with high ratios. The number 30 came as an output from the Figure_S1 code

# Remove samples from the field data that don't meet the threshold requirement (Ratio between NPOC Field/NPOC INC larger than 2)
all.field.new = field.data
  #all.data.v2
#all.field.new$Sample_ID = all.field.new$ID
for (i in 1:nrow(all.field.new)){
  if (length(which(filter.2.threshold$Sample_ID == all.field.new$Sample_ID[i])) > 0){
    print (i)
    temp = all.field.new[-grep(all.field.new$Sample_ID[i],all.field.new$Sample_ID),]
    all.field.new = temp
  }
  
}

for (i in 1:nrow(all.field.new)){
  if (length(which(filter.2.threshold$Sample_ID == all.field.new$Sample_ID[i])) > 0){
    print (i)
    temp = all.field.new[-grep(all.field.new$Sample_ID[i],all.field.new$Sample_ID),]
    all.field.new = temp
  }
  
}
################################
field.data = all.field.new

field.data$Sample_ID = gsub("Sed","SED",field.data$Sample_ID)
field.data$Sample_ID = gsub("Field_ICR","INC",field.data$Sample_ID)
field.data$Category = "Field"



all.field = merge(respiration.data,field.data, by = "Sample_ID")

all.field$ID = substr(all.field$Sample_ID_ICR,1,25)

names(npoc.data)[2] = "ID"

all.data.v2 = merge(all.field,npoc.data, by = "ID")

all.data.v2$overNPOC = 1/all.data.v2$X00681_NPOC_mg_per_L_as_C
##############################################################################

# Density plots

scaleFUN <- function(x) sprintf("%.2f", x)

p1 =  ggplot(all.data.v2, aes(x = Total_number_of_peaks))+
    geom_density(aes(fill = Category), alpha = 0.5)+ xlim(2180,5000)+
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    theme_bw() + theme(legend.title = element_blank(), 
                       text = element_text(size=12, color="black"),
                       axis.text = element_text(color = "black"),
                       axis.ticks = element_line(color = "black"),
                       panel.background = element_blank(),
                      aspect.ratio = 1,
                       legend.position = "",
                       panel.grid = element_blank())+
                      labs(x = "OM Richness", y = "")+
  geom_text(x=2180,y=0.00125, label="C", size = 5)+ 
  scale_fill_manual(values=  c("#32287d", "#d95f02"))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p2 = ggplot(all.data.v2, aes(x = rate_mg_per_L_per_h))+
  geom_density(alpha = 0.5)+ xlim(-2,42)+
  #scale_y_continuous(labels=scaleFUN) +
  theme_bw() + theme(legend.title = element_blank(),
                     text = element_text(size=12, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     aspect.ratio = 1,
                     panel.grid = element_blank())+ 
                    labs(x = expression(atop("Respiration rate",(mg~L^{-1}~h^{-1}))))+                scale_fill_brewer(palette="Dark2")+
  geom_text(x=-2,y=0.1, label="A", size = 5)+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p3 = ggplot(all.data.v2, aes(x = rate_mass))+
  geom_density(alpha = 0.5)+ xlim(-0.25,3)+
  theme_bw() + theme(legend.title = element_blank(),
                     text = element_text(size=12, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     aspect.ratio = 1,
                     panel.grid = element_blank())+ 
  #+ scale_y_continuous(labels=scaleFUN)+
  labs(x = expression(atop("Respiration rate",(mg~L^{-1}~h^{-1}~g^{-1}))), y = "")+                  scale_fill_brewer(palette="Dark2")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
geom_text(x=-0.25,y=1.95, label="B", size = 5)
  

p = p2+p3+p1
# p = plot_grid(p2,p3,p1, nrow = 1,align = "hv")

ggsave(paste(output.dir,"Figure_2_",Sys.Date(),".pdf", sep = ""),p, dpi = 500)

ggsave(paste(output.dir,"Figure_2_",Sys.Date(),".png", sep = ""),p, dpi = 500)
