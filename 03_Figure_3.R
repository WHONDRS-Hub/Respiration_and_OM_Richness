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
output.dir = ("Output/")


# read in the data

respiration.data = read.csv("Data/WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv")
respiration.mass = read.csv("Data/WHONDRS_S19S_Respiration_IncubationData.csv")
field.data = read.csv("Data/Processed_ICR_Field_Hawkes_Peaks.csv")
npoc.data = read.csv("Data/WHONDRS_S19S_Sediment_NPOC.csv")
npoc.to.filter = read.csv("Data/NPOC_ordered_to_filter.csv")
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

#############################################################################
# Figure 3. Regressions
##############################################################################
ideal.num.segments = 10

#####################################################################
all.data.v2 = all.data.v2[order(all.data.v2$Total_number_of_peaks),]

field.step = (max(all.data.v2$Total_number_of_peaks)- min(all.data.v2$Total_number_of_peaks))/(ideal.num.segments)

peak.segments = as.data.frame(matrix(NA,ncol = 4,nrow = ideal.num.segments))
colnames(peak.segments)= c("Low.boundary","High.boundary","Max.resp","Mid.peaks")

temp.low = min(all.data.v2$Total_number_of_peaks)


for (i in 1:ideal.num.segments){
  
  if (i == 10){
    
    temp.high = temp.low + field.step
    temp.high = plyr::round_any(temp.high, 0.0001, f = ceiling)
    temp.dat = all.data.v2[which(all.data.v2$Total_number_of_peaks >= temp.low & all.data.v2$Total_number_of_peaks <= temp.high),]
    
  }else{
    temp.high = temp.low + field.step
    temp.dat = all.data.v2[which(all.data.v2$Total_number_of_peaks >= temp.low & all.data.v2$Total_number_of_peaks <= temp.high),]
  }
  
  if (nrow(temp.dat) > 0){
    temp.rate = max(all.data.v2$rate_mg_per_L_per_h[which(all.data.v2$Total_number_of_peaks >= temp.low & all.data.v2$Total_number_of_peaks <= temp.high)])
    
    
    temp.mid =  temp.dat$Total_number_of_peaks[which(temp.dat$rate_mg_per_L_per_h %in% temp.rate)]
    
    peak.segments$Low.boundary[i] = temp.low 
    peak.segments$High.boundary[i] = temp.high
    peak.segments$Max.resp[i] = temp.rate
    peak.segments$Mid.peaks[i] = temp.mid
    
    temp.low = temp.high  
    
  } else if (nrow(temp.dat) == 0){
    temp.rate = NA
    
    
    temp.mid =  NA
    
    peak.segments$Low.boundary[i] = temp.low 
    peak.segments$High.boundary[i] = temp.high
    peak.segments$Max.resp[i] = temp.rate
    peak.segments$Mid.peaks[i] = temp.mid
    
    temp.low = temp.high  
    
  }
}

## Plot quadratic regression

x = peak.segments$Mid.peaks
y = peak.segments$Max.resp

label.x =  expression(OM~Richness)
label.y =  expression(Respiration~rate~(mg~L^{-1}~h^{-1}))

dfi=as.data.frame(cbind(y,x))
names(dfi)[1]="y";names(dfi)[2]="x"

dfi=na.omit(dfi)

mod2 = summary(lm(y ~ poly(x, 2)))

mod3 = summary(lm(all.data.v2$rate_mg_per_L_per_h ~ all.data.v2$Total_number_of_peaks))


#Plots
library (ggpmisc)
my.formula <- y ~ poly(x, 2)


p1 = ggplot(dfi, aes(x = x, y = y)) + coord_cartesian(ylim = c(0,60))+
  geom_point(size = 3,color = "#32287d") + 
  geom_point(data=all.data.v2, aes(x=Total_number_of_peaks, y= rate_mg_per_L_per_h), size = 2, color= "#a29bd5") +
  geom_point(size = 3,color = "#32287d") + 
  geom_smooth(method = "lm", se=F, formula = my.formula,color = "#32287d") +
  # geom_vline(xintercept = peak.segments$Low.boundary[1])+
  # geom_vline(xintercept = peak.segments$High.boundary)+
  geom_text(x=2400, y=60, label="A", size = 5)+
  geom_text(x=4500, y=60, label=paste0("R^2 == ",round(mod2$r.squared,2)), color = "#32287d", parse = T)+
   #geom_text(x=4500, y=55, label=paste0("p = ",round(mod2$coefficients[2,4],2)), color = "#32287d")+
  geom_text(x=4500, y=55, label="p = 0.16 ",color = "#32287d")+
  geom_text(x=3000, y=60, label=paste0("R^2 == ",round(mod3$r.squared,2)), color = "#a29bd5", parse = T)+
  geom_text(x=3000, y=55, label=paste0("p = 0.05"), color = "#a29bd5")+
  geom_abline(slope = mod3$coefficients[2,1]
, intercept = mod3$coefficients[1,1], color = "#a29bd5", size = 1)+
  theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=12,face="bold"))+labs(y = label.y, x = label.x)+ 
  #theme(legend.position=c(0.7,0.2))+ 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.y=element_text(colour = c("black","black")))+
  theme(aspect.ratio=1)+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x =element_text(size=12,face="bold"))+
  theme(axis.title =element_text(size=12,face="bold"))+
  theme(axis.title.y =element_text(size=12,face="bold"))
#+
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 

####################################
# rates normalized
all.data.v2 = all.data.v2[order(all.data.v2$Total_number_of_peaks),]

field.step = (max(all.data.v2$Total_number_of_peaks)- min(all.data.v2$Total_number_of_peaks))/(ideal.num.segments)

peak.segments = as.data.frame(matrix(NA,ncol = 4,nrow = ideal.num.segments))
colnames(peak.segments)= c("Low.boundary","High.boundary","Max.resp","Mid.peaks")

temp.low = min(all.data.v2$Total_number_of_peaks)


for (i in 1:ideal.num.segments){
  
  if (i == 10){
    
    temp.high = temp.low + field.step
    temp.high = plyr::round_any(temp.high, 0.0001, f = ceiling)
    temp.dat = all.data.v2[which(all.data.v2$Total_number_of_peaks >= temp.low & all.data.v2$Total_number_of_peaks <= temp.high),]
    
  }else{
    temp.high = temp.low + field.step
    temp.dat = all.data.v2[which(all.data.v2$Total_number_of_peaks >= temp.low & all.data.v2$Total_number_of_peaks <= temp.high),]
  }
  
  if (nrow(temp.dat) > 0){
    temp.rate = max(all.data.v2$rate_mass[which(all.data.v2$Total_number_of_peaks >= temp.low & all.data.v2$Total_number_of_peaks <= temp.high)])
    
    
    temp.mid =  temp.dat$Total_number_of_peaks[which(temp.dat$rate_mass %in% temp.rate)]
    
    peak.segments$Low.boundary[i] = temp.low 
    peak.segments$High.boundary[i] = temp.high
    peak.segments$Max.resp[i] = temp.rate
    peak.segments$Mid.peaks[i] = temp.mid
    
    temp.low = temp.high  
    
  } else if (nrow(temp.dat) == 0){
    temp.rate = NA
    
    
    temp.mid =  NA
    
    peak.segments$Low.boundary[i] = temp.low 
    peak.segments$High.boundary[i] = temp.high
    peak.segments$Max.resp[i] = temp.rate
    peak.segments$Mid.peaks[i] = temp.mid
    
    temp.low = temp.high  
    
  }
}

## Plot quadratic regression

x = peak.segments$Mid.peaks
y = peak.segments$Max.resp

label.x =  expression(OM~Richness)
label.y =  expression(Respiration~rate~(mg~L^{-1}~h^{-1}~g^{-1}))

dfi=as.data.frame(cbind(y,x))
names(dfi)[1]="y";names(dfi)[2]="x"

dfi=na.omit(dfi)

mod22 = summary(lm(y ~ poly(x, 2)))

mod33 = summary(lm(all.data.v2$rate_mass ~ all.data.v2$Total_number_of_peaks))


#Plots
library (ggpmisc)
my.formula <- y ~ poly(x, 2)


p11 = ggplot(dfi, aes(x = x, y = y)) + coord_cartesian(ylim = c(0,2.5))+
  geom_point(size = 3,color = "#32287d") + 
  geom_point(data=all.data.v2, aes(x=Total_number_of_peaks, y= rate_mass ), size = 2, color= "#a29bd5") +
  geom_point(size = 3,color = "#32287d") + 
  geom_smooth(method = "lm", se=F, formula = my.formula,color = "#32287d") +
  # geom_vline(xintercept = peak.segments$Low.boundary[1])+
  # geom_vline(xintercept = peak.segments$High.boundary)+
  geom_text(x=2400, y=2.5, label="B", size = 5)+
  geom_text(x=4500, y=2.5, label=paste0("R^2 == ",round(mod22$r.squared,2)), color = "#32287d", parse = T)+
  #geom_text(x=4500, y=55, label=paste0("p = ",round(mod2$coefficients[2,4],2)), color = "#32287d")+
  geom_text(x=4500, y=2.3, label="p = 0.13 ",color = "#32287d")+
  geom_text(x=3000, y=2.5, label=paste0("R^2 == ",round(mod33$r.squared,2)), color = "#a29bd5", parse = T)+
  geom_text(x=3000, y=2.3, label=paste0("p = 0.07"), color = "#a29bd5")+
  geom_abline(slope = mod33$coefficients[2,1]
              , intercept = mod33$coefficients[1,1], color = "#a29bd5", size = 1)+
  theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=12,face="bold"))+labs(y = label.y, x = label.x)+ 
  #theme(legend.position=c(0.7,0.2))+ 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.y=element_text(colour = c("black","black")))+
  theme(aspect.ratio=1)+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x =element_text(size=12,face="bold"))+
  theme(axis.title =element_text(size=12,face="bold"))+
  theme(axis.title.y =element_text(size=12,face="bold"))
#+
 # theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 


 

 p = (p1|p11)
 
 ggsave(paste(output.dir,"Figure_3_",Sys.Date(),".pdf", sep = ""),p, dpi = 500)
 ggsave(paste(output.dir,"Figure_3_",Sys.Date(),".png", sep = ""),p, dpi = 500)
 
 write.csv(all.data.v2,paste0(output.dir,"all_data.csv"))