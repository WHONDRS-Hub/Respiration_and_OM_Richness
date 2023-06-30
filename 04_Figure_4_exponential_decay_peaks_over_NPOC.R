#Loading packages 
require(reshape2); require(ggplot2); require(easycsv)
library(TPD); library(ggplot2); library(reshape2);library(vegan)
library(ggthemes);library(ggpubr);library(dplyr)
library(stringr);library(e1071) 
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(cowplot)

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
#####Field Data################################################
# Sample ID formatting

field.data$Sample_ID = field.data$Sample_ID_ICR
field.data$Sample_ID = gsub("_p15","",field.data$Sample_ID)
field.data$Sample_ID = gsub("_p2","",field.data$Sample_ID)
field.data$Sample_ID = gsub("_p1","",field.data$Sample_ID)
field.data$Sample_ID = gsub("Sed","SED",field.data$Sample_ID)
field.data$Sample_ID = gsub("Field_ICR","INC",field.data$Sample_ID)
field.data$Category = "Field"
field.all = merge(field.data,respiration.data, by = "Sample_ID")
field.all$Sample_ID = field.all$Sample_ID_ICR
field.all$Sample_ID = gsub("_p15","",field.all$Sample_ID)
field.all$Sample_ID = gsub("_p2","",field.all$Sample_ID)
field.all$Sample_ID = gsub("_p1","",field.all$Sample_ID)

all.field = merge(field.all,npoc.data, by = "Sample_ID")


all.field$peaksovernpoc = all.field$Total_number_of_peaks/all.field$X00681_NPOC_mg_per_L_as_C

###############################################################
# Removing the samples in the flagged NPOC list to filter
npoc.to.filter$Sample_ID = paste0(substr(npoc.to.filter$Sample,1,9),"_Sed_Field_ICR",substr(npoc.to.filter$Sample,10,11))

filter.2.threshold = npoc.to.filter[c(1:pts.to.remove),] # List of samples to be removed. Remove the first 30 pts with high ratios. The number 30 came as an output from the Figure_S1 code


all.field.new = all.field
for (i in 1:nrow(all.field.new)){
  if (length(which(filter.2.threshold$Sample_ID == all.field.new$Sample_ID[i])) > 0){
    print (i)
    temp = all.field.new[-grep(all.field.new$Sample_ID[i],all.field.new$Sample_ID),]
    all.field.new = temp
  }

}
all.field.old = all.field
all.field = all.field.new

#####################################################################
#### Defining the segments
###################
ideal.num.segments = 10

#####################################################################
# Number of peaks normalized by NPOC
  
  all.field.peaks = all.field[order(all.field$peaksovernpoc),]

  field.step = (max(all.field.peaks$peaksovernpoc)-min(all.field.peaks$peaksovernpoc))/(ideal.num.segments)
  
  peak.segments = as.data.frame(matrix(NA,ncol = 4,nrow = ideal.num.segments))
  colnames(peak.segments)= c("Low.boundary","High.boundary","Max.resp","Mid.peaks")
  
  temp.low = (min(all.field.peaks$peaksovernpoc))
  
  
  for (i in 1:ideal.num.segments){

    if (i == 10){
     
      temp.high = temp.low + field.step
      temp.high = plyr::round_any(temp.high, 0.0001, f = ceiling)
      temp.dat = all.field.peaks[which(all.field.peaks$peaksovernpoc >= temp.low & all.field.peaks$peaksovernpoc <= temp.high),]
      
    }else{
    temp.high = temp.low + field.step
    temp.dat = all.field.peaks[which(all.field.peaks$peaksovernpoc >= temp.low & all.field.peaks$peaksovernpoc <= temp.high),]
    }
    
    if (nrow(temp.dat) > 0){
      temp.rate = max(all.field.peaks$rate_mg_per_L_per_h[which(all.field.peaks$peaksovernpoc >= temp.low & all.field.peaks$peaksovernpoc <= temp.high)])
      
    
      temp.mid =  temp.dat$peaksovernpoc[which(temp.dat$rate_mg_per_L_per_h %in% temp.rate)]
      
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
  
  
  
  ## Plot exp relationships
  
  x = peak.segments$Mid.peaks
  y = peak.segments$Max.resp
  
  label.x =  expression(OM~Richness/NPOC~(mg~C~L^{-1}))
  label.y =  expression(Respiration~rate~(mg~L^{-1}~h^{-1}))
  
  dfi=as.data.frame(cbind(y,x))
  names(dfi)[1]="y";names(dfi)[2]="x"
  
  dfi=na.omit(dfi)
  
  mod = summary(lm(log(y)~x))
  
  #Plots
  library (ggpmisc)
  my.formula <- y ~ exp(mod$coefficients[1,1] + (mod$coefficients[2,1]*x))
  
  # Calculating stats for the regressions that use all the points
  mod11 = summary(lm(log(rate_mg_per_L_per_h)~peaksovernpoc, data = all.field.peaks))
  
p1 = ggplot(dfi, aes(x = x, y = y)) + #coord_cartesian(ylim = c(-5,82))+
    geom_point(size = 3,color = "#32287d") + 
    geom_point(data=all.field.peaks, aes(x=peaksovernpoc, y= rate_mg_per_L_per_h), size = 2, color= "#a29bd5") +
    geom_point(size = 3,color = "#32287d") + 
  geom_smooth(method = "lm", se=TRUE, formula = my.formula) +
    geom_text(x=80, y=63, label="A", size = 5)+
    geom_text(x=985, y=61.7, label=paste0("R^2 == ",round(mod$r.squared,2)), color = "#32287d", parse = T)+
    # geom_text(x=985, y=75, label=paste0("p = ",round(mod$coefficients[2,4],8)), color = "#32287d")+
 # geom_vline(xintercept = peak.segments$Low.boundary[1])+
  #geom_vline(xintercept = peak.segments$High.boundary)+
    geom_text(x=985, y=58, label="p < 0.001 ",color = "#32287d")+
    theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=12,face="bold"))+labs(y = label.y, x = label.x)+ 
    #theme(legend.position=c(0.7,0.2))+ 
    theme(axis.text.x=element_text(size=12))+theme(axis.text.x=element_text(colour = c("black","black")))+
    theme(aspect.ratio=1)+
    theme(axis.text.y=element_text(size=12))+
    theme(axis.title.x =element_text(size=12,face="bold"))+
    theme(axis.title =element_text(size=12,face="bold"))+
    theme(axis.title.y =element_text(size=12,face="bold"))
#+
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
  # ggsave(paste0(output.dir,regression.folder,"Field_Peaks_over_NPOC_Max_respiration_rate_",ideal.num.segments,"_",Sys.Date(),".png"),width = 7.38, height = 7.38)
  # ggsave(paste0(output.dir,regression.folder,"Field_Peaks_over_NPOC_Max_respiration_rate_",ideal.num.segments,"_",Sys.Date(),".pdf"),width = 7.38, height = 7.38)
  # 
##########################################################################
#########################################################################
# Normalized respiration 

all.field.peaks = all.field[order(all.field$peaksovernpoc),]

field.step = (max(all.field.peaks$peaksovernpoc)-min(all.field.peaks$peaksovernpoc))/(ideal.num.segments)

peak.segments = as.data.frame(matrix(NA,ncol = 4,nrow = ideal.num.segments))
colnames(peak.segments)= c("Low.boundary","High.boundary","Max.resp","Mid.peaks")

temp.low = (min(all.field.peaks$peaksovernpoc))


for (i in 1:ideal.num.segments){
  
  if (i == 10){
    
    temp.high = temp.low + field.step
    temp.high = plyr::round_any(temp.high, 0.0001, f = ceiling)
    temp.dat = all.field.peaks[which(all.field.peaks$peaksovernpoc >= temp.low & all.field.peaks$peaksovernpoc <= temp.high),]
    
  }else{
    temp.high = temp.low + field.step
    temp.dat = all.field.peaks[which(all.field.peaks$peaksovernpoc >= temp.low & all.field.peaks$peaksovernpoc <= temp.high),]
  }
  
  if (nrow(temp.dat) > 0){
    temp.rate = max(all.field.peaks$rate_mass[which(all.field.peaks$peaksovernpoc >= temp.low & all.field.peaks$peaksovernpoc <= temp.high)])
    
    
    temp.mid =  temp.dat$peaksovernpoc[which(temp.dat$rate_mass %in% temp.rate)]
    
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



## Plot exp relationships

x = peak.segments$Mid.peaks
y = peak.segments$Max.resp

label.x =  expression(OM~Richness/NPOC~(mg~C~L^{-1}))
label.y =  expression(Respiration~rate~(mg~L^{-1}~h^{-1}~g^{-1}))

dfi=as.data.frame(cbind(y,x))
names(dfi)[1]="y";names(dfi)[2]="x"

dfi=na.omit(dfi)

mod2 = summary(lm(log(y)~x))

# Calculating stats for the regressions that use all the points
mod22 = summary(lm(log(rate_mass)~peaksovernpoc, data = all.field.peaks))

#Plots
library (ggpmisc)
my.formula2 <- y ~ exp(mod2$coefficients[1,1] + (mod2$coefficients[2,1]*x))


p2 = ggplot(dfi, aes(x = x, y = y)) + #coord_cartesian(ylim = c(-5,82))+
  geom_point(size = 3,color = "#32287d") + 
  geom_point(data=all.field.peaks, aes(x=peaksovernpoc, y= rate_mass), size = 2, color= "#a29bd5") +
  geom_point(size = 3,color = "#32287d") + 
  geom_smooth(method = "lm", se=TRUE, formula = my.formula2) +
  geom_text(x=80, y=3.3, label="B", size = 5)+
  geom_text(x=985, y=3.1, label=paste0("R^2 == ",round(mod2$r.squared,2)), color = "#32287d", parse = T)+
  # geom_text(x=985, y=75, label=paste0("p = ",round(mod$coefficients[2,4],8)), color = "#32287d")+
  # geom_vline(xintercept = peak.segments$Low.boundary[1])+
  #geom_vline(xintercept = peak.segments$High.boundary)+
  geom_text(x=985, y=2.9, label="p < 0.001 ",color = "#32287d")+
  theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=12,face="bold"))+labs(y = label.y, x = label.x)+ 
  #theme(legend.position=c(0.7,0.2))+ 
  theme(axis.text.x=element_text(size=12))+theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(aspect.ratio=1)+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x =element_text(size=12,face="bold"))+
  theme(axis.title =element_text(size=12,face="bold"))+
  theme(axis.title.y =element_text(size=12,face="bold"))
  p = plot_grid(p1,p2)
  
  ggsave(paste(output.dir,"Figure_4_",Sys.Date(),".pdf", sep = ""),p, dpi = 500)
  
  ggsave(paste(output.dir,"Figure_4_",Sys.Date(),".png", sep = ""),p)
  
  