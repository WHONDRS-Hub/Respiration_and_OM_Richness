#Loading packages 
require(reshape2); require(ggplot2); require(easycsv)
library(TPD); library(ggplot2); library(reshape2);library(vegan)
library(ggthemes);library(ggpubr);library(dplyr)
library(stringr);library(e1071)  

#options(digits=10)
rm(list=ls());graphics.off()

#####################################################################################
# Set directories
setwd("//PNL/Projects/SBR_SFA/RC4/Respiration_ms/Input/")
output.dir = ("//PNL/Projects/SBR_SFA/RC4/Respiration_ms/Output/")
regression.folder = "Poly_regressions/"


# read in the data

respiration.data = read.csv("WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv")
field.data = read.csv("Processed_ICR_Field_Peaks.csv")
incubation.data = read.csv("Processed_ICR_INC_Peaks.csv")
npoc.data = read.csv("WHONDRS_S19S_Sediment_NPOC.csv")


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


#####################################################################
#### Defining the segments
###################
ideal.num.segments = 10

#####################################################################
# Number of peaks normalized by NPOC
  
  all.field.peaks = all.field[order(all.field$peaksovernpoc),]

  field.step = (max(all.field.peaks$peaksovernpoc)- min(all.field.peaks$peaksovernpoc))/(ideal.num.segments)
  
  peak.segments = as.data.frame(matrix(NA,ncol = 4,nrow = ideal.num.segments))
  colnames(peak.segments)= c("Low.boundary","High.boundary","Max.resp","Mid.peaks")
  
  temp.low = min(all.field.peaks$peaksovernpoc)
  
  
  for (i in 1:ideal.num.segments){
    
    temp.high = temp.low + field.step
    temp.dat = all.field.peaks[which(all.field.peaks$peaksovernpoc >= temp.low & all.field.peaks$peaksovernpoc <= temp.high),]
    
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
  
  
  
  ## Plot polynomial relationships
  
  x = peak.segments$Mid.peaks
  y = peak.segments$Max.resp
  
  label.x = "Peaks_over_NPOC"
  label.y = "Max_Resp_rate_mg_per_L_per_h"
  
  dfi=as.data.frame(cbind(y,x))
  names(dfi)[1]="y";names(dfi)[2]="x"
  
  dfi=na.omit(dfi)
  
  
  #Plots
  library (ggpmisc)
  my.formula <- y ~ poly(x, 2, raw=TRUE)
  
  ggplot(dfi, aes(x = x, y = y)) + coord_cartesian(xlim = c(80,1100), ylim = c(-5,82))+
    geom_point(size = 7,color = "#32287d") + 
    geom_point(data=all.field.peaks, aes(x=peaksovernpoc, y= rate_mg_per_L_per_h), size = 5, color= "#a29bd5") +
    geom_point(size = 7,color = "#32287d") + 
    geom_smooth(method = "lm", se=TRUE, formula = my.formula) +
    stat_poly_eq(formula = my.formula,label.y = "top",label.x = "left", aes(label = paste( ..rr.label..,..eq.label.., sep = "~~~"),size=4),
                 parse = TRUE)+stat_fit_glance(data=dfi, method = 'lm',
                                               method.args = list(formula = my.formula),
                                               geom = 'text',
                                               aes(label =paste("p = ",signif(..p.value.., digits = 2), sep = ""),size=4),label.y = "top",label.x = "right")+
    theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=12,face="bold"))+labs(y = label.y, x = label.x)+ 
    #theme(legend.position=c(0.7,0.2))+ 
    theme(axis.text.x=element_text(size=12,face="bold"))+theme(axis.text.x=element_text(colour = c("black","black")))+
    #theme(aspect.ratio=1)+
    theme(axis.text.y=element_text(size=12,face="bold"))+
    theme(axis.title.x =element_text(size=12,face="bold"))+
    theme(axis.title =element_text(size=12,face="bold"))+
    theme(axis.title.y =element_text(size=12,face="bold"))
  ggsave(paste0(output.dir,regression.folder,"Field_Peaks_over_NPOC_Max_respiration_rate_",ideal.num.segments,"_",Sys.Date(),".png"),width = 7.38, height = 7.38)
  ggsave(paste0(output.dir,regression.folder,"Field_Peaks_over_NPOC_Max_respiration_rate_",ideal.num.segments,"_",Sys.Date(),".pdf"),width = 7.38, height = 7.38)
  
##########################################################################
#########################################################################
  #Incubation data
  
  # Sample ID formatting
  
  incubation.data$Sample_ID = incubation.data$Sample_ID_ICR
  incubation.data$Sample_ID = gsub("_p15","",incubation.data$Sample_ID)
  incubation.data$Sample_ID = gsub("_p2","",incubation.data$Sample_ID)
  incubation.data$Sample_ID = gsub("_p1","",incubation.data$Sample_ID)
  incubation.data$Sample_ID = gsub("Sed","SED",incubation.data$Sample_ID)
  incubation.data$Sample_ID = gsub("INC_ICR","INC",incubation.data$Sample_ID)
  incubation.data$Category = "Incubation"
  
  incubation.all = merge(incubation.data,respiration.data, by = "Sample_ID")
  incubation.all$Sample_ID = incubation.all$Sample_ID_ICR
  incubation.all$Sample_ID = gsub("_p15","",incubation.all$Sample_ID)
  incubation.all$Sample_ID = gsub("_p2","",incubation.all$Sample_ID)
  incubation.all$Sample_ID = gsub("_p1","",incubation.all$Sample_ID)
  
  all.incubation = merge(incubation.all,npoc.data, by = "Sample_ID")
  

  all.incubation$peaksovernpoc = all.incubation$Total_number_of_peaks/all.incubation$X00681_NPOC_mg_per_L_as_C
  
#####################################################################
#### Defining the segments
###################
  ideal.num.segments = 10
##################################################################
#####################################################################
  # Number of peaks normalized by NPOC
  
  all.incubation.peaks = all.incubation[order(all.incubation$peaksovernpoc),]
  incubation.step = (max(all.incubation.peaks$peaksovernpoc)- min(all.incubation.peaks$peaksovernpoc))/(ideal.num.segments)
  
  peak.segments = as.data.frame(matrix(NA,ncol = 4,nrow = ideal.num.segments))
  colnames(peak.segments)= c("Low.boundary","High.boundary","Max.resp","Mid.peaks")
  
  temp.low = min(all.incubation.peaks$peaksovernpoc)
  
  
  for (i in 1:ideal.num.segments){
    
    temp.high = temp.low + incubation.step
    temp.dat = all.incubation.peaks[which(all.incubation.peaks$peaksovernpoc >= temp.low & all.incubation.peaks$peaksovernpoc <= temp.high),]
    
    if (nrow(temp.dat) > 0){
      temp.rate = max(all.incubation.peaks$rate_mg_per_L_per_h[which(all.incubation.peaks$peaksovernpoc >= temp.low & all.incubation.peaks$peaksovernpoc <= temp.high)])
      
      #temp.mid = temp.low + (incubation.step/2)
      temp.mid =  temp.dat$peaksovernpoc[which(temp.dat$rate_mg_per_L_per_h %in% temp.rate)]
      
      peak.segments$Low.boundary[i] = temp.low 
      peak.segments$High.boundary[i] = temp.high
      peak.segments$Max.resp[i] = temp.rate
      peak.segments$Mid.peaks[i] = temp.mid
      
      temp.low = temp.high  
      
    } else if (nrow(temp.dat) == 0){
      temp.rate = NA
      
      #temp.mid = temp.low + (incubation.step/2)
      temp.mid =  NA
      
      peak.segments$Low.boundary[i] = temp.low 
      peak.segments$High.boundary[i] = temp.high
      peak.segments$Max.resp[i] = temp.rate
      peak.segments$Mid.peaks[i] = temp.mid
      
      temp.low = temp.high  
      
    }
  }
  
  
  
  ## Plot polynomial relationships
  
  x = peak.segments$Mid.peaks
  y = peak.segments$Max.resp
  
  label.x = "Peaks_over_NPOC"
  label.y = "Max_Resp_rate_mg_per_L_per_h"
  
  dfi=as.data.frame(cbind(y,x))
  names(dfi)[1]="y";names(dfi)[2]="x"
  
  dfi=na.omit(dfi)
  
  
  #Plots
  library (ggpmisc)
  my.formula <- y ~ poly(x, 2, raw=TRUE)
  
  ggplot(dfi, aes(x = x, y = y)) + coord_cartesian(xlim = c(80,1100), ylim = c(-5,82))+
    geom_point(size = 7, color = "#d95f02") +
    geom_point(data=all.incubation.peaks, aes(x=peaksovernpoc, y= rate_mg_per_L_per_h), size = 5, color= '#f2b78a') +
    geom_point(size = 7, color = "#d95f02") +
     #geom_vline(data = peak.segments, aes(xintercept = as.numeric(Low.boundary)))+
     #geom_vline(data = peak.segments, aes(xintercept = as.numeric(High.boundary[ideal.num.segments])))+
    geom_smooth(method = "lm", se=TRUE, formula = my.formula) +
    stat_poly_eq(formula = my.formula,label.y = "top",label.x = "left", aes(label = paste( ..rr.label..,..eq.label.., sep = "~~~"),size=4),
                 parse = TRUE)+stat_fit_glance(data=dfi, method = 'lm',
                                               method.args = list(formula = my.formula),
                                               geom = 'text',
                                               aes(label =paste("p = ",signif(..p.value.., digits = 2), sep = ""),size=4),label.y = "top",label.x = "right")+
    theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=12,face="bold"))+labs(y = label.y, x = label.x)+ 
    #theme(legend.position=c(0.7,0.2))+ 
    theme(axis.text.x=element_text(size=12,face="bold"))+theme(axis.text.x=element_text(colour = c("black","black")))+
    #theme(aspect.ratio=1)+
    theme(axis.text.y=element_text(size=12,face="bold"))+
    theme(axis.title.x =element_text(size=12,face="bold"))+
    theme(axis.title =element_text(size=12,face="bold"))+
    theme(axis.title.y =element_text(size=12,face="bold"))
  ggsave(paste0(output.dir,regression.folder,"Incubation_Peaks_over_NPOC_Max_respiration_rate_",ideal.num.segments,"_",Sys.Date(),".png"),width = 7.38, height = 7.38)
  ggsave(paste0(output.dir,regression.folder,"Incubation_Peaks_over_NPOC_Max_respiration_rate_",ideal.num.segments,"_",Sys.Date(),".pdf"),width = 7.38, height = 7.38)
  