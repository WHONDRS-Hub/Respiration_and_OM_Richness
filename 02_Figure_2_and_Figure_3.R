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
density.folder = "Density/"
regression.folder = "Regressions/"
histogram.folder = "Histogram/"

# read in the data

respiration.data = read.csv("WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv")
field.data = read.csv("Processed_ICR_Field_Peaks.csv")
incubation.data = read.csv("Processed_ICR_INC_Peaks.csv")

# Sample ID formatting

field.data$Sample_ID = field.data$Sample_ID_ICR
field.data$Sample_ID = gsub("_p15","",field.data$Sample_ID)
field.data$Sample_ID = gsub("_p2","",field.data$Sample_ID)
field.data$Sample_ID = gsub("_p1","",field.data$Sample_ID)
field.data$Sample_ID = gsub("Sed","SED",field.data$Sample_ID)
field.data$Sample_ID = gsub("Field_ICR","INC",field.data$Sample_ID)
field.data$Category = "Field"


incubation.data$Sample_ID = incubation.data$Sample_ID_ICR
incubation.data$Sample_ID = gsub("_p15","",incubation.data$Sample_ID)
incubation.data$Sample_ID = gsub("_p2","",incubation.data$Sample_ID)
incubation.data$Sample_ID = gsub("_p1","",incubation.data$Sample_ID)
incubation.data$Sample_ID = gsub("_p3","",incubation.data$Sample_ID)
incubation.data$Sample_ID = gsub("_ICR","",incubation.data$Sample_ID)
incubation.data$Category = "Incubation"

##############################################################################

# Density plots

all.data = rbind(field.data,incubation.data)

  ggplot(all.data, aes(x = Total_number_of_peaks , group = Category))+
    geom_density(aes(fill = Category), alpha = 0.5)+
    theme_bw() + theme(legend.title = element_blank(), 
                       text = element_text(size=12, color="black"),
                       axis.text = element_text(color = "black"),
                       axis.ticks = element_line(color = "black"),
                       panel.background = element_blank(),
                       aspect.ratio = 1,
                       legend.position = "top",
                       panel.grid = element_blank())+
                      labs(x = colnames(all.data[2]))+
                       scale_fill_manual(values=  c("#32287d", "#d95f02")) 
  ggsave(paste0(output.dir,density.folder,"png/",colnames(all.data[2]),"_",Sys.Date(),".png"),width = 7.38, height = 7.38)
  ggsave(paste0(output.dir,density.folder,colnames(all.data[2]),"_",Sys.Date(),".pdf"))



ggplot(respiration.data, aes(x = rate_mg_per_L_per_h))+
  geom_density(alpha = 0.5)+
  theme_bw() + theme(legend.title = element_blank(),
                     text = element_text(size=12, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     aspect.ratio = 1,
                     panel.grid = element_blank())+ 
                    labs(x = "rate_mg_per_L_per_h")+
                    scale_fill_brewer(palette="Dark2") 
ggsave(paste0(output.dir,density.folder,"png/","rate_mg_per_L_per_h","_",Sys.Date(),".png"),width = 7.38, height = 7.38)
ggsave(paste0(output.dir,density.folder,"rate_mg_per_L_per_h","_",Sys.Date(),".pdf"))




#############################################################################
# Figure 3. Regressions
##############################################################################


# Both on the same plot
all.field = merge(respiration.data,field.data, by = "Sample_ID")
all.incubation = merge(respiration.data,incubation.data, by = "Sample_ID")

all.data.v2 = rbind(all.incubation,all.field)

m <- matrix(NA, ncol = 5, nrow=1) #ncol has to be equal to the number of exported variables
regression=data.frame(m)
names(regression)[1]="Peak.metric" ;names(regression)[2]="regression.slope" ; names(regression)[3]="R2.adj" ;names(regression)[4]="R2" ; names(regression)[5]="pvalue"  

  
  x = all.data.v2$Total_number_of_peaks
  y = all.data.v2$rate_mg_per_L_per_h
  z = as.data.frame(all.data.v2$Category)
  
  dfi=cbind(y,x,z)
  name=colnames(dfi)
  names(dfi)[1]="y";names(dfi)[2]="x";names(dfi)[3]="Category"

  dfi=na.omit(dfi)
  
  fit=lm(dfi$y~dfi$x)
  u=fit$coefficients
  b=u[[1]] #Intercept
  c=u[[2]] #regression slope
  #r=summary(lm(dfi$y~dfi$x))$r.squared
  r.adj=summary(lm(dfi$y~dfi$x))$adj.r.squared
  r=summary(lm(dfi$y~dfi$x))$r.squared
  p=summary(fit)$coefficients[4]  
  
  name=names(all.data.v2[10])
  
  
  #This is what makes the plots
  library (ggpmisc)
  my.formula <- y ~ x
  
  ggplot(dfi, aes(x = x, y = y, color = Category)) + coord_cartesian()+
    geom_point(size = 5) + scale_color_manual(values=  c("#32287d", "#d95f02")) +
    geom_smooth(method = "lm", se=TRUE, formula = my.formula) +
    stat_poly_eq(formula = my.formula,label.y = "top",label.x = "left", aes(label = paste( ..rr.label.., sep = "~~~"),size=4),
                 parse = TRUE)+stat_fit_glance(data=dfi, method = 'lm',
                                               method.args = list(formula = my.formula),
                                               geom = 'text',
                                               aes(label =paste("p = ",signif(..p.value.., digits = 2), sep = ""),size=4),label.y = "top",label.x = "right")+
    theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=12,face="bold"))+labs(y = "rate_mg_per_L_per_h", x = colnames(all.data.v2[10]))+ 
    theme(legend.position="top")+ 
    theme(axis.text.x=element_text(size=12,face="bold"))+theme(axis.text.x=element_text(colour = c("black","black")))+
    #theme(aspect.ratio=1)+
    theme(axis.text.y=element_text(size=12,face="bold"))+
    theme(axis.title.x =element_text(size=12,face="bold"))+
    theme(axis.title =element_text(size=12,face="bold"))+
    theme(axis.title.y =element_text(size=12,face="bold"))
  
  ggsave(paste0(output.dir,regression.folder,"Both/png/",colnames(all.data.v2[10]),"_Both_respiration_rate_",Sys.Date(),".png"),width = 7.38, height = 7.38)
  ggsave(paste0(output.dir,regression.folder,"Both/",colnames(all.data.v2[10]),"_Both_respiration_rate_",Sys.Date(),".pdf"),width = 7.38, height = 7.38)
  
  regression = rbind(regression,c(
    name[1],
    as.numeric((c)),
    as.numeric(abs(r.adj)),
    as.numeric(abs(r)),
    p))


# regression=regression[-1,]
# write.csv(regression,paste0(output.dir,regression.folder,"Regression_metrics_Both_ICR_",Sys.Date(),".csv"), row.names = F)
# rm(regression); rm(i); rm(all.data.v2)

