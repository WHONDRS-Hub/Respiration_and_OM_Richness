rm(list=ls());graphics.off()
library(ggplot2)
#####################################################################################
# Set directories
output.dir = ("Output/")

# read in the data

npoc.data = read.csv("Data/WHONDRS_S19S_Sediment_NPOC.csv")

npoc.data$Type = NA

for (i in 1:nrow(npoc.data)){
  if (length(grep(pattern = "Sed_Field", npoc.data$Sample_ID[i]))>0){  
    npoc.data$Type[i] = "Field"
} else{
  npoc.data$Type[i] = "Incubation"
  }
}

npoc.inc = subset(npoc.data,npoc.data$Type == "Incubation")
npoc.field = subset(npoc.data,npoc.data$Type == "Field")

names(npoc.inc)[3] = "NPOC_INC_mg_per_L"
names(npoc.field)[3] = "NPOC_Field_mg_per_L"

npoc.inc$Site = npoc.inc$Sample_ID
npoc.inc$Site = gsub("_SED_INC_ICR","",npoc.inc$Site)

npoc.field$Site = npoc.field$Sample_ID
npoc.field$Site = gsub("_Sed_Field_ICR","",npoc.field$Site)



npoc.all.big = merge(npoc.field,npoc.inc, by = "Site")

npoc.all = cbind.data.frame(npoc.all.big$Site, npoc.all.big$NPOC_Field_mg_per_L, npoc.all.big$NPOC_INC_mg_per_L)
names(npoc.all) = c("Sample","NPOC_Field","NPOC_INC")

npoc.all$Ratio_field_inc = npoc.all$NPOC_Field/npoc.all$NPOC_INC 

# Invert the ratio so it is always greater than 1
for (i in 1:nrow(npoc.all)){
  if (npoc.all$Ratio_field_inc[i]<1){
    npoc.all$Ratio_field_inc[i] = 1/npoc.all$Ratio_field_inc[i]
  }
}

npoc.all$Site = substr(npoc.all$Sample,1,9)


###a ###############################
# Make plots as function of points removed

npoc.ordered = npoc.all[order(npoc.all$Ratio_field_inc, decreasing = T),]
saved.stats = as.data.frame(matrix(NA,nrow = nrow(npoc.ordered)-2, ncol = 4))
colnames(saved.stats) = c("Pts.removed","Rsq","slope","pvalue")

for (i in 1:(nrow(npoc.ordered)-2)){
  nums = seq(1,i)
  df = npoc.ordered[-nums,]
  mod = summary(lm(log10(df$NPOC_INC)~log10(df$NPOC_Field)))
  saved.stats$Pts.removed[i] = i
  saved.stats$Rsq[i] = format(mod$r.squared, digits=4)
  saved.stats$slope[i]=format(mod$coefficients[2,1], digits=4)
  saved.stats$pvalue[i] = format(mod$coefficients[2,4], digits=4)
}

# Calculate half saturation constant and number of points to remove

half.sat = round(((as.numeric(max(saved.stats$Rsq)) - as.numeric(min(saved.stats$Rsq)))/2) + as.numeric(min(saved.stats$Rsq)),2)
  
pts.to.remove = nrow(subset(saved.stats,saved.stats$Rsq<half.sat))

# plot Rsq as a function of pts removed
ggplot(saved.stats, aes(x = Pts.removed, y = as.numeric(Rsq))) + coord_cartesian()+
  geom_point(size = 3) + 
  geom_vline(xintercept = pts.to.remove, colour = 'red')+labs(y = expression(R^{2}), x = "Number of points removed")+ theme_bw()+theme(legend.title = element_blank(),  legend.background = element_rect(fill = 'NA'), legend.text = element_text(size=16,face="bold"))+ 
  theme(axis.text.x=element_text(size=16))+theme(axis.text.x=element_text(colour = c("black","black")))+
  theme(axis.text.y=element_text(colour = c("black","black")))+
  theme(aspect.ratio=1)+
  theme(axis.text.y=element_text(size=16))+
  theme(axis.title.x =element_text(size=16))+
  theme(axis.title =element_text(size=16))+
  theme(axis.title.y =element_text(size=16, color = "black"))

ggsave(paste0(output.dir,"Figure_S1_",Sys.Date(),".png"),width = 7.38, height = 7.38)
ggsave(paste0(output.dir,"Figure_S1_",Sys.Date(),".pdf"))


write.csv(npoc.ordered,"NPOC_ordered_to_filter.csv", row.names = F)
