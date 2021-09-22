#Loading packages 
require(reshape2); require(ggplot2); require(easycsv)
library(TPD); library(ggplot2); library(reshape2);library(vegan)
library(ggthemes);library(ggpubr);library(dplyr)
library(stringr);library(e1071)  

#options(digits=10)
rm(list=ls());graphics.off()

#####################################################################################
# Set directories
basedatafolder = "//PNL/Projects/SBR_SFA/RC4/Respiration_ms/Input/"

input.data.folders = c("Processed_ICR_INC","Processed_ICR_Field")

for (i in 1:length(input.data.folders)){
  
  setwd (paste0(basedatafolder,input.data.folders[i]))
  
  # read in the data
  # read in molecular file
  data = read.csv(list.files(pattern = "Processed.*Data.csv"), row.names = 1) # Keeping data and mol-data separate to ensure they are unaltered
  mol = read.csv(list.files(pattern = "Processed.*Mol.csv"), row.names = 1)
  
  #Fixing the names from the data so we can remove the poorly calibrated samples
  #data.names= colnames(data)
  colnames(data)=gsub("ICR.U","ICR-U",colnames(data))
  colnames(data)=gsub("ICR.M","ICR-M",colnames(data))
  colnames(data)=gsub("ICR.D","ICR-D",colnames(data))
  
  
  # Loading in poorly calibrated samples and removing them if necessary
  poor.cal = read.csv(list.files(pattern = "*Poorly_Calibrated_Samples.csv"), stringsAsFactors = F)
  
  
  data = data[,-which(colnames(data) %in% poor.cal$samples)]
  
  if(length(which(rowSums(data) == 0)) > 0){
    mol = mol[-which(rowSums(data) == 0),]
    data = data[-which(rowSums(data) == 0),]
    
  }
  
  # Changing all the peaks for presence absence
  data[data > 0] = 1

  # Removing data that doesn't have a Mol Form
  data.v2 = data[!is.na(mol$MolForm),]
  mol.v2 = mol[!is.na(mol$MolForm),]

  
  
  # Calculating number of peaks present in each sample

  el.comp = unique(mol.v2$El_comp)
  chem.class = unique(mol.v2$Class)
  class = c(el.comp,chem.class)
  
  summary = as.data.frame(matrix(NA,nrow = (ncol(data)),ncol = (3+length(el.comp)+length(chem.class))))
  colnames(summary) = c("Sample_ID_ICR","Total_number_of_peaks","Total_Peaks_with_Molecular_Formulas")
  
  for (j in 1:ncol(data)){
    
    summary$Sample_ID_ICR[j] = colnames(data)[j] 
    summary$Total_number_of_peaks[j] = sum(data[,j])
    
    summary$Total_Peaks_with_Molecular_Formulas[j] = sum(data.v2[,j])
  }
  
 for (v in 1:length(chem.class)){
   for (k in 1:ncol(data)){
     df = merge(mol,data, by=0, all=TRUE)
     colnames(df)[1] = "Mass"
     df1 = subset(df,df$Class == chem.class[v])
     col = 3+v
     summary[k,col] = sum(df1[,k+38])
 }
   colnames(summary)[col] = chem.class[v] 

 }
  
  rm(v)
  rm(k)
  rm(df1)
  
  for (v in 1:length(el.comp)){
    for (k in 1:ncol(data)){
      df = merge(mol,data, by=0, all=TRUE)
      df1 = subset(df,df$El_comp == el.comp[v])
      col = 12+v
      summary[k,col] = sum(df1[,k+38])
    }
    colnames(summary)[col] = el.comp[v] 
  
  }
  
  write.csv(summary,paste0(basedatafolder,input.data.folders[i],"_Peaks.csv"), row.names = F)

}

  