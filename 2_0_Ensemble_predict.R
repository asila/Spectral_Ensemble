library(caret)

#options(java.parameters = "-Xmx5g")

library(bartMachine)

#library(soil.spec)

library(dplyr)

library(xgboost)

library(prospectr)

library(reshape2)

# Set directory with ensemble models

setwd("~/Studies/AfSIS_models")

options(java.parameters = "-Xmx30000m")

# Read new spectra from invenio S
 raw <- read.csv('~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/To_share/NRT/Averaged_raw_spectra.csv')
 
 
source("/Users/andrewsila/Dropbox/Diby/Code/make.compatible.R")

# Raw table to be made compatible to match model spectra (AfSIS + NRT)

raw.c <- make.compatible(spectra, raw, 2,2)

write.csv(raw.c , file = '~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/To_share/NRT/made_compatible.csv', row.names = FALSE)

raw <- read.csv('~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/To_share/NRT/made_compatible.csv')

raw0 <- savitzkyGolay(raw[,-1],m = 1, p = 2, w = 23)

rawsnv <- standardNormalVariate(raw0)


SSN <- as.vector(raw.c[,1])

raw.c <- cbind(SSN,rawsnv)

# Store the processed table and read it back for predictions

write.csv(raw.c, file = '~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/To_share/NRT/dsnv_made_compatible.csv', row.names = FALSE)

deriv <- read.csv('~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/To_share/NRT/dsnv_made_compatible.csv')

SSN <- as.vector(deriv[,1])

deriv <- deriv[,-1]

# Always ensure SSN is not part of the prediction object passed to the models because BART is very particular it requires only those features used in model development are found in the prediction set.

# To determine modeling features in a model e.g bart model use: md1 <- readRDS('SOC.RDS'; md1$training_data_features)

mods <-list.files("~/ICRAF SPECTRAL DATA Dropbox/Handover_Data/Prediction_ Models/HTS-xt_Tensor_27/AfSIS+NRT", pattern = "*RDS")#~/Dropbox/TNC/Results

 # Get the models
models <-  gsub(".RDS", "", mods)

pats <- c("_cubist","_pls","_RFO","_xgb","_bart")#xgb

soil <- gsub(pats[1], "",models)

soil <- gsub(pats[2], "",soil)

soil <- gsub(pats[3], "",soil)

soil <- gsub(pats[4], "",soil)

soil <- gsub(pats[5], "",soil)

soil <- gsub(".RDS", "",soil)

soil <- unique(soil)

# Exclude TC

#b <- which(soil == 'TC')

#soil <- soil[-b]

# Identify the extreme gradient boosted regression tree models.

j <- grep("_xgb", pats)

# update the model list to read

pats <- pats[-j]

pred.all <- NULL

for( s in 1:length(soil)){ 	  
  predc <- NULL
  
  for(p in 1:length(pats)){
    

    predc <- cbind(predc,predict(readRDS(paste0(soil[s],pats[p],".RDS")),deriv))
        
  }
  
  colnames(predc) <- gsub("_","",toupper(pats))
  
  # sapply the xgb model
  predc <- cbind(predc,predict(xgb.load(paste0(soil[s],"_xgb.RDS")),as.matrix(deriv)))
  
  colnames(predc) <- c(gsub("_","",toupper(pats)),"XGB")
  
  p <- which(colnames(predc)=="PLS")
  
  predc[,p] <- exp(predc[,p])
  
  # Read full model
  ens <- readRDS(paste0(soil[s],".RDS"))
  
  pred.a <- round(predict(ens, predc),3)
  
  pred.all <- cbind(pred.all,pred.a)
  
}


predcssn <- cbind(SSN, pred.all)

colnames(predcssn) <- c("SSN",soil)

write.csv(as.data.frame(predcssn), "~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/To_share/NRT/To_report/MIR_predicted CN.csv", row.names = FALSE)

