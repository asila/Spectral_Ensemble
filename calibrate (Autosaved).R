library(prospectr)

spectra <- read.csv('~/India_final_report/Data/Alpha_znse_spectra.csv')

ref <- read.csv('~/India_final_report/Data/Compiled_reference.csv')

b <- which(ref$OC==0)

ref$OC[b] = NA

b <- which(ref$pH>14)

ref$pH[b] = NA

 # Do first derivativative then SNV

raw0 <- savitzkyGolay(spectra[,-1],m = 1, p = 2, w = 23)

rawsnv <- standardNormalVariate(raw0)

SSN <- as.vector(spectra[,1])

spectra <- cbind(SSN,rawsnv)

# Save preprocessed spectra
write.csv(spectra,file  = "~/India_final_report/Data/dsnv_preprocessed.csv", row.names = FALSE)

spectra <- read.csv('~/India_final_report/Data/dsnv_preprocessed.csv')

# Provide ensemble script which requires spectra to be preprocessed using first derivative and snv
#source('~/India_final_report/Ensemble_Models/0_ensemble_dsnv.R', chdir = TRUE)

source('~/Code/0_ensemble_dsnv.R', chdir = TRUE)

wd <- "~/India_final_report/Models/ZnSe"

j <- which(spectra[,"SSN"] %in% ref$SSN)

# Specify hout-out set for testing the models

hout <- 0.3

ensemble(wd = wd ,infrared.data = spectra,reference.data = ref, hout = hout , process = "none")

# Full_models

wd <- "~/India_final_report/Models/ZnSe/Full_models"

j <- which(spectra[,"SSN"] %in% ref$SSN)

# Specify hout-out set for testing the models

hout <- 0

ensemble(wd = wd ,infrared.data = spectra,reference.data = ref, hout = hout , process = "none")
