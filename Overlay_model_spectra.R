library(soil.spec)
library(prospectr)

# Read Alpha spectra
alp <- read.csv('/Volumes/Data/AfSIS_models/logins/Raw_spectra alpha 2.csv')

# Read Model data; if not available extract this from the model

calib <- read.csv('/Volumes/Data/AfSIS_models/none processed spectra.csv')

############. Alpha 2
alp2 <- alp[,-c(2:17)]

 # Do first derivativative then SNV; if first derivative was done for the calibration only then skip SNV
wave <- as.numeric(substr(colnames(alp2[,-1]),2,19))

colnames(alp2) <- c("SSN",wave)

raw0 <- as.matrix(alp2[,-1])

colnames (raw0) <- wave

#Savitzky_Golay

#der <- as.data.frame(trans(raw0)$trans)

de1 <- trans(raw0,tr = "derivative",order = 1,gap = 23)
  
der <- rev(as.data.frame(de1$trans))
  
colnames(der) <- paste0("X",colnames(der))

derssn <- cbind(as.vector(alp2[,1]),standardNormalVariate(der))

colnames(derssn) <- c('SSN', colnames(der))

derssn <- rbind(calib,derssn)

pcs <- prcomp(derssn[,-1])

pcss <- pcs$x[,1:10]

pcss[1:6,]

plot(pcss)

points(pcss[1:2434,1:2], col = "red")

points(pcss[-c(1:2434),1:2], col = "blue")

var <- round(summary(pcs)$importance[2,] * 100, 1)

scores <- cbind("Calib",as.data.frame(pcs$x[,1:5])) # get first 5 principal components

names(scores) <- c("set", colnames(scores[,-1]))

scores <- as.matrix(scores)

scores[-c(1:2434),1] <- "new spectra"

scores <- as.data.frame(scores)

write.csv(scores, file = "/Volumes/Data/AfSIS_models/logins/alp2 Calib and Pred scores.csv", row.names = FALSE)


scores <- read.csv("/Volumes/Data/AfSIS_models/logins/alp2 Calib and Pred scores.csv")

#sp <- sp +  labs(color = "set")
sp <- ggplot(scores, aes(x = PC1, y =PC2, colour = set)) +

geom_point(size = 1.5, alpha = 0.85 ) +
 
ggtitle("Alpha2 and Calibration PCA scores plot") +
  
      
    xlab(paste0("PC1 explains ", var[1], "% total variance")) +
  
    ylab(paste0("PC2 explains ", var[2], "% total variance")) +

    theme_bw() +

    theme(
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
    )
sp <- sp + theme(plot.title = element_text(hjust = 0.5))

sp <- sp + scale_color_manual(values =c("orange","grey"))

ggsave(filename  = "/Volumes/Data/AfSIS_models/logins/Calibration and pred scores.png", height = 6, width = 6,sp)

# Save preprocessed alp2nios spectra
write.csv(derssn[-c(1:2434),],file  = "/Volumes/Data/AfSIS_models/logins/pred_preprocessed.csv", row.names = FALSE)











