library(caret)

suppressMessages(library(gridExtra))

suppressMessages(library(yardstick))

library(ggplot2)


setwd("~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/Bertin/Gambia/Data/holdout_models/FUll_models")

source('~/Dropbox/Scripts/multiplot.R', chdir = TRUE)
# read cross-validated predictions with reference data

p.ls <- list.files(pattern = "_pmirv.csv")

params <- gsub("_pmirv.csv","",p.ls)

mss <- NULL

ensma <- NULL

for (k in 1:length(p.ls)){
	
	pds <- read.csv(p.ls[k])
	
	colnames(pds) <- c(params[k],colnames(pds[,-1]))

	msummary <- NULL
	
	# Get model summary and save
	for (p in 2:ncol(pds)){
		
		hout <- na.omit (cbind(pds[,c(1,p)]))
		
		testing.parameters <- round(postResample(hout[,1],hout[,2]),2)
		
		RSQ <- testing.parameters[2]
  	
  		RMSEP <- testing.parameters[1]
  	
  		model.summary <- c(colnames(pds)[1],colnames(pds)[p],RSQ,RMSEP)
		
		msummary <- rbind(msummary, model.summary)
		
	}

	# Get plots
	
	hout <- na.omit (cbind(pds[,c(1,2)]))
	
	testing.parameters <- round(postResample(hout[,1],hout[,2]),2)
		
	RSQ <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
	pm <- as.data.frame(hout)
  	
  	p <- ggplot(pm, aes(x = pm[,2],y = pm[,1]))+
	
	geom_point(col = "blue",size = 3,alpha = 0.8)+
	
	ggtitle(colnames(pm)[1]) +
	
	ylab("Measured")+
	
	xlim(range(pm))+
	
	ylim(range(pm))+
	
	xlab(colnames(pm)[2])
	
	p <- p+stat_smooth(method = lm, se=FALSE, color='red',alpha = 0.5)
	
	p <- p+theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
	
	p <- p + theme(text = element_text(size = 20)) # this will change all text size 
	
	p <- p + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 1.5, vjust = -2.8)
	
	p <- p + annotate('text', label=paste('RMSEP==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.3, vjust = -1.4)
	
	p <- p + theme(plot.title = element_text(hjust =0.5))
	
	# GBM plot
	hout <- na.omit (cbind(pds[,c(1,3)]))
	
	testing.parameters <- round(postResample(hout[,1],hout[,2]),2)
		
	RSQ <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
	pm <- as.data.frame(hout)
  	
  	p.gbm <- ggplot(pm, aes(x = pm[,2],y = pm[,1]))+
	
	geom_point(col = "blue",size = 3,alpha = 0.8)+
	
	ggtitle(colnames(pm)[1]) +
	
	ylab("Measured")+
	
	xlim(range(pm))+
	
	ylim(range(pm))+
	
	xlab(colnames(pm)[2])
	
	p.gbm <- p.gbm+stat_smooth(method = lm, se=FALSE, color='red',alpha = 0.5)
	
	p.gbm <- p.gbm+theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
	
	p.gbm <- p.gbm + theme(text = element_text(size = 20)) # this will change all text size 
	
	p.gbm <- p.gbm + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 1.5, vjust = -2.8)
	
	p.gbm <- p.gbm + annotate('text', label=paste('RMSEP==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.3, vjust = -1.4)
	
	p.gbm <- p.gbm + theme(plot.title = element_text(hjust =0.5))
	
	# PLS plot
	hout <- na.omit (cbind(pds[,c(1,4)]))
	
	testing.parameters <- round(postResample(hout[,1],hout[,2]),2)
		
	RSQ <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
	pm <- as.data.frame(hout)
  	
  	p.pls <- ggplot(pm, aes(x = pm[,2],y = pm[,1]))+
	
	geom_point(col = "blue",size = 3,alpha = 0.8)+
	
	ggtitle(colnames(pm)[1]) +
	
	ylab("Measured")+
	
	xlim(range(pm))+
	
	ylim(range(pm))+
	
	xlab(colnames(pm)[2])
	
	p.pls <- p.pls+stat_smooth(method = lm, se=FALSE, color='red',alpha = 0.5)
	
	p.pls <- p.pls+theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
	
	p.pls <- p.pls + theme(text = element_text(size = 20)) # this will change all text size 
	
	p.pls <- p.pls + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 1.5, vjust = -2.8)
	
	p.pls <- p.pls + annotate('text', label=paste('RMSEP==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.3, vjust = -1.4)
	
	p.pls <- p.pls + theme(plot.title = element_text(hjust =0.5))
	
	# bart
	hout <- na.omit (cbind(pds[,c(1,5)]))
	
	testing.parameters <- round(postResample(hout[,1],hout[,2]),2)
		
	RSQ <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
	
	pm <- as.data.frame(hout)
  	
  	p.bart <- ggplot(pm, aes(x = pm[,2],y = pm[,1]))+
	
	geom_point(col = "blue",size = 3,alpha = 0.8)+
	
	ggtitle(colnames(pm)[1]) +
	
	ylab("Measured")+
	
	xlim(range(pm))+
	
	ylim(range(pm))+
	
	xlab(colnames(pm)[2])
	
	p.bart <- p.bart+stat_smooth(method = lm, se=FALSE, color='red',alpha = 0.5)
	
	p.bart <- p.bart+theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
	
	p.bart <- p.bart + theme(text = element_text(size = 20)) # this will change all text size 
	
	p.bart <- p.bart + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 1.5, vjust = -2.8)
	
	p.bart <- p.bart + annotate('text', label=paste('RMSEP==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.3, vjust = -1.4)
	
	p.bart <- p.bart + theme(plot.title = element_text(hjust =0.5))
	
	
		# Cubist
	hout <- na.omit (cbind(pds[,c(1,6)]))
	
	testing.parameters <- round(postResample(hout[,1],hout[,2]),2)
		
	RSQ <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
	
	pm <- as.data.frame(hout)
  	
  	p.cubist <- ggplot(pm, aes(x = pm[,2],y = pm[,1]))+
	
	geom_point(col = "blue",size = 3,alpha = 0.8)+
	
	ggtitle(colnames(pm)[1]) +
	
	ylab("Measured")+
	
	xlim(range(pm))+
	
	ylim(range(pm))+
	
	xlab(colnames(pm)[2])
	
	p.cubist <- p.cubist+stat_smooth(method = lm, se=FALSE, color='red',alpha = 0.5)
	
	p.cubist <- p.cubist+theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
	
	p.cubist <- p.cubist + theme(text = element_text(size = 20)) # this will change all text size 
	
	p.cubist <- p.cubist + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 1.5, vjust = -2.8)
	
	p.cubist <- p.cubist + annotate('text', label=paste('RMSEP==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.3, vjust = -1.4)
	
	p.cubist <- p.cubist + theme(plot.title = element_text(hjust =0.5))
	
	#Ensemble
	hout <- na.omit (cbind(pds[,c(1,7)]))
	
	soilp <- colnames(hout)[1]
	
	colnames(hout) <- c('measured','predicted')
	
	testing.parameters <- round(postResample(hout[,1],hout[,2]),2)
		
	RSQ <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
  	RPD <- round(as.numeric(rpd(hout, measured, predicted)[3]),2)
  	
  	m_set <- metric_set(rsq,rmse,rpd,rpiq)
	
	mmmm <- m_set(hout, truth = measured, estimate = predicted)
	
	mmmm <- t(mmmm)[3,]
		
	ensm <- c(params[k],round(as.numeric(mmmm),2))
  	
  	names(ensm) <- c('Property','R_squared','RMSEP',"RPD","RPIQ")
  	
  	#ensm <- c(soilp,RSQ,RMSEP,RPD)
  	
  	#names(ensm) <- c('Property','R_squared','RMSEP',"RPD")
	
	pm <- as.data.frame(hout)
  	
  	p.ens <- ggplot(pm, aes(x = pm[,2],y = pm[,1]))+
	
	geom_point(col = "purple",size = 3,alpha = 0.3)+
	
	ggtitle(colnames(pds)[1]) +
	
	ylab("Measured")+
	
	xlab("Predicted")+

	
	xlim(range(pm))+
	
	ylim(range(pm))
	
	#xlab(colnames(pm)[2])
	
	p.ens <- p.ens+stat_smooth(method = lm, se=FALSE, color='purple',alpha = 0.1)#red
	
	p.ens <- p.ens+theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
	
	p.ens <- p.ens + theme(text = element_text(size = 20)) # this will change all text size 
	
	p.ens <- p.ens + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 1.5, vjust = -2.8)
	
	p.ens <- p.ens + annotate('text', label=paste('RMSEP==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.3, vjust = -1.4)
	
	p.ens <- p.ens + theme(plot.title = element_text(hjust =0.5))
	
   png(file = paste0("~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/Bertin/Gambia/Data/Figures/",params[k],".png"), height = 600, width = 900)

	multiplot(p,p.gbm, p.pls,p.cubist,p.bart,p.ens, cols=3)
	
	dev.off()
	
	#Save ens only
p.ens.o <- p.ens + xlab("Predicted values") + theme(text = element_text(size = 12))

ggsave(file = paste0("~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/Bertin/Gambia/Data/Figures/EnSe_only/",params[k],".png"),p.ens.o, height = 4, width = 4)

mss <- rbind(mss, msummary)

ensma <- rbind(ensma,ensm)

	}
	
write.csv(mss,  file = "~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/Bertin/Gambia/Data/Figures/Ensemble models summary.csv", row.names = FALSE)

write.csv(ensma,  file = "~/ICRAF SPECTRAL DATA Dropbox/Andrew Sila/Bertin/Gambia/Data/Figures/Stacked model summary.csv", row.names = FALSE)

