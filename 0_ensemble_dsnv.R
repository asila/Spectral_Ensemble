# Script for fitting calibration models for infrared spectroscopy data. ----------------------------------------------------------------
# Author: Andrew Sila , April 2016. --------------------------------

# Begin by checking the required packages are installed. --------------------------------
is.installed <- function(anypkg){
  
  is.element(anypkg, installed.packages()[,1])
  
}

required.packages <- c("caret","prospectr","doParallel","reshape")

installp <- which(!is.installed(required.packages) == TRUE)

#Install missing packages

if (length (installp) > 0){
  
  install.packages(required.packages[installp])
}

library(caret)
library(prospectr)
library(reshape)
library(doParallel)
library(readr)

registerDoParallel()
getDoParWorkers()
ensemble <- function(wd,infrared.data,reference.data,hout,process){
  setwd(wd)
  
  si <- infrared.data
  
  ref <- reference.data
  
  #0. Raw
  mir <- si
  
  raw <- mir
  
  hout <- hout
  
  # 'NIR: Si spectra
  
  prefx <- substr(colnames(si),1,1)[90]
  
  wavenumbers <- round(as.numeric(gsub(prefx,"",colnames(si[,-1]))),1)
  
  colnames(si) <- c("SSN",wavenumbers)
  
  #preprocess sensor data
  
  # set preprocessing methods.
  
  if(process == "none"){
    
    #0. Raw
    mir <- si
    
    raw <- mir
    
    
    colnames(raw)<-c("SSN",colnames(raw[,-1]))
    
    write.table(raw,file= paste0( process, " processed spectra.csv"),sep=",",row.names=FALSE) # same as what was read, creates a duplicate!
  }
  
  
  # set preprocessing methods.
  
  if(process == "msc"){
    
    #1. MSC
    mir.msc<-cbind(as.vector(mir[,1]),msc(as.matrix(mir[,-1])))
    
    colnames(mir.msc)<-c("SSN",colnames(mir[,-1]))
    
    write.table(mir.msc,file= paste0( process, " processed spectra.csv"),sep=",",row.names=FALSE)
  }
  
  #2. SNV
  
  if(process == "snv"){
    
    mir.snv<-cbind(as.vector(mir[,1]),standardNormalVariate(as.matrix(mir[,-1])))
    
    colnames(mir.snv) <- c("SSN",colnames(mir.snv[,-1]))
    
    write.table(mir.snv,file= paste0( process, " processed spectra.csv"),sep=",",row.names=FALSE)
    
  }
  
  
  #3. SNV-detrending after snv a second order polynomial is fit on the spectra and substracted from it
  
  if(process =="detrend"){
    
    mir.de<-cbind(as.vector(mir[,1]),detrend(mir[,-1],as.numeric(substr(colnames(mir[,-1]),2,19))))
    
    colnames(mir.de)<-c("SSN",colnames(mir.de[,-1]))
    
    write.table(mir.de,file= paste0( process, " processed spectra.csv"),sep=",",row.names=FALSE)
  }
  
  #4. diff
  if(process == "diff"){
    gd1 <- cbind(as.vector(mir[,1]),t(diff(t(mir[,-1]),differences=1,lag=10)))
    
    colnames(gd1) <- c("SSN", colnames(gd1[,-1]))
    
    write.table(gd1,file= paste0( process, " processed spectra.csv"),sep=",",row.names=FALSE)
    
  }
  
  
  #5. First derivative
  
  mir1 <- as.matrix(mir[,-1])
  
  wave <- as.numeric(substr(colnames(mir1),2,19))
  
  prefx <- substr(colnames(mir1),1,1)[900]
  
  colnames(mir1) <- wave
  
  if(process == "derivative"){
    
    de1 <- trans(mir1,tr = "derivative",order = 1,gap = 23)
    
    der1 <- rev(as.data.frame(de1$trans))
    
    colnames(der1) <- paste0(prefx,wave)
    
    # Save derivative spectra.
    der1.ssn <- as.data.frame(cbind(as.vector(mir[,1]),der1))
    
    colnames(der1.ssn) <- c("SSN",colnames(der1))
    
    write.table(der1.ssn,file = paste0( process, " processed spectra.csv"),sep = ",",row.names = FALSE)
  }
  
  # Use preprocessed table
  
  der1.ssn<-as.data.frame(read_csv(paste0( process, " processed spectra.csv")))
  
  # Add SNV
  #der1.ssn <-cbind(as.vector(der1.ssn[,1]),standardNormalVariate(as.matrix(der1.ssn[,-1])))
  
  #colnames(der1.ssn) <- c("SSN",colnames(der1.ssn[,-1]))
  
  #write.table(der1.ssn,file= paste0( 'First_deriv + SNV', " processed spectra.csv"),sep=",",row.names=FALSE)
  
  #der1.ssn<-as.data.frame(read_csv(paste0('First_deriv + SNV', " processed spectra.csv")))
  
  # Merge with preprocessed spectra.
  si.r <- na.omit(merge(ref,der1.ssn))
  
  #Select a random set for calibration and validation
  if(hout!=0){
    m <- round(hout*nrow(si.r))
    
    set.seed(98752)
    intest <- sample(1:nrow(si.r),m)
    
    cal <- si.r [-intest, ]
    
    val <- si.r [intest, ]
  }
  
  if(hout==0){
    cal <- si.r 
    
    val <- si.r
  }
  
  ## check potential labels for response variable
  
  ref.hd <- colnames(ref)[-1] # Remove metadata
  
  for (q in 1 : length(ref.hd)){
    
    lt <- cal[,ref.hd[q]]
    
    lv <- val[,ref.hd[q]]
    
    
    
    # Soil spectral features
    
    mirt <- cal[,-c(1:ncol(ref))] # soil MIR
    
    mirv <- val[,-c(1:ncol(ref))] # soil MIR
    
    mir.a <- si[,-1] # all mir spectra for prediction 
    
    mir.a <- rbind(mirt,mirv)
    
    # bartMachine models ------------------------------------------------------
    
    options(java.parameters = "-Xmx100g")
    
    library(bartMachine)
    
    # Start doParallel to parallelize model fitting
    
    mc <- makeCluster(detectCores())
    
    registerDoParallel(mc)
    
    k <- round(0.002*length(lt))
    
    k <- ifelse(k < 2, 3, k)
    
    myFolds <- createFolds(lt, k=)
    
    # My control --------------------------------------------------------------
    tc <- trainControl(
      method = "boot",
      classProbs = FALSE,
      verboseIter = FALSE,
      savePredictions = TRUE,
      index = myFolds,
      adaptive = list(min = 5, alpha = 0.05, method ="gls", complete = TRUE),
      search = "random"
    )
    
    # Control setup
    
    #tc <- trainControl(method = "cv", returnResamp = "all", allowParallel = T)
    
    # Fit model
    options(java.parameters = "-Xmx100g")
    
    require(bartMachine)
    
    #mir.bar <- train(mirt, lt,
                     
                    # method = 'bartMachine', 
                     
                     #preProc = c("center", "scale"),
                     
                     #trControl = tc,
                     
                     #tuneLength = 2,
                     
                     #serialize = TRUE,
                     
                     #seed = 123)
    
    #print(bartcv)
    
    mir.bar <- bartMachine(mirt,lt, serialize = TRUE)
        
    bar_mir <- predict(mir.bar, mirv)
    
    bar_mir.a <- predict(mir.bar, mir.a) ## predict full set
    
    
    # Save the model 
    
    saveRDS(mir.bar, file = paste0(ref.hd[q],"_bart.RDS"))
    
    
    rm("mir.bar")
    
    
    gcinfo(TRUE)
    
    stopCluster(mc)
    
    detach("package:bartMachine", unload=TRUE)
    
    # RF models ---------------------------------------------------------------
    
    library(doParallel)
    
    library(randomForest)
    
    
    # Start doParallel to parallelize model fitting
    
    mc <- makeCluster(detectCores())
    
    registerDoParallel(mc)
    
    
    # Control setup
    
    set.seed(1385321)
    
    tc <- trainControl(method = "cv", allowParallel = T)
    
    
    # Tuning parameters
    
    #tg <- expand.grid(mtry=seq(10, 150, by=10))
    
    tg <- expand.grid(mtry=seq(5, 30, by=5))
    
    # Fit model
    
    mir.rfo <- train(mirt, lt,
                     
                     preProc = c("center", "scale"),
                     
                     method = "rf",
                     
                     ntree = 5,
                     
                     tuneGrid = tg,
                     
                     trControl = tc)
    
    print(mir.rfo)
    
    rfo_mir <- predict(mir.rfo, mirv) ## predict validation set
    
    ref_mir.a <- predict(mir.rfo,mir.a)
    rfo_mir.a <- predict(mir.rfo, mirv) ## predict validation set; replace with mir.a
    
    # Save the model 
    
    saveRDS(mir.rfo, file = paste0(ref.hd[q],"_RFO.RDS"))
    
    
    rm("mir.rfo")
    
    
    stopCluster(mc)
    
    #detach("package:randomForest", unload=TRUE)
    
    # GBM models --------------------------------------------------------------
    
    library(plyr)
    
    library(gbm)
    
    # Start doParallel to parallelize model fitting
    
    mc <- makeCluster(detectCores())
    
    registerDoParallel(mc)
    
    # Control setup
    
    set.seed(1385321)
    
    tc <- trainControl(method = "repeatedcv", repeats=5, allowParallel = T)
    
    
    
    # Tuning parameters
    
    tg <- expand.grid(.n.trees=seq(1, 7, by=1), 
                      
                      .interaction.depth = 2,
                      
                      .shrinkage = 0.1,
                      
                      .n.minobsinnode = 5)
    
    metric <- "RMSE"
    trainControl <- trainControl(method="cv", number=3)
    
    
    
    # Fit model
    
    #mir.gbm <- train(mirt, lt, 
    
    # method = "gbm",#,
    
    #trControl = tc,
    
    #metric = metric)#,
    
    #tuneGrid = tg)
    
    #print(mir.gbm)
    
    #gbm_mir <- predict(mir.gbm, mirv) ## predict validation set
    
    #gbm_mir.a <- predict(mir.gbm, mir.a) ## predict full set
    
    
    # Save the model 
    
    #saveRDS(mir.gbm, file = paste0(ref.hd[q],"_gbm.RDS"))
    
    
    #rm("mir.gbm")
    
    
    
    
    stopCluster(mc)
    
    detach("package:gbm", unload=TRUE)
    
    # xgbm
    # Use  extreme gradient boosting
    
    library(xgboost)
    
    # Start doParallel to parallelize model fitting
    
    mc <- makeCluster(detectCores())
    
    registerDoParallel(mc)
    
    # Control setup
    
    set.seed(1385321)
    
    maxTrees <- 0.5*nrow(mirt)
    
    shrinkage <- 0.1
    
    gamma <- 2
    
    earlyStopRound <- 8
    
    # Matrix
    dtrain <- xgb.DMatrix(as.matrix(mirt), label = lt)
    
    dtest <- xgb.DMatrix(as.matrix(mirv))
    
    dtest.a <- xgb.DMatrix(as.matrix(mir.a))
    
    xgbCv <- xgb.cv(params=list(eta = shrinkage, gamma = gamma, objective = "reg:linear"), data = dtrain, nrounds = maxTrees, eval_metric = "rmse", nfold =3) 
    
    numTrees <- min(which(xgbCv$evaluation_log$test_rmse_mean==min(xgbCv$evaluation_log$test_rmse_mean)))
    
    mir.xgb <- xgboost(params=list(eta = shrinkage, gamma = gamma, objective = "reg:linear"), data = dtrain, nrounds = numTrees, eval_metric = "rmse")
    
    print(mir.xgb)
    
    xgb_mir <- predict(mir.xgb, dtest) ## predict validation set
    
    xgb_mir.a <- predict(mir.xgb, dtest.a) ## predict full set
    
    
    # Save the model 
    
    
    #saveRDS(mir.xgb, file = paste0(ref.hd[q],"_xgb.RDS"))
    
    xgb.save(mir.xgb, paste0(ref.hd[q],"_xgb.RDS"))
    
    rm("mir.xgb")
    
    stopCluster(mc)
    
    detach("package:xgboost", unload=TRUE)
    
    # PLS models --------------------------------------------------------------
    
    library(pls)
    
    
    # Start doParallel to parallelize model fitting
    
    mc <- makeCluster(detectCores())
    
    registerDoParallel(mc)
    
    
    
    # Control setup
    
    set.seed(1385321)
    
    tc <- trainControl(method = "repeatedcv", repeats = 5, allowParallel = TRUE)
    
    # Fit models
    
    mir.pls <- train(mirt, log(ifelse(lt==0,NA,lt)),
                     
                     preProc = c("center", "scale"),
                     
                     method = "pls",
                     
                     tuneGrid = expand.grid(ncomp=seq(2, 0.2*nrow(mirt), by=1)),
                     
                     trControl = tc)
    
    print(mir.pls)
    
    pls_mir <- exp(predict(mir.pls, mirv)) ## predict validation set
    
    pls_mir.a <- exp(predict(mir.pls, mir.a)) ## predict full set
    
    
    # Save the model 
    
    saveRDS(mir.pls ,file = paste0(ref.hd[q],"_pls.RDS"))
    
    
    rm("mir.pls")
    
    
    stopCluster(mc)
    
    #detach("package:pls", unload=TRUE)
    
    
    
    # cubist models 
    library(Cubist)
    
    # Start doParallel to parallelize model fitting
    
    mc <- makeCluster(detectCores())
    
    registerDoParallel(mc)
    # Control setup
    
    set.seed(1385321)
    
    ctc <- cubistControl(unbiased = TRUE, rules = NA, sample = 0.3)
    
    
    
    # Fit models
    #committees <- round(c(0.25*nrow(mirt),0.35*nrow(mirt),0.45*nrow(mirt),0.5*nrow(mirt),0.75*nrow(mirt)),0)
    
    mir.cubist <- train(mirt, lt,
                        
                        preProc = c("center", "scale"),
                        
                        method = "cubist")
    
    #tuneGrid = expand.grid(.committees = committees,.neighbors = 0))
    
    #.committees = committees)
    
    #trControl = ctc)
    
    print(mir.cubist)
    
    cubist_mir <- predict(mir.cubist, mirv) ## predict validation set
    
    cubist_mir.a <- predict(mir.cubist, mir.a) ## predict full set
    
    # Save the model 
    
    saveRDS(mir.cubist, file = paste0(ref.hd[q],"_cubist.RDS"))
    
    rm("mir.cubist")
    
    
    stopCluster(mc)
    
    detach("package:Cubist", unload=TRUE)
    
    # Model stacking setup ----------------------------------------------------
    
    pmirv <- as.data.frame(cbind(lv, rfo_mir, xgb_mir, pls_mir,bar_mir,cubist_mir))
    
    
    pmirv.a <- as.data.frame(cbind(rfo_mir.a, xgb_mir.a, pls_mir.a,bar_mir.a,cubist_mir.a))
    
    names(pmirv) <- c("L", "RFO", "XGB", "PLS","BART","CUBIST")
    
    names(pmirv.a) <- c("RFO", "XGB", "PLS","BART","CUBIST")
    
    
    
    # Remove extraneous objects from memory -----------------------------------
    
    # rm(list=setdiff(ls(), pmirv"))
    
    
    
    
    # Model stacking ----------------------------------------------------------
    
    library(glmnet)
    
    
    
    
    # Start doParallel to parallelize model fitting
    
    mc <- makeCluster(detectCores())
    
    registerDoParallel(mc)
    
    
    
    
    # Control setup
    
    set.seed(1385321)
    
    tc <- trainControl(method = "cv", allowParallel = T)
    
    
    
    
    # MIR model stack
    
    set.seed(1385321)
    
    mir.ens <- train(L ~ ., data = pmirv,
                     
                     method = "glmnet",
                     
                     family = "gaussian",
                     
                     trControl = tc)
    
    print(mir.ens)
    
    saveRDS(mir.ens, file = paste0(ref.hd[q],".RDS"))
    
    ens_mir <- as.data.frame(predict(mir.ens, pmirv))
    
    names(ens_mir) <- c("ENS")
    
    ens_mir.a <- as.data.frame(predict(mir.ens, pmirv.a))
    
    names(ens_mir.a) <- c("ENS")
    
    pmirv <- cbind(pmirv, ens_mir)
    
    pmirv.a <- cbind(pmirv.a, ens_mir.a)
    
    
    
    stopCluster(mc)
    
    
    
    
    # Write data files --------------------------------------------------------
    
    write.csv(pmirv, paste0(ref.hd[q],"_pmirv.csv"), row.names=F)
    
    write.csv(pmirv.a, paste0(ref.hd[q],"_pmirv.a.csv"), row.names=F)
    
    png(file = paste0(ref.hd[q],".png"), height = 400, width = 600)
    # Prediction plots --------------------------------------------------------
    
    par(mfrow=c(2,3), mar=c(5,4.5,1,1))
    
    
    
    
    # MIR predictions # note that x & y axis limits will need to be adjusted
    
    lmin <- 0
    
    lmax <- max(pmirv$L)
    
    plot(L ~ RFO, pmirv, xlim=c(min(pmirv$L,pmirv$RFO), max(pmirv$L,pmirv$RFO)), ylim=c(min(pmirv$L,pmirv$RFO), max(pmirv$L,pmirv$RFO)), xlab = "RFO prediction", ylab = "Observed", cex.lab=1.3)
    
    abline(c(0,1), col="red")
    
    plot(L ~ XGB, pmirv, xlim=c(min(pmirv$L,pmirv$XGB), max(pmirv$L,pmirv$XGB)), ylim=c(min(pmirv$L,pmirv$XGB), max(pmirv$L,pmirv$XGB)), xlab = "XGB prediction", ylab = "Observed", cex.lab=1.3)
    
    abline(c(0,1), col="red")
    
    plot(L ~ PLS, pmirv, xlim=c(min(pmirv$L,pmirv$PLS), max(pmirv$L,pmirv$PLS)), ylim=c(min(pmirv$L,pmirv$PLS), max(pmirv$L,pmirv$PLS)), xlab = "PLS prediction", ylab = "Observed", cex.lab=1.3)
    
    abline(c(0,1), col="red")
    
    plot(L ~ CUBIST, pmirv, xlim=c(min(pmirv$L,pmirv$CUBIST), max(pmirv$L,pmirv$CUBIST)), ylim=c(min(pmirv$L,pmirv$CUBIST), max(pmirv$L,pmirv$CUBIST)), xlab = "CUBIST prediction", ylab = "Observed", cex.lab=1.3)
    
    abline(c(0,1), col="red")
    
    
    # Ensemble predictions 
    
    plot(L ~ ENS, pmirv, xlim=c(min(pmirv$L,pmirv$ENS), max(pmirv$L,pmirv$ENS)), ylim=c(min(pmirv$L,pmirv$ENS), max(pmirv$L,pmirv$ENS)), xlab = "Model ensemble prediction", ylab = "Observed", cex.lab=1.3,col = "blue",pch = 16)
    
    abline(c(0,1), col="red")
    dev.off()
  }
}



