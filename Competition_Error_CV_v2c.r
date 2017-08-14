###Competition using K matrix - inter-clonal and inter-plot - demonstration ##
###Ani A. Elias #### April, 2017 ####

library(rrBLUP) # to calculate K matrix
library(nlme) # to fit the first base model
library(LMERConvenienceFunctions) # remove outliers 
library(regress) #to fit models 

#this functon assumes that field coordinates are used for spatial analysis. It is required to change the 
  #method of distance calculation etc. if geocoordinates are used. 
#data is expected to have no rows with all NA values 
#genotype is the name of column containing genotypic predictor variable
#trait is a vector of names of traits that needs to be tested - response variable
#K is the relationship matrix from 'rrBLUP' OR snp is the snp data 
  #snp <- load(".../V6_IITA_BeagleDosages_V100715_MAF01.Rdata")
#covariates is a vector of covariates to be used as predictor variables. eg. c("a","b","a+b")
#assuming that range is on y axis and column on x axix as field coordinates; so for AR1, Range is the row and column is the column
#plotSize is a vector of plotdimensions, plotSize[1] is the length of plot in the range direction, 
         #plotSize[2] is width of plot in the column direction

#remove outliers (2.5 times the sd) before proceeding to analysis
#get an influence matrix (from distance matrix and Z matrix for genotype)
#check for competition using this influence matrix and K matrix


geno_compete_snp_compE <- function(data1, trait, genotype, plotSize, loc,s, mkfold){
for(i in c(1:length(trait))){  

  test.geno <- vector("list",mkfold)

  base.prmse <- base.pcor <- base.fail <- model1.prmse <- model1.pcor <- model1.fail <- matrix(NA,mkfold,s)
  model2.prmse <- model2.pcor <- model2.fail <- model3.prmse <- model3.pcor <- model3.fail <- matrix(NA,mkfold,s)
 
   
  data1a <- na.omit(subset(data1,select = c("Range","Column",genotype,trait[i])))
  
  fixed.effect <- formula(paste(trait[i],"1", sep="~"))
  
  #K matrix aligning with data1a
  #subsetting K matrix  - helps to remove clones that are not genotyped but phenotyped 
  table(data1a[,genotype] %in% rownames(snps)) 
  snp2.names <- intersect(data1a[,genotype], rownames(snps))
  snps2 <- subset(snps, rownames(snps) %in% snp2.names)

  #calculate K matrix
  K <- A.mat(snps2-1)
     

  data2 <- data1a


  #computing the distance matrix
  rowVec<- data2$Range
  colVec <- data2$Column

  rowVec2 <- rowVec * plotSize[1] 
  rowDistMat2 <- sapply(rowVec2, function(rowNum) abs(rowNum - rowVec2))

  colVec2 <- colVec * plotSize[2]
  colDistMat2 <- sapply(colVec2, function(colNum) abs(colNum - colVec2))

  distMat <- sqrt(rowDistMat2^2 + colDistMat2^2)

  data2$Slno <- c(1:nrow(data2))
  rownames(distMat) <- colnames(distMat) <- data2$Slno


#Z.base - to facilitate competition
  id <- factor(as.character(data2[,genotype]), levels = rownames(K))
  Z.base <- model.matrix(~id - 1)


#computing inverse of distance matrix for the competition matrix
  distMat.inv <- 1/distMat
  diag(distMat.inv) <- 0 
  rownames(distMat.inv) <- colnames(distMat.inv) <- data2$Slno
  
#k-fold cross-validation  
  #splitting data - same genotypes not present in training and test data
  data.geno <- unique(data2[,genotype])

  for(m in c(1:mkfold)){
  idx <- sample(rep(1:s, length.out=length(data.geno)))
  test.geno[[m]] <- split(data.geno, idx) 
  } 

# CV for base model
for(m in c(1:length(test.geno))){
for (s in c(1:length(test.geno[[m]]))){ 
  test.data <- data2[(data2[,genotype] %in% test.geno[[m]][[s]]),]
  train.data <-data2[!(data2[,genotype] %in% test.geno[[m]][[s]]),]

  y <- train.data[,trait[i]]

  #X design matrix 
  X.train <- model.matrix(fixed.effect,data=train.data)
  X.test <- model.matrix(fixed.effect,data=test.data)

  #Identity matrix for training data
  Identity <- diag(nrow(train.data))

  idg <- factor(as.character(train.data[,genotype]), levels = rownames(K))
    Z.geno <- model.matrix(~idg - 1)

  idgt <- factor(as.character(test.data[,genotype]), levels = rownames(K))
    Ztest.geno <- model.matrix(~idgt - 1)

  G.geno <- Z.geno%*%K%*%t(Z.geno) 

#base model
  base <- try(regress(y ~ X.train, ~G.geno, pos= rep(TRUE,2), tol = 1e-4, data = train.data),silent = TRUE)

  if(class(base) != "try-error"){

  R.residual <- Identity * base$sigma[[2]]

  Khat <- K * base$sigma[[1]] 
  Vhat <- G.geno * base$sigma[[1]] + R.residual

  gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - base$fitted)

  gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% base$beta)
    yhat.test <- gamma.test[,1] + gamma.test[,2]


base.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
base.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
}else{
  base.fail[m,s] <- 1
 }#end of try
} # end of one CV
} # end of multiple CV 

#for comparison and plotting
base.prmse.1 <- mean(base.prmse, na.rm = TRUE)
base.pcor.1 <- mean(base.pcor, na.rm = TRUE)
base.fail.1 <- sum(base.fail, na.rm = TRUE)
 

#model1, 2, and 3 demonstration
#model 1 here is Base + inter-genotypic competition 
#model 2 here is Base + inter-genotypic + competition error
#model 3 here is Base + competition error

  for(m in c(1:length(test.geno))){
    
      # to demonstration inter-clonal competition
      distMat.2 <- distMat
      distMat.inv.2 <- 1/distMat.2
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.1 <- 1/(0.009 + 0.001^distMat.inv.2)
      a <- 1/plotSize[2]
      distMat.1[distMat.1 < (1/(0.009 + 0.001^a))] <- 0
      distMat.1[is.na(distMat.1)] <- 0
      diag(distMat.1) <- 0

      Z.comp <- distMat.1 %*% Z.base


      #competition error
      distMat.4 <- distMat
      distMat.4[distMat.4 > plotSize[1]] <- 0
      distMat.4[distMat.4 == plotSize[1]] <- 1
      distMat.4[distMat.4 > 1] <- 0
    
      #competition error
      Identityall <- diag(nrow(data2))
      #distMat.4 %*% Identiti matrix is distMat.4; so ignoring that line   
      G.compE <- distMat.4 %*% Identityall %*% t(distMat.4) #can ignore multiplying with Identity but just ZZ'

      for(s in c(1:length(test.geno[[m]]))){        
        test.data <- data2[(data2[,genotype] %in% test.geno[[m]][[s]]),]
        train.data <-data2[!(data2[,genotype] %in% test.geno[[m]][[s]]),]
        
        y <- train.data[,trait[i]]

        #X design matrix 
        X.train <- model.matrix(fixed.effect,data=train.data)
        X.test <- model.matrix(fixed.effect,data=test.data)
        
        #Identity matrix for training data
        Identity <- diag(nrow(train.data))

        #genotypic design matrix
        idg <- factor(as.character(train.data[,genotype]), levels = rownames(K))
          Z.geno <- model.matrix(~idg - 1)

        idgt <- factor(as.character(test.data[,genotype]), levels = rownames(K))
          Ztest.geno <- model.matrix(~idgt - 1)

        G.geno <- Z.geno%*%K%*%t(Z.geno) 

        Z.comp.train <- Z.comp[rownames(Z.comp) %in% train.data$Slno,]
        Z.comp.test <- Z.comp[!(rownames(Z.comp) %in% train.data$Slno),]

        G.comp <- Z.comp.train%*%K%*%t(Z.comp.train) 

        Z.compE.train <- distMat.4[rownames(distMat.4) %in% train.data$Slno,]
        Z.compE.test <- distMat.4[!(rownames(distMat.4) %in% train.data$Slno),]

        G.compE <- Z.compE.train%*%Identityall%*%t(Z.compE.train) 

        model1 <- try(regress(y ~ X.train, ~G.geno + G.comp, pos= rep(TRUE,3), tol = 1e-4, data = train.data),silent = TRUE)

        if(class(model1) != "try-error"){
          R.residual <- Identity * model1$sigma[[3]]

          Khat.g <- K * model1$sigma[[1]]
          Khat.c <- K * model1$sigma[[2]]
          Vhat <- G.geno * model1$sigma[[1]] + G.comp * model1$sigma[[2]] + R.residual

          gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1$fitted)
          zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model1$fitted)

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1$beta)
          zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
            yhat.test <- zeta.test[,1] + zeta.test[,2] + zeta.test[,3] 

            model1.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))  
            model1.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
          }else{
            model1.fail[m,s] <- 1
          } # end of try

        model2 <- try(regress(y ~ X.train, ~G.geno + G.comp + G.compE, pos= rep(TRUE,4), tol = 1e-4, data = train.data),silent = TRUE)

        if(class(model2) != "try-error"){
          R.residual <- Identity * model2$sigma[[4]]
          Rc.error <- Identityall * model2$sigma[[3]]

          Khat.g <- K * model2$sigma[[1]]
          Khat.c <- K * model2$sigma[[2]]

          Vhat <- G.geno * model2$sigma[[1]] + G.comp * model2$sigma[[2]] + G.compE * model2$sigma[[3]] + R.residual

          gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model2$fitted)
          zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model2$fitted)
          delta <- Rc.error %*% t(Z.compE.train) %*% solve(Vhat) %*% (y - model2$fitted)

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model2$beta)
          zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
          delta.test <- cbind(Z.compE.test %*% delta, zeta.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3] + delta.test[,4] 

            model2.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))  
            model2.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
          }else{
            model2.fail[m,s] <- 1
          } # end of try   

        model3 <- try(regress(y ~ X.train, ~G.geno +  G.compE, pos= rep(TRUE,3), tol = 1e-4, data = train.data),silent = TRUE)

        if(class(model3) != "try-error"){
          R.residual <- Identity * model3$sigma[[3]]
          Rc.error <- Identityall * model3$sigma[[2]]

          Khat.g <- K * model3$sigma[[1]]
        
          Vhat <- G.geno * model3$sigma[[1]] + G.compE * model3$sigma[[2]] + R.residual

          gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model3$fitted)
          delta <- Rc.error %*% t(Z.compE.train) %*% solve(Vhat) %*% (y - model3$fitted)

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model3$beta)
          delta.test <- cbind(Z.compE.test %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3] 

            model3.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))  
            model3.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
          }else{
            model3.fail[m,s] <- 1
          } # end of try     
       } #end of one CV       
    } # end of multiple CV 

    model1.prmse.1 <- mean(model1.prmse, na.rm = TRUE)
    model1.pcor.1 <- mean(model1.pcor, na.rm = TRUE)
    model1.fail.1 <- sum(model1.fail, na.rm = TRUE)

    model2.prmse.1 <- mean(model2.prmse, na.rm = TRUE)
    model2.pcor.1 <- mean(model2.pcor, na.rm = TRUE)
    model2.fail.1 <- sum(model2.fail, na.rm = TRUE)

    model3.prmse.1 <- mean(model3.prmse, na.rm = TRUE)
    model3.pcor.1 <- mean(model3.pcor, na.rm = TRUE)
    model3.fail.1 <- sum(model3.fail, na.rm = TRUE)

print("model end")  
   

#best model from each set
  models.prmse <- cbind(base.prmse.1,model1.prmse.1,model2.prmse.1, model3.prmse.1)
  models.pcor <- cbind(base.pcor.1, model1.pcor.1, model2.pcor.1, model3.pcor.1)
  models.fail <- cbind(base.fail.1, model1.fail.1, model2.fail.1, model3.fail.1)

#patch work
  if (nrow(models.prmse) >1){
    models.prmse <- models.prmse[1,]    
  }

  if (nrow(models.pcor) >1){
    models.pcor <- models.pcor[1,]    
  }

  if (nrow(models.fail) >1){
    models.fail <- models.fail[1,]    
  }


  model.win.no <- which(models.prmse == min(models.prmse))[1]


#output tables
  write.csv(models.prmse, paste0(trait[i],"-",loc,"_K_models_prmse_comp_err.csv")) 
  write.csv(models.pcor, paste0(trait[i],"-",loc,"_K_models_pcor_comp_err.csv")) 
  write.csv(models.fail, paste0(trait[i],"-",loc,"_K_models_fail_comp_err.csv")) 

#correlation and design matrices for final model

  y <- data2[,trait[i]]

  #X design matrix 
  X.mat <- model.matrix(fixed.effect,data=data2)
 
  #genotypic design matrix
  idg <- factor(as.character(data2[,genotype]), levels = rownames(K))
    Z.geno <- model.matrix(~idg - 1)
    geno <- Z.geno%*%K%*%t(Z.geno) 

  distMat.4 <- distMat
    distMat.4[distMat.4 > plotSize[1]] <- 0
    distMat.4[distMat.4 == plotSize[1]] <- 1

  #competition error
  Identity <- diag(nrow(data2))
  #distMat.4 %*% Identiti matrix is distMat.4; so ignoring that line   
  G.compE <- distMat.4 %*% Identity %*% t(distMat.4) #can ignore multiplying with Identity but just ZZ'

     
  distMat.7 <- distMat
      distMat.inv.2 <- 1/distMat.7
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.7a <- 1/(0.9 + 0.1^distMat.inv.2)
      a <- 1/plotSize[2]
      distMat.7a[distMat.7a < (1/(0.9 + 0.1^a))] <- 0
      distMat.7a[is.na(distMat.7a)] <- 0
      diag(distMat.7a) <- 0       
      
      Z.comp.7 <- distMat.7a %*% Z.base
      G.comp.7 <- Z.comp.7%*%K%*%t(Z.comp.7)    
         
         
#choosing the final model

 model <- switch(model.win.no,
    "1" <- regress(y ~ X.mat, ~geno, pos= rep(TRUE,2), tol = 1e-4, data = data2),
    "2" <- regress(y ~ X.mat, ~geno + G.comp.7, pos= rep(TRUE,3), tol = 1e-4, data = data2),
    "3" <- regress(y ~ X.mat, ~geno + G.comp.7 + G.compE, pos= rep(TRUE,4), tol = 1e-4, data = data2),
    "4" <- regress(y ~ X.mat, ~geno + G.compE, pos= rep(TRUE,3), tol = 1e-4, data = data2),       
    regress(y ~ X.mat, ~geno, pos= rep(TRUE,2), tol = 1e-4, data = data2)
  ) 

  base <- regress(y ~ X.mat, ~geno, pos= rep(TRUE,2), tol = 1e-4, data = data2)

#summary of models
  base.formula <- capture.output(base$formula)
  base.Vformula <- capture.output(base$Vformula)    
  base.out <- capture.output(summary(base))

  model.formula <- capture.output(model$formula)
  model.Vformula <- capture.output(model$Vformula)
  model.out <- capture.output(summary(model))  


#output summary of final model
 cat(base.formula, file = paste0(trait[i],"-",loc,"_K_base_summary_comp_err.txt"), sep = "\n", append=TRUE)
 cat(base.Vformula, file = paste0(trait[i],"-",loc,"_K_base_summary_comp_err.txt"), sep = "\n", append=TRUE )
 cat(base.out, file= paste0(trait[i],"-",loc,"_K_base_summary_comp_err.txt"), sep="\n", append=TRUE)

 cat(model.formula, file = paste0(trait[i],"-",loc,"_K_model_summary_comp_err.txt"), sep = "\n", append=TRUE)
 cat(model.Vformula, file = paste0(trait[i],"-",loc,"_K_model_summary_comp_err.txt"), sep = "\n", append=TRUE )
 cat(model.out, file= paste0(trait[i],"-",loc,"_K_model_summary_comp.txt"), sep="\n", append=TRUE)

} # end of for w.r.t. trait  
} # end of function


  
