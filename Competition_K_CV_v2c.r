###Inter-clonal competition using K matrix ##
###Ani A. Elias #### March, 2016 ####

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


geno_compete_snp_test <- function(data1, trait, genotype, plotSize, loc,s, mkfold){
for(i in c(1:length(trait))){  

  test.geno <- vector("list",mkfold)

  base.prmse <- base.pcor <- base.fail <- model2.prmse <- model2.pcor <- model2.fail <- matrix(NA,mkfold,s)
  model3.prmse <- model3.pcor <- model3.fail <- model4.prmse <- model4.pcor <- model4.fail <- matrix(NA,mkfold,s)

  kvalue <- seq(0.1,2,0.3) # to compute incidence matrix using sigmoid function
  b.c.values <- matrix(0,3,2)
  b.c.values[1,] <- c(0.9,0.1)
  b.c.values[2,] <- c(0.99,0.01)
  b.c.values[3,] <- c(0.999,0.001)

  model1.prmse.1 <- model1.pcor.1 <- model1.fail.1 <- matrix(NA, mkfold, length(kvalue))
  model5.prmse.1 <- model5.pcor.1 <- model5.fail.1 <- matrix(NA,mkfold,nrow(b.c.values))
 # model6.prmse.1 <- model6.pcor.1 <- model6.fail.1 <- matrix(NA,mkfold,nrow(b.c.values))
  model7.prmse.1 <- model7.pcor.1 <- model7.fail.1 <- matrix(NA,mkfold,nrow(b.c.values))
  model8.prmse.1 <- model8.pcor.1 <- model8.fail.1 <- matrix(NA,mkfold,nrow(b.c.values))
  model9.prmse.1 <- model9.pcor.1 <- model9.fail.1 <- matrix(NA,mkfold,nrow(b.c.values))

   
  data1a <- na.omit(subset(data1,select = c("Range","Column",genotype,trait[i])))
  
  fixed.effect <- formula(paste(trait[i],"1", sep="~"))
  base <- lme(fixed = fixed.effect, random = as.formula(paste0("~ 1|", genotype)), data = na.omit(data1a), method = "ML")
  data1a <- (romr.fnc(base, na.omit(data1a), trim = 2.5))$data #outlier removal

  #K matrix aligning with data1a
  #subsetting K matrix  - helps to remove clones that are not genotyped but phenotyped 
  table(data1a[,genotype] %in% rownames(snps)) 
  snp2.names <- intersect(data1a[,genotype], rownames(snps))
  snps2 <- subset(snps, rownames(snps) %in% snp2.names)

  #calculate K matrix
  K <- A.mat(snps2-1)
     

  data2 <- subset(data1a,data1a[,genotype] %in% snp2.names)
    write.csv(data2, "data2.csv")
    data2 <- read.csv("data2.csv", h=T, sep=",")   


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
base.prmse.2 <- mean(base.prmse, na.rm = TRUE)
base.pcor.2 <- mean(base.pcor, na.rm = TRUE)
base.fail.2 <- sum(base.fail, na.rm = TRUE)


#model 1 - updating base model - with sigmoid function
 #computing influence matrix, Zc for competition from distance matrix and adding it to model 

  for(m in c(1:length(test.geno))){
    model1.prmse <- model1.pcor <- model1.fail <- matrix(NA,s,length(kvalue))

    for(c in c(1:length(kvalue))){
      distMat.1 <- kvalue[c]*(distMat.inv)/(kvalue[c]-(distMat.inv) + 1)
      #distMat.1[distMat.1 < 0.25] <- 0   
      distMat.1[distMat.1 == "NaN"] <- 0

      Z.comp <- distMat.1 %*% Z.base

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

            model1.prmse[s,c] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))  
            model1.pcor[s,c] <- cor(yhat.test, test.data[,trait[i]])
          }else{
            model1.fail[s,c] <- 1
          } # end of try
        } #end of one CV
      } # end of k-value
      model1.prmse.1[m,] <- apply(model1.prmse,2,mean, na.rm = TRUE)
      model1.pcor.1[m,] <- apply(model1.pcor,2,mean, na.rm = TRUE)
      model1.fail.1[m,] <- apply(model1.fail,2,sum, na.rm = TRUE)
    } # end of multiple CV 
print("model1 end")  
  model1.prmse.2 <- apply(model1.prmse.1,2, mean)
  k.1a <- which(model1.prmse.2 == min(model1.prmse.2))[1]
  k.win <- kvalue[k.1a]
  model1.prmse.3 <- model1.prmse.2[k.1a]

  model1.pcor.2 <- apply(model1.pcor.1,2,mean)
  model1.pcor.3 <- model1.pcor.2[k.1a]

  model1.fail.2 <- apply(model1.fail.1,2,sum) 
  
  

#model2 - with inverse distance matrix
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

      distMat.2 <- distMat.inv
      #distMat.2[distMat.2 < 0.25] <- 0
      distMat.2[distMat.2 == "Inf"] <- 0
    
      Z.comp <- distMat.2 %*% Z.base

      Z.comp.train <- Z.comp[rownames(Z.comp) %in% train.data$Slno,]
      Z.comp.test <- Z.comp[!(rownames(Z.comp) %in% train.data$Slno),]

      G.comp <- Z.comp.train%*%K%*%t(Z.comp.train) 
      model2 <- try(regress(y ~ X.train, ~G.geno + G.comp, pos= rep(TRUE,3), tol = 1e-4, data = train.data))

        if(class(model2) != "try-error"){

        R.residual <- Identity * model2$sigma[[3]]

        Khat.g <- K * model2$sigma[[1]]
        Khat.c <- K * model2$sigma[[2]]
        Vhat <- G.geno * model2$sigma[[1]] + G.comp * model2$sigma[[2]] + R.residual

        gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model2$fitted)
        zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model2$fitted)

        gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model2$beta)
        zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
          yhat.test <- zeta.test[,1] + zeta.test[,2] + zeta.test[,3] 

        model2.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
        model2.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
        }else{
          model2.fail[m,s] <- 1
         }#end of try
        } # end of one CV
        } # end of multiple CV 

      #for comparison and plotting
      model2.prmse.2 <- mean(model2.prmse, na.rm = TRUE)
      model2.pcor.2 <- mean(model2.pcor, na.rm = TRUE)
      model2.fail.2 <- sum(model2.fail, na.rm = TRUE)
print("model2 end")


#model3 - nearest neighbours compete irrespective of distance from four sides and diagonal
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

      distMat.3 <- distMat
      distMat.3[distMat.3 > sqrt(plotSize[1]^2 + plotSize[2]^2)] <- 0
      distMat.3[distMat.3 == plotSize[1]] <- 1
      distMat.3[distMat.3 == plotSize[2]] <- 0.5
      distMat.3[distMat.3 == sqrt(plotSize[1]^2 + plotSize[2]^2)] <- 0.2
      distMat.3[distMat.3 > 1] <- 0
    
      Z.comp <- distMat.3 %*% Z.base

      Z.comp.train <- Z.comp[rownames(Z.comp) %in% train.data$Slno,]
      Z.comp.test <- Z.comp[!(rownames(Z.comp) %in% train.data$Slno),]

      G.comp <- Z.comp.train%*%K%*%t(Z.comp.train) 
      model3 <- try(regress(y ~ X.train, ~G.geno + G.comp, pos= rep(TRUE,3), tol = 1e-4, data = train.data))

        if(class(model3) != "try-error"){

        R.residual <- Identity * model3$sigma[[3]]

        Khat.g <- K * model3$sigma[[1]]
        Khat.c <- K * model3$sigma[[2]]
        Vhat <- G.geno * model3$sigma[[1]] + G.comp * model3$sigma[[2]] + R.residual

        gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model3$fitted)
        zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model3$fitted)

        gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model3$beta)
        zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
          yhat.test <- zeta.test[,1] + zeta.test[,2] + zeta.test[,3] 

        model3.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
        model3.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
        }else{
          model3.fail[m,s] <- 1
         }#end of try
        } # end of one CV
        } # end of multiple CV 

      #for comparison and plotting
      model3.prmse.2 <- mean(model3.prmse, na.rm = TRUE)
      model3.pcor.2 <- mean(model3.pcor, na.rm = TRUE)
      model3.fail.2 <- sum(model3.fail, na.rm = TRUE)
print("model3 end")


#model4 - nearest neighbours compete irrespective of distance along the long side of plots
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

      distMat.4 <- distMat
      distMat.4[distMat.4 > plotSize[1]] <- 0
      distMat.4[distMat.4 == plotSize[1]] <- 1
      distMat.4[distMat.4 > 1] <- 0
    
      Z.comp <- distMat.4 %*% Z.base

      Z.comp.train <- Z.comp[rownames(Z.comp) %in% train.data$Slno,]
      Z.comp.test <- Z.comp[!(rownames(Z.comp) %in% train.data$Slno),]

      G.comp <- Z.comp.train%*%K%*%t(Z.comp.train) 
      model4 <- try(regress(y ~ X.train, ~G.geno + G.comp, pos= rep(TRUE,3), tol = 1e-4, data = train.data))

        if(class(model4) != "try-error"){

        R.residual <- Identity * model4$sigma[[3]]

        Khat.g <- K * model4$sigma[[1]]
        Khat.c <- K * model4$sigma[[2]]
        Vhat <- G.geno * model4$sigma[[1]] + G.comp * model4$sigma[[2]] + R.residual

        gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model4$fitted)
        zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model4$fitted)

        gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model4$beta)
        zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
          yhat.test <- zeta.test[,1] + zeta.test[,2] + zeta.test[,3] 

        model4.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
        model4.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
        }else{
          model4.fail[m,s] <- 1
         }#end of try
        } # end of one CV
        } # end of multiple CV 

      #for comparison and plotting
      model4.prmse.2 <- mean(model4.prmse, na.rm = TRUE)
      model4.pcor.2 <- mean(model4.pcor, na.rm = TRUE)
      model4.fail.2 <- sum(model4.fail, na.rm = TRUE)
print("model4 end")


#model 5 - updating base model - with modified sigmoid function - high competition among neighbours up to 3 ranges
 #computing influence matrix, Zc for competition from distance matrix and adding it to model 

  for(m in c(1:length(test.geno))){
    model5.prmse <- model5.pcor <- model5.fail <- matrix(NA,s,nrow(b.c.values))

    for(c in c(1:nrow(b.c.values))){
      distMat.2 <- distMat
      distMat.inv.2 <- 1/distMat.2
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.1 <- 1/(b.c.values[c,1] + b.c.values[c,2]^distMat.inv.2)
      distMat.1[is.na(distMat.1)] <- 0 
      diag(distMat.1) <- 0      
      
      Z.comp <- distMat.1 %*% Z.base

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

        model5 <- try(regress(y ~ X.train, ~G.geno + G.comp, pos= rep(TRUE,3), tol = 1e-4, data = train.data),silent = TRUE)

        if(class(model5) != "try-error"){
          R.residual <- Identity * model5$sigma[[3]]

          Khat.g <- K * model5$sigma[[1]]
          Khat.c <- K * model5$sigma[[2]]
          Vhat <- G.geno * model5$sigma[[1]] + G.comp * model5$sigma[[2]] + R.residual

          gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model5$fitted)
          zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model5$fitted)

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model5$beta)
          zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
            yhat.test <- zeta.test[,1] + zeta.test[,2] + zeta.test[,3] 

            model5.prmse[s,c] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))  
            model5.pcor[s,c] <- cor(yhat.test, test.data[,trait[i]])
          }else{
            model5.fail[s,c] <- 1
          } # end of try
        } #end of one CV
      } # end of b.c-value
      model5.prmse.1[m,] <- apply(model5.prmse,2,mean, na.rm = TRUE)
      model5.pcor.1[m,] <- apply(model5.pcor,2,mean, na.rm = TRUE)
      model5.fail.1[m,] <- apply(model5.fail,2,sum, na.rm = TRUE)
    } # end of multiple CV 
print("model5 end")  
  model5.prmse.2 <- apply(model5.prmse.1,2, mean)
  bc.1a <- which(model5.prmse.2 == min(model5.prmse.2))[1]
  bc.win <- b.c.values[bc.1a,]
  model5.prmse.3 <- model5.prmse.2[bc.1a]


  model5.pcor.2 <- apply(model5.pcor.1,2,mean)
  model5.pcor.3 <- model5.pcor.2[bc.1a]

  model5.fail.2 <- apply(model5.fail.1,2,sum)


#model 7 - updating base model - with modified sigmoid function - high competition among neighbours up to a column
 #computing influence matrix, Zc for competition from distance matrix and adding it to model 

  for(m in c(1:length(test.geno))){
    model7.prmse <- model7.pcor <- model7.fail <- matrix(NA,s,nrow(b.c.values))

    for(c in c(1:nrow(b.c.values))){
      distMat.2 <- distMat
      distMat.inv.2 <- 1/distMat.2
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.1 <- 1/(b.c.values[c,1] + b.c.values[c,2]^distMat.inv.2)
      a <- 1/plotSize[2]
      distMat.1[distMat.1 < (1/(b.c.values[c,1] + b.c.values[c,2]^a))] <- 0
      distMat.1[is.na(distMat.1)] <- 0
      diag(distMat.1) <- 0

      Z.comp <- distMat.1 %*% Z.base

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

        model7 <- try(regress(y ~ X.train, ~G.geno + G.comp, pos= rep(TRUE,3), tol = 1e-4, data = train.data),silent = TRUE)

        if(class(model7) != "try-error"){
          R.residual <- Identity * model7$sigma[[3]]

          Khat.g <- K * model7$sigma[[1]]
          Khat.c <- K * model7$sigma[[2]]
          Vhat <- G.geno * model7$sigma[[1]] + G.comp * model7$sigma[[2]] + R.residual

          gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model7$fitted)
          zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model7$fitted)

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model7$beta)
          zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
            yhat.test <- zeta.test[,1] + zeta.test[,2] + zeta.test[,3] 

            model7.prmse[s,c] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))  
            model7.pcor[s,c] <- cor(yhat.test, test.data[,trait[i]])
          }else{
            model7.fail[s,c] <- 1
          } # end of try
        } #end of one CV
      } # end of b.c-value
      model7.prmse.1[m,] <- apply(model7.prmse,2,mean, na.rm = TRUE)
      model7.pcor.1[m,] <- apply(model7.pcor,2,mean, na.rm = TRUE)
      model7.fail.1[m,] <- apply(model7.fail,2,sum, na.rm = TRUE)
    } # end of multiple CV 
print("model7 end")  
  model7.prmse.2 <- apply(model7.prmse.1,2, mean)
  bc.3a <- which(model7.prmse.2 == min(model7.prmse.2))[1]
  bc.win.3 <- b.c.values[bc.3a,]
  model7.prmse.3 <- model7.prmse.2[bc.3a]

  model7.pcor.2 <- apply(model7.pcor.1,2,mean)
  model7.pcor.3 <- model7.pcor.2[bc.3a]

  model7.fail.2 <- apply(model7.fail.1,2,sum)  


#model 8 - updating base model - with modified sigmoid function - high competition limited to 3*plotSize[1] 

  for(m in c(1:length(test.geno))){
    model8.prmse <- model8.pcor <- model8.fail <- matrix(NA,s,nrow(b.c.values))

    for(c in c(1:nrow(b.c.values))){
      distMat.2 <- distMat
      distMat.inv.2 <- 1/distMat.2
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.1 <- 1/(b.c.values[c,1] + b.c.values[c,2]^distMat.inv.2)
      a <- 1/(3*plotSize[1])
      distMat.1[distMat.1 < (1/(b.c.values[c,1] + b.c.values[c,2]^a))] <- 0
      distMat.1[is.na(distMat.1)] <- 0
      diag(distMat.1) <- 0
      
      Z.comp <- distMat.1 %*% Z.base

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

        model8 <- try(regress(y ~ X.train, ~G.geno + G.comp, pos= rep(TRUE,3), tol = 1e-4, data = train.data),silent = TRUE)

        if(class(model8) != "try-error"){
          R.residual <- Identity * model8$sigma[[3]]

          Khat.g <- K * model8$sigma[[1]]
          Khat.c <- K * model8$sigma[[2]]
          Vhat <- G.geno * model8$sigma[[1]] + G.comp * model8$sigma[[2]] + R.residual

          gamma <- Khat.g %*% t(Z.geno) %*% solve(Vhat) %*% (y - model8$fitted)
          zeta <- Khat.c %*% t(Z.comp.train) %*% solve(Vhat) %*% (y - model8$fitted)

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model8$beta)
          zeta.test <- cbind(Z.comp.test %*% zeta, gamma.test) # directly merge into gamma.test which is now in order of Slno.
            yhat.test <- zeta.test[,1] + zeta.test[,2] + zeta.test[,3] 

            model8.prmse[s,c] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))  
            model8.pcor[s,c] <- cor(yhat.test, test.data[,trait[i]])
          }else{
            model8.fail[s,c] <- 1
          } # end of try
        } #end of one CV
      } # end of b.c-value
      model8.prmse.1[m,] <- apply(model8.prmse,2,mean, na.rm = TRUE)
      model8.pcor.1[m,] <- apply(model8.pcor,2,mean, na.rm = TRUE)
      model8.fail.1[m,] <- apply(model8.fail,2,sum, na.rm = TRUE)
    } # end of multiple CV 
print("model8 end")  
  model8.prmse.2 <- apply(model8.prmse.1,2, mean)
  bc.4a <- which(model8.prmse.2 == min(model8.prmse.2))[1]
  bc.win.4 <- b.c.values[bc.4a,]
  model8.prmse.3 <- model8.prmse.2[bc.4a]

  model8.pcor.2 <- apply(model8.pcor.1,2,mean)
  model8.pcor.3 <- model8.pcor.2[bc.4a]

  model8.fail.2 <- apply(model8.fail.1,2,sum)  


#best model from each set
  models.prmse <- cbind(base.prmse.2,model1.prmse.3,model2.prmse.2, model3.prmse.2, model4.prmse.2, model5.prmse.3,
                          model7.prmse.3, model8.prmse.3)
  models.pcor <- cbind(base.pcor.2, model1.pcor.3, model2.pcor.2, model3.pcor.2, model4.pcor.2, model5.pcor.3,
                        model7.pcor.3, model8.pcor.3)
  models.fail <- cbind(base.fail.2, model1.fail.2, model2.fail.2, model3.fail.2, model4.fail.2, model5.fail.2,
                        model7.fail.2, model8.fail.2)

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
  write.csv(models.prmse, paste0(trait[i],"-",loc,"_K_models_prmse_comp.csv")) 
  write.csv(models.pcor, paste0(trait[i],"-",loc,"_K_models_pcor_comp.csv")) 
  write.csv(models.fail, paste0(trait[i],"-",loc,"_K_models_fail_comp.csv")) 

#correlation and design matrices for final model

  y <- data2[,trait[i]]

  #X design matrix 
  X.mat <- model.matrix(fixed.effect,data=data2)
 
  #genotypic design matrix
  idg <- factor(as.character(data2[,genotype]), levels = rownames(K))
    Z.geno <- model.matrix(~idg - 1)
    geno <- Z.geno%*%K%*%t(Z.geno) 

  #competition matrices
   distMat.1 <- k.win*(distMat.inv)/(k.win-(distMat.inv) + 1)
    distMat.1[distMat.1 == "NaN"] <- 0

      Z.comp.1 <- distMat.1 %*% Z.base
      G.comp.1 <- Z.comp.1%*%K%*%t(Z.comp.1) 

   distMat.2 <- distMat.inv
    distMat.2[distMat.2 == "Inf"] <- 0
    
      Z.comp.2 <- distMat.2 %*% Z.base
      G.comp.2 <- Z.comp.2%*%K%*%t(Z.comp.2)

   distMat.3 <- distMat
    distMat.3[distMat.3 > sqrt(plotSize[1]^2 + plotSize[2]^2)] <- 0
    distMat.3[distMat.3 == plotSize[1]] <- 1
    distMat.3[distMat.3 == plotSize[2]] <- 0.5
    distMat.3[distMat.3 == sqrt(plotSize[1]^2 + plotSize[2]^2)] <- 0.2

      Z.comp.3 <- distMat.3 %*% Z.base
      G.comp.3 <- Z.comp.3%*%K%*%t(Z.comp.3)

   distMat.4 <- distMat
    distMat.4[distMat.4 > plotSize[1]] <- 0
    distMat.4[distMat.4 == plotSize[1]] <- 1

      Z.comp.4 <- distMat.4 %*% Z.base
      G.comp.4 <- Z.comp.4%*%K%*%t(Z.comp.4)

   distMat.5 <- distMat
      distMat.inv.2 <- 1/distMat.5
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.5a <- 1/(b.c.values[bc.1a,1] + b.c.values[bc.1a,2]^distMat.inv.2)
      distMat.5a[is.na(distMat.5a)] <- 0 
      diag(distMat.5a) <- 0  
      
      Z.comp.5 <- distMat.5a %*% Z.base
      G.comp.5 <- Z.comp.5%*%K%*%t(Z.comp.5)

   
   distMat.7 <- distMat
      distMat.inv.2 <- 1/distMat.7
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.7a <- 1/(b.c.values[bc.3a,1] + b.c.values[bc.3a,2]^distMat.inv.2)
      a <- 1/plotSize[2]
      distMat.7a[distMat.7a < (1/(b.c.values[bc.3a,1] + b.c.values[bc.3a,2]^a))] <- 0
      distMat.7a[is.na(distMat.7a)] <- 0
      diag(distMat.7a) <- 0       
      
      Z.comp.7 <- distMat.7a %*% Z.base
      G.comp.7 <- Z.comp.7%*%K%*%t(Z.comp.7)   


   distMat.8 <- distMat
      distMat.inv.2 <- 1/distMat.8
      distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
      distMat.8a <- 1/(b.c.values[bc.4a,1] + b.c.values[bc.4a,2]^distMat.inv.2)
      a <- 1/(3*plotSize[1])
      distMat.8a[distMat.8a < (1/(b.c.values[bc.4a,1] + b.c.values[bc.4a,2]^a))] <- 0
      distMat.8a[is.na(distMat.8a)] <- 0
      diag(distMat.8a) <- 0
              
      Z.comp.8 <- distMat.8a %*% Z.base
      G.comp.8 <- Z.comp.8%*%K%*%t(Z.comp.8)   

 
       
         
#choosing the final model

 model <- switch(model.win.no,
    "1" <- regress(y ~ X.mat, ~geno, pos= rep(TRUE,2), tol = 1e-4, data = data2),
    "2" <- regress(y ~ X.mat, ~geno + G.comp.1, pos= rep(TRUE,3), tol = 1e-4, data = data2),
    "3" <- regress(y ~ X.mat, ~geno + G.comp.2, pos= rep(TRUE,3), tol = 1e-4, data = data2),
    "4" <- regress(y ~ X.mat, ~geno + G.comp.3, pos= rep(TRUE,3), tol = 1e-4, data = data2),
    "5" <- regress(y ~ X.mat, ~geno + G.comp.4, pos= rep(TRUE,3), tol = 1e-4, data = data2),
    "6" <- regress(y ~ X.mat, ~geno + G.comp.5, pos= rep(TRUE,3), tol = 1e-4, data = data2),
    "7" <- regress(y ~ X.mat, ~geno + G.comp.7, pos= rep(TRUE,3), tol = 1e-4, data = data2),
    "8" <- regress(y ~ X.mat, ~geno + G.comp.8, pos= rep(TRUE,3), tol = 1e-4, data = data2),
       
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
 cat(base.formula, file = paste0(trait[i],"-",loc,"_K_base_summary_comp.txt"), sep = "\n", append=TRUE)
 cat(base.Vformula, file = paste0(trait[i],"-",loc,"_K_base_summary_comp.txt"), sep = "\n", append=TRUE )
 cat(base.out, file= paste0(trait[i],"-",loc,"_K_base_summary_comp.txt"), sep="\n", append=TRUE)

 cat(model.formula, file = paste0(trait[i],"-",loc,"_K_model_summary_comp.txt"), sep = "\n", append=TRUE)
 cat(model.Vformula, file = paste0(trait[i],"-",loc,"_K_model_summary_comp.txt"), sep = "\n", append=TRUE )
 cat(paste0("k.value =", k.win), file = paste0(trait[i],"-",loc,"_K_model_summary_comp.txt"), sep = "\n", append=TRUE)
 cat(paste0("b.c.value =", bc.win, bc.win.3, bc.win.4), file = paste0(trait[i],"-",loc,"_K_model_summary_comp.txt"), sep = "\n", append=TRUE)
 #cat(paste0("phi value =", bc.win.2, bc.win.6, bc.win.7), file = paste0(trait[i],"-",loc,"_K_model_summary_comp.txt"), sep = "\n", append=TRUE)
 cat(model.out, file= paste0(trait[i],"-",loc,"_K_model_summary_comp.txt"), sep="\n", append=TRUE)

} # end of for w.r.t. trait  
} # end of function



###to get correct competitive BLUP from model - demonstration

#blup from base and models
  #baseBLUP <- BLUP(base)$Mean
  #modelBLUP <- BLUP(model)$Mean
  #model2BLUP <- BLUP(model2)$Mean
  

  #separating spatial and entry blups for models
   #model.inflBLUP <- modelBLUP[grep("comp", names(modelBLUP), fixed=TRUE)]
      #data2.sp <- cbind(data2, model.inflBLUP)

  #base.EntryBLUP <- baseBLUP[grep("geno", names(baseBLUP), fixed=TRUE)]
    
  #model.EntryBLUP <- modelBLUP[grep("geno", names(modelBLUP), fixed=TRUE)]

 #data2.entry <- cbind(data2.sp,base.EntryBLUP, model.EntryBLUP) 

  #modelcomp <- ginv(Z.comp.1) %*% model.inflBLUP
    #rownames(modelcomp) <- colnames(Z.comp.1)
    #modelcomp <- as.data.frame(modelcomp)
    #modelcomp <- setDT(modelcomp, keep.rownames = TRUE)[]
    #colnames(modelcomp) <- c("CLONE","modelcomp")
    #modelcomp$CLONE <- gsub("id","",modelcomp$CLONE)
    #data2.entry <- merge(data2.entry,modelcomp)
       
 #output BLUP
   #write.csv(data2.entry, paste0(trait[i],"-",loc,"_K_genotype_comp_BLUP.csv"))
  