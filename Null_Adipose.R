###Analysis code for Metformin mouse dataset
#######This code fits GLM by only accounting for technical variability 

library(dplyr)
library(pscl)
library(MASS)
library(dplyr)
library(data.table)
library(emdbook)
library(future.apply)
library(Matrix)


###functions for running the models

#function for calculating library size the library size
Library_size <- function(y){
  apply(y,2, function(x) sum(x))
}


#function for calculating mouse ID
Mouse_ID <- function(z){
  mou.1 <- apply( z[ , grepl(as.vector(colnames(z)), pattern = "1F_")] , 2 , function(x) rep(1) )
  mou.2 <- apply( z[ , grepl(as.vector(colnames(z)), pattern = "2F_")] , 2 , function(x) rep(2) )
  mou.3 <- apply( z[ , grepl(as.vector(colnames(z)), pattern = "3F_")] , 2 , function(x) rep(3) )
  mou.4 <- apply( z[ , grepl(as.vector(colnames(z)), pattern = "4F_")] , 2 , function(x) rep(4) )
  
  mou.id <- as.factor(c(mou.1, mou.2, mou.3, mou.4)) 
}

#function for performing KS test to select genes belonging to the family of ZINB distributions
KS_ZINB <- function(x, lib.size, mou.id){
  
  library(pscl)
  m1 <- try(zeroinfl(x ~ 1 | mou.id, offset=log(lib.size), dist = "negbin"), silent = TRUE)
  
  if(!(class(m1) == "try-error")){
    pi_ML = predict(m1, type = "zero")
    theta_ML = m1$theta
    mean_ML = predict(m1, type = "count")
    var_ML = mean_ML + (mean_ML ^ 2 / theta_ML)
    ccc = rbind(pi_ML, theta_ML, mean_ML)
  }
  else {
    ccc = "NA"
  }  
  
  
  if(!(class(ccc) == "character")){
    library(VGAM)
    pp <- try(rzinegbin(n = length(x), size = ccc[2,], munb = ccc[3,], pstr0 = ccc[1,]), silent = TRUE)
  }
  else {
    pp <- "NA"
  }  
  
  if(!(class(pp) == "character")){ 
    
    library(dgof)  
    D <- try(ks.test(x, ecdf(pp), simulate.p.value = TRUE)$p.value, silent = TRUE)
  }
  else {
    D <- "NA"
  }
  
  if(!(class(D) == "character")){
    
    p_value <- D
  }
  else {
    p_value <- "NA"
  }
  
  return(p_value) 
}



#functions for fitting the 4 distributions
#Fit Poisson distribution
model_poi <- function(x, lib.size, mou.id){
  
  stats::glm(x ~ mou.id, offset=log(lib.size), family = "poisson")
}


model_nb <- function(x, lib.size, mou.id){
  library(MASS)
  try(glm.nb(x ~ mou.id + offset(log(lib.size))), silent = TRUE)
}


model_zip <- function(x, lib.size, mou.id){
  library(pscl)
  try(zeroinfl(x ~ 1 | mou.id, offset=log(lib.size), dist = "poisson"), silent = TRUE)
}


model_zinb <- function(x, lib.size, mou.id){
  library(pscl)
  try(zeroinfl(x ~ 1 | mou.id, offset=log(lib.size), dist = "negbin"), silent = TRUE)
}


#Calculate BIC
model_BIC <- function(z){
  if(!(class(z)=="try-error")){
    stats::BIC(z)
  }
  else{
    "NA"
  }
}


#function for selecting the model with the least BIC value
#Selecting best model and subsetting genes for GOF test
Best_model <- function(y,x){
  
  
  #return column name of min value for each row
  BIC_colmax <- colnames(y)[apply(y,1,function(x) which.min(x[x>0]))]
  
  #create data frame
  BIC_dataframe <- as.data.frame(BIC_colmax)
  
  #Get the total count of genes follwoing each distribution
  BIC_dataframe %>% group_by(BIC_colmax) %>% tally()
  
  
  #Create data frame of minimum BIC values with gene names
  BIC_genename <- as.data.frame(BIC_dataframe)
  row.names(BIC_genename) <- names(x)
  
  #grep the names of genes following each distribution
  poi_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^Poi_BIC")]
  nb_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^NB_BIC")]
  zip_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^ZIP_BIC")]
  zinb_genes <- rownames(BIC_genename)[grep(BIC_genename[,1], pattern = "^ZINB_BIC")]
  
  poi_MGenes <- x[names(x) %in% poi_genes]
  nb_MGenes <- x[names(x) %in% nb_genes]
  zip_MGenes <- x[names(x) %in% zip_genes]
  zinb_MGenes <- x[names(x) %in% zinb_genes]
  
  return(list(poi_MGenes, nb_MGenes, zip_MGenes, zinb_MGenes))
}


#function for GOF test using LR test and mixed LR test
#Fit the models
#Poisson distribution
model_Chipoi <- function(x, lib.size, mou.id){
  
  library(stats)
  
  1-pchisq(summary(glm(x ~ mou.id, offset=log(lib.size), family = "poisson"))$deviance, 
           df= summary(glm(x ~ mou.id, offset=log(lib.size), family = "poisson"))$df.residual)
}

#NB distribution
model_Chinb <- function(x, lib.size, mou.id){
  
  library(MASS)
  
  1-pchisq(summary(glm.nb(x ~ mou.id + offset(log(lib.size))))$deviance, 
           df= summary(glm.nb(x ~ mou.id + offset(log(lib.size))))$df.residual)
}

#ZIP distribution
model_Chizip <- function(x, lib.size, mou.id){
  library(pscl)
  library(emdbook)
  
  1-pchibarsq(2*(logLik(zeroinfl(x ~ 1 | mou.id, offset=log(lib.size), dist = "poisson"))-logLik(glm(x ~ mou.id, offset=log(lib.size), family = "poisson"))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)
  
}


#ZINB distribution
model_Chizinb <- function(x, lib.size, mou.id){
  library(pscl)
  library(emdbook)
  
  1-pchibarsq(2*(logLik(zeroinfl(x ~ 1 | mou.id, offset=log(lib.size), dist = "negbin"))-logLik(glm.nb(x ~ mou.id + offset(log(lib.size))))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)
}




#########################################################################################################################
####run the script for the three conditions seperately
####load cell-type rds object
load("METF/Cell-types/Adipose_celltypes.rds")

#convert to matrix form
O_fat_Cells <- subset(Adipose_celltypes, subset = condition == "Old")
O_fat_count <- as.matrix(O_fat_Cells@assays$RNA@counts)

###convert data matrix to sparse matrix
O_fat_count <- Matrix(O_fat_count, sparse = TRUE)

#filter genes not expressed in at least 10% of the cells
O_fat_filt <- O_fat_count[apply(O_fat_count, 1, function(x){sum(x == 0)} < ncol(O_fat_count)*0.90),]


########Old
### calculate library size
lib.size <- Library_size(O_fat_filt)

###calculte mouse ID
mou.id <- Mouse_ID(O_fat_filt)

plan(multisession, workers = 24)

#perform the KS test
O_fat_KS <- future_apply(O_fat_filt, MARGIN = 1L, FUN = KS_ZINB, lib.size <- lib.size, future.seed = 0xBEEF)
save(O_fat_KS, file = "METF/Null_data/Adipose/O_fat/O_fat_KS")

#select genes passing the KS test
O_fat_KS <- stack(O_fat_KS)
O_fat_KS$values <- as.numeric(O_fat_KS$values)
O_fat_KS <- O_fat_KS[!is.na(O_fat_KS$values), ]

O_fat_KSig <- O_fat_KS[O_fat_KS$values > 0.05, ]
O_fat_genes <- O_fat_KSig$ind

O_fat_data <- O_fat_filt[rownames(O_fat_filt) %in% O_fat_genes,]


#fit models and save the models
O_fat_poi <- future_apply(O_fat_data, MARGIN = 1L, FUN = model_poi, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
O_fat_nb <- future_apply(O_fat_data, MARGIN = 1L, FUN = model_nb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
O_fat_zip <- future_apply(O_fat_data, MARGIN = 1L, FUN = model_zip, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
O_fat_zinb <- future_apply(O_fat_data, MARGIN = 1L, FUN = model_zinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)

#Calculate BIC
O_fat_poi_BIC <- t(as.data.frame(lapply(O_fat_poi, model_BIC)))
O_fat_nb_BIC <- t(as.data.frame(lapply(O_fat_nb, model_BIC)))
O_fat_zip_BIC <- t(as.data.frame(lapply(O_fat_zip, model_BIC)))
O_fat_zinb_BIC <- t(as.data.frame(lapply(O_fat_zinb, model_BIC)))


O_fat_BIC <- cbind(O_fat_poi_BIC, O_fat_nb_BIC, O_fat_zip_BIC, O_fat_zinb_BIC)
colnames(O_fat_BIC) <- c("Poi_BIC", "NB_BIC", "ZIP_BIC", "ZINB_BIC")
save(O_fat_BIC, file = "METF/Null_data/Adipose/O_fat/O_fat_BIC.rds")


##convert data matrix to list for selecting the model with the least BIC value
O_fat_list <- lapply(as.list(1:dim(O_fat_data)[1]), function(x) O_fat_data[x[1],])
names(O_fat_list) <- rownames(O_fat_data)

#Select the model with the least BIC value
BM_genes <- Best_model(O_fat_BIC, O_fat_list)


O_fat_Chipoi_pval <- future_lapply(BM_genes[[1]], model_Chipoi, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
O_fat_Chinb_pval <- future_lapply(BM_genes[[2]], model_Chinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
O_fat_Chizip_pval <- future_lapply(BM_genes[[3]], model_Chizip, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
O_fat_Chizinb_pval <- future_lapply(BM_genes[[4]], model_Chizinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)


save(O_fat_Chipoi_pval, file = "METF/Null_data/Adipose/O_fat/O_fat_Chipoi_pval.rds")
save(O_fat_Chinb_pval, file = "METF/Null_data/Adipose/O_fat/O_fat_Chinb_pval.rds")
save(O_fat_Chizip_pval, file = "METF/Null_data/Adipose/O_fat/O_fat_Chizip_pval.rds")
save(O_fat_Chizinb_pval, file = "METF/Null_data/Adipose/O_fat/O_fat_Chizinb_pval.rds")



##################################################################################################
########Young
load("METF/Cell-types/Adipose_celltypes.rds")

#convert to matrix form
Y_fat_Cells <- subset(Adipose_celltypes, subset = condition == "Young")
Y_fat_count <- as.matrix(Y_fat_Cells@assays$RNA@counts)

###convert data matrix to sparse matrix
Y_fat_count <- Matrix(Y_fat_count, sparse = TRUE)

#filter genes not expressed in at least 10% of the cells
Y_fat_filt <- Y_fat_count[apply(Y_fat_count, 1, function(x){sum(x == 0)} < ncol(Y_fat_count)*0.90),]


########Young
### calculate library size
lib.size <- Library_size(Y_fat_filt)

###calculte mouse ID
mou.id <- Mouse_ID(Y_fat_filt)

plan(multisession, workers = 24)

#perform the KS test
Y_fat_KS <- future_apply(Y_fat_filt, MARGIN = 1L, FUN = KS_ZINB, lib.size <- lib.size, future.seed = 0xBEEF)
save(Y_fat_KS, file = "METF/Null_data/Adipose/Y_fat/Y_fat_KS")

#select genes passing the KS test
Y_fat_KS <- stack(Y_fat_KS)
Y_fat_KS$values <- as.numeric(Y_fat_KS$values)
Y_fat_KS <- Y_fat_KS[!is.na(Y_fat_KS$values), ]

Y_fat_KSig <- Y_fat_KS[Y_fat_KS$values > 0.05, ]
Y_fat_genes <- Y_fat_KSig$ind

Y_fat_data <- Y_fat_filt[rownames(Y_fat_filt) %in% Y_fat_genes,]


#fit models and save the models
Y_fat_poi <- future_apply(Y_fat_data, MARGIN = 1L, FUN = model_poi, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
Y_fat_nb <- future_apply(Y_fat_data, MARGIN = 1L, FUN = model_nb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
Y_fat_zip <- future_apply(Y_fat_data, MARGIN = 1L, FUN = model_zip, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
Y_fat_zinb <- future_apply(Y_fat_data, MARGIN = 1L, FUN = model_zinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)

#Calculate BIC
Y_fat_poi_BIC <- t(as.data.frame(lapply(Y_fat_poi, model_BIC)))
Y_fat_nb_BIC <- t(as.data.frame(lapply(Y_fat_nb, model_BIC)))
Y_fat_zip_BIC <- t(as.data.frame(lapply(Y_fat_zip, model_BIC)))
Y_fat_zinb_BIC <- t(as.data.frame(lapply(Y_fat_zinb, model_BIC)))


Y_fat_BIC <- cbind(Y_fat_poi_BIC, Y_fat_nb_BIC, Y_fat_zip_BIC, Y_fat_zinb_BIC)
colnames(Y_fat_BIC) <- c("Poi_BIC", "NB_BIC", "ZIP_BIC", "ZINB_BIC")
save(Y_fat_BIC, file = "METF/Null_data/Adipose/Y_fat/Y_fat_BIC.rds")


##convert data matrix to list for selecting the model with the least BIC value
Y_fat_list <- lapply(as.list(1:dim(Y_fat_data)[1]), function(x) Y_fat_data[x[1],])
names(Y_fat_list) <- rownames(Y_fat_data)

#Select the model with the least BIC value
BM_genes <- Best_model(Y_fat_BIC, Y_fat_list)


Y_fat_Chipoi_pval <- future_lapply(BM_genes[[1]], model_Chipoi, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
Y_fat_Chinb_pval <- future_lapply(BM_genes[[2]], model_Chinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
Y_fat_Chizip_pval <- future_lapply(BM_genes[[3]], model_Chizip, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
Y_fat_Chizinb_pval <- future_lapply(BM_genes[[4]], model_Chizinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)


save(Y_fat_Chipoi_pval, file = "METF/Null_data/Adipose/Y_fat/Y_fat_Chipoi_pval.rds")
save(Y_fat_Chinb_pval, file = "METF/Null_data/Adipose/Y_fat/Y_fat_Chinb_pval.rds")
save(Y_fat_Chizip_pval, file = "METF/Null_data/Adipose/Y_fat/Y_fat_Chizip_pval.rds")
save(Y_fat_Chizinb_pval, file = "METF/Null_data/Adipose/Y_fat/Y_fat_Chizinb_pval.rds")



##################################################################################################
########Treated
load("METF/Cell-types/Adipose_celltypes.rds")

#convert to matrix form
T_fat_Cells <- subset(Adipose_celltypes, subset = condition == "Treated")
T_fat_count <- as.matrix(T_fat_Cells@assays$RNA@counts)

###convert data matrix to sparse matrix
T_fat_count <- Matrix(T_fat_count, sparse = TRUE)

#filter genes not expressed in at least 10% of the cells
T_fat_filt <- T_fat_count[apply(T_fat_count, 1, function(x){sum(x == 0)} < ncol(T_fat_count)*0.90),]


########Treated
### calculate library size
lib.size <- Library_size(T_fat_filt)

###calculte mouse ID
mou.id <- Mouse_ID(T_fat_filt)

plan(multisession, workers = 24)

#perform the KS test
T_fat_KS <- future_apply(T_fat_filt, MARGIN = 1L, FUN = KS_ZINB, lib.size <- lib.size, future.seed = 0xBEEF)
save(T_fat_KS, file = "METF/Null_data/Adipose/T_fat/T_fat_KS")

#select genes passing the KS test
T_fat_KS <- stack(T_fat_KS)
T_fat_KS$values <- as.numeric(T_fat_KS$values)
T_fat_KS <- T_fat_KS[!is.na(T_fat_KS$values), ]

T_fat_KSig <- T_fat_KS[T_fat_KS$values > 0.05, ]
T_fat_genes <- T_fat_KSig$ind

T_fat_data <- T_fat_filt[rownames(T_fat_filt) %in% T_fat_genes,]


#fit models and save the models
T_fat_poi <- future_apply(T_fat_data, MARGIN = 1L, FUN = model_poi, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
T_fat_nb <- future_apply(T_fat_data, MARGIN = 1L, FUN = model_nb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
T_fat_zip <- future_apply(T_fat_data, MARGIN = 1L, FUN = model_zip, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
T_fat_zinb <- future_apply(T_fat_data, MARGIN = 1L, FUN = model_zinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)

#Calculate BIC
T_fat_poi_BIC <- t(as.data.frame(lapply(T_fat_poi, model_BIC)))
T_fat_nb_BIC <- t(as.data.frame(lapply(T_fat_nb, model_BIC)))
T_fat_zip_BIC <- t(as.data.frame(lapply(T_fat_zip, model_BIC)))
T_fat_zinb_BIC <- t(as.data.frame(lapply(T_fat_zinb, model_BIC)))


T_fat_BIC <- cbind(T_fat_poi_BIC, T_fat_nb_BIC, T_fat_zip_BIC, T_fat_zinb_BIC)
colnames(T_fat_BIC) <- c("Poi_BIC", "NB_BIC", "ZIP_BIC", "ZINB_BIC")
save(T_fat_BIC, file = "METF/Null_data/Adipose/T_fat/T_fat_BIC.rds")


##convert data matrix to list for selecting the model with the least BIC value
T_fat_list <- lapply(as.list(1:dim(T_fat_data)[1]), function(x) T_fat_data[x[1],])
names(T_fat_list) <- rownames(T_fat_data)

#Select the model with the least BIC value
BM_genes <- Best_model(T_fat_BIC, T_fat_list)


T_fat_Chipoi_pval <- future_lapply(BM_genes[[1]], model_Chipoi, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
T_fat_Chinb_pval <- future_lapply(BM_genes[[2]], model_Chinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
T_fat_Chizip_pval <- future_lapply(BM_genes[[3]], model_Chizip, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)
T_fat_Chizinb_pval <- future_lapply(BM_genes[[4]], model_Chizinb, lib.size <- lib.size, mou.id <- mou.id, future.seed = 0xBEEF)


save(T_fat_Chipoi_pval, file = "METF/Null_data/Adipose/T_fat/T_fat_Chipoi_pval.rds")
save(T_fat_Chinb_pval, file = "METF/Null_data/Adipose/T_fat/T_fat_Chinb_pval.rds")
save(T_fat_Chizip_pval, file = "METF/Null_data/Adipose/T_fat/T_fat_Chizip_pval.rds")
save(T_fat_Chizinb_pval, file = "METF/Null_data/Adipose/T_fat/T_fat_Chizinb_pval.rds")

