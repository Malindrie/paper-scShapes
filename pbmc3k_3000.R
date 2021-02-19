library(dplyr)
library(pscl)
library(MASS)
library(dplyr)
library(data.table)
library(emdbook)
library(future.apply)
library(Matrix)
library(sprof)

load("pbmc3k/pbmc_Chipoi_pval.rds")
load("pbmc3k/pbmc_Chinb_pval.rds")
load("pbmc3k/pbmc_Chizip_pval.rds")
load("pbmc3k/pbmc_Chizinb_pval.rds")


pbmc_data <- readRDS("pbmc3k/pbmc_data.rds")

pbmc_data <- pbmc_data[,apply(pbmc_data, 2, function(x) !all(x ==0))]
pbmc_data <- pbmc_data[apply(pbmc_data, 1, function(x) {sum(x == 0) < ncol(pbmc_data)*0.90}),]
pbmc_data <- as.matrix(pbmc_data)


###convert data matrix to sparse matrix
pbmc_filt <- Matrix(pbmc_data, sparse = TRUE)


poi_genes <- names(pbmc_Chipoi_pval)[pbmc_Chipoi_pval > 0.05]
nb_genes <- names(pbmc_Chinb_pval)[pbmc_Chinb_pval > 0.05]
zip_genes <- names(pbmc_Chizip_pval)[pbmc_Chizip_pval < 0.05]
zinb_genes <- names(pbmc_Chizinb_pval)[pbmc_Chizinb_pval < 0.05]


set.seed(0xBEEF)

poi_matrix <- pbmc_filt[rownames(pbmc_filt) %in% poi_genes,]
nb_matrix <- pbmc_filt[rownames(pbmc_filt) %in% nb_genes,]
zip_matrix <- pbmc_filt[rownames(pbmc_filt) %in% zip_genes,]
zinb_matrix <- pbmc_filt[rownames(pbmc_filt) %in% zinb_genes,]


#Compute the library size
Library_size <- function(y){
  apply(y,2, function(x) sum(x))
}

lib.size <- Library_size(pbmc_filt)

N <- 3000
n.umi <- sort(sample(x = lib.size, size = N, replace = TRUE))

#Fit the model and simulate
#poi
model_poi <- function(x, lib.size){
  
  stats::glm(x ~ 1, offset=log(lib.size), family = "poisson")
}


pbmc_poi <- future_apply(poi_matrix, 1L, FUN = model_poi, lib.size <- lib.size, future.seed = 0xBEEF)
poi_coefficients <- lapply(pbmc_poi, function(x) x$coefficients)
poi_coeff <- as.matrix(unlist(poi_coefficients))
rownames(poi_coeff) <- names(pbmc_poi)

poi.regg <- n.umi
poi.mu <- exp(poi_coeff) %*% poi.regg

poi.sim <- t(sapply(rownames(poi.mu), function(gene) {
  gene.mu <- poi.mu[gene, ]
  return(stats::rpois(n = length(gene.mu), lambda = gene.mu))
}))


#nb
model_nb <- function(x, lib.size){
  library(MASS)
  try(glm.nb(x ~ 1 + offset(log(lib.size))), silent = TRUE)
}

pbmc_nb <- future_apply(nb_matrix, 1L, FUN = model_nb, lib.size <- lib.size, future.seed = 0xBEEF)
nb_coefficients <- lapply(pbmc_nb, function(x) x$coefficients)
nb_coeff <- as.matrix(unlist(nb_coefficients))
rownames(nb_coeff) <- names(pbmc_nb)

nb.regg <- n.umi
nb.mu <- exp(nb_coeff) %*% nb.regg

nb_theta <- lapply(pbmc_nb, function(x) x$theta)
nb_theta <- as.matrix(unlist(nb_theta))


nb.sim <- t(sapply(rownames(nb.mu), function(gene) {
  gene.mu <- nb.mu[gene, ]
  gene.theta <- nb_theta[gene,]
  
  return(MASS::rnegbin(n = length(gene.mu), mu = gene.mu, theta = gene.theta))
  
}))



#zip
model_zip <- function(x, lib.size){
  library(pscl)
  try(zeroinfl(x ~ 1 | 1, offset=log(lib.size), dist = "poisson"), silent = TRUE)
}

pbmc_zip <- future_apply(zip_matrix, 1L, FUN = model_zip, lib.size <- lib.size, future.seed = 0xBEEF)
zip_coefficients <- lapply(pbmc_zip, function(x) x$coefficients)

zip_coeff <- NULL
zip_coeff_01 <- NULL
for (i in 1:length(zip_coefficients)) {
  zip_coeff_01 <- as.matrix(unlist(zip_coefficients[[i]]$count))  
  zip_coeff <- rbind(zip_coeff, zip_coeff_01)
}


rownames(zip_coeff) <- names(pbmc_zip)

zip_pi_01 <- NULL
zip.pi_ML <- NULL

for (i in 1:length(pbmc_zip)) {
  zip_pi_01 <- predict(pbmc_zip[[i]], type = "zero") 
  zip.pi_ML <- rbind(zip.pi_ML, zip_pi_01)
}

rownames(zip.pi_ML) <- names(pbmc_zip)


zip.pi_ML = lapply(pbmc_zip, function(x) predict(x, type = "zero")[1])
zip.pi_ML <- as.matrix(unlist(zip.pi_ML))
rownames(zip.pi_ML) <- names(pbmc_zip)

zip.regg <- n.umi
zip.mu <- exp(zip_coeff) %*% zip.regg

zip.sim <- t(sapply(rownames(zip.mu), function(gene) {
  gene.mu <- zip.mu[gene, ]
  gene.pi <- zip.pi_ML[gene,]
  
  return(VGAM::rzipois(n = length(gene.mu), lambda = gene.mu, pstr0 = gene.pi))
}))



#zinb
model_zinb <- function(x, lib.size){
  library(pscl)
  try(zeroinfl(x ~ 1 | 1, offset=log(lib.size), dist = "negbin"), silent = TRUE)
}

pbmc_zinb <- future_apply(zinb_matrix, 1L, FUN = model_zinb, lib.size <- lib.size, future.seed = 0xBEEF)
zinb_coefficients <- lapply(pbmc_zinb, function(x) x$coefficients)

zinb_coeff <- NULL
zinb_coeff_01 <- NULL
for (i in 1:length(zinb_coefficients)) {
  zinb_coeff_01 <- as.matrix(unlist(zinb_coefficients[[i]]$count))  
  zinb_coeff <- rbind(zinb_coeff, zinb_coeff_01)
}


rownames(zinb_coeff) <- names(pbmc_zinb)

zinb_pi_01 <- NULL
zinb.pi_ML <- NULL

for (i in 1:length(pbmc_zinb)) {
  zinb_pi_01 <- predict(pbmc_zinb[[i]], type = "zero") 
  zinb.pi_ML <- rbind(zinb.pi_ML, zinb_pi_01)
}

rownames(zinb.pi_ML) <- names(pbmc_zinb)


zinb.pi_ML = lapply(pbmc_zinb, function(x) predict(x, type = "zero")[1])
zinb.pi_ML <- as.matrix(unlist(zinb.pi_ML))
rownames(zinb.pi_ML) <- names(pbmc_zinb)

zinb.regg <- n.umi
zinb.mu <- exp(zinb_coeff) %*% zinb.regg

zinb_theta <- lapply(pbmc_zinb, function(x) x$theta)
zinb_theta <- as.matrix(unlist(zinb_theta))

zinb.sim <- t(sapply(rownames(zinb.mu), function(gene) {
  gene.mu <- zinb.mu[gene, ]
  gene.pi <- zinb.pi_ML[gene,]
  gene.theta <- zinb_theta[gene,]
  
  return(VGAM::rzinegbin(n = length(gene.mu), size = gene.theta, munb = gene.mu, pstr0 = gene.pi))
  
}))

result <- as.matrix(rbind(poi.sim, nb.sim, zip.sim, zinb.sim))

save(result, file = "pbmc3k/n3000/result.rds")

lib.size <- Library_size(result)


#Zero Infalted Negative Binomial
KS_ZINB <- function(x, lib.size){
  
  library(pscl)
  m1 <- try(zeroinfl(x ~ 1 | 1, offset=log(lib.size), dist = "negbin"), silent = TRUE)
  
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

plan(multisession, workers = 24)

result_KS <- future_apply(result, MARGIN = 1L, FUN = KS_ZINB, lib.size <- lib.size, future.seed = 0xBEEF)


save(result_KS, file = "pbmc3k/n3000/result_KS.rds")



##########################################################################################################
#Fit the models
model_poi <- function(x, lib.size){
  
  stats::glm(x ~ 1, offset=log(lib.size), family = "poisson")
}


model_nb <- function(x, lib.size){
  library(MASS)
  try(glm.nb(x ~ 1 + offset(log(lib.size))), silent = TRUE)
}


model_zip <- function(x, lib.size){
  library(pscl)
  try(zeroinfl(x ~ 1 | 1, offset=log(lib.size), dist = "poisson"), silent = TRUE)
}


model_zinb <- function(x, lib.size){
  library(pscl)
  try(zeroinfl(x ~ 1 | 1, offset=log(lib.size), dist = "negbin"), silent = TRUE)
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

result_KS <- stack(result_KS)
result_KS$values <- as.numeric(result_KS$values)
result_KS <- result_KS[!is.na(result_KS$values), ]

result_KSig <- result_KS[result_KS$values >0.01, ]
result_KSigenes <- result_KSig$ind

stim_data <- result[rownames(result) %in% result_KSigenes,]


#Convert the matrix to sparse matrix
sim <- Matrix(stim_data, sparse = TRUE)



sim_poi <- future_apply(sim, MARGIN = 1L, FUN = model_poi, lib.size <- lib.size, future.seed = 0xBEEF)
sim_nb <- future_apply(sim, MARGIN = 1L, FUN = model_nb, lib.size <- lib.size, future.seed = 0xBEEF)
sim_zip <- future_apply(sim, MARGIN = 1L, FUN = model_zip, lib.size <- lib.size, future.seed = 0xBEEF)
sim_zinb <- future_apply(sim, MARGIN = 1L, FUN = model_zinb, lib.size <- lib.size, future.seed = 0xBEEF)


sim_poi_BIC <- t(as.data.frame(lapply(sim_poi, model_BIC)))
sim_nb_BIC <- t(as.data.frame(lapply(sim_nb, model_BIC)))
sim_zip_BIC <- t(as.data.frame(lapply(sim_zip, model_BIC)))
sim_zinb_BIC <- t(as.data.frame(lapply(sim_zinb, model_BIC)))


result_BIC <- cbind(sim_poi_BIC, sim_nb_BIC, sim_zip_BIC, sim_zinb_BIC)
colnames(result_BIC) <- c("Poi_BIC", "NB_BIC", "ZIP_BIC", "ZINB_BIC")

save(result_BIC, file = "pbmc3k/n3000/result_BIC.rds")

#####################################################################################################
stim_data_list <- lapply(as.list(1:dim(stim_data)[1]), function(x) stim_data[x[1],])
names(stim_data_list) <- rownames(stim_data)

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

#Fit the models
#Poisson distribution
model_Chipoi <- function(x, lib.size){
  
  library(stats)
  
  1-pchisq(summary(glm(x ~ 1, offset=log(lib.size), family = "poisson"))$deviance, 
           df= summary(glm(x ~ 1, offset=log(lib.size), family = "poisson"))$df.residual)
}

#NB distribution
model_Chinb <- function(x, lib.size){
  
  library(MASS)
  
  1-pchisq(summary(glm.nb(x ~ 1 + offset(log(lib.size))))$deviance, 
           df= summary(glm.nb(x ~ 1 + offset(log(lib.size))))$df.residual)
}

#ZIP distribution
model_Chizip <- function(x, lib.size){
  library(pscl)
  library(emdbook)
  
  1-pchibarsq(2*(logLik(zeroinfl(x ~ 1 | 1, offset=log(lib.size), dist = "poisson"))-logLik(glm(x ~ 1, offset=log(lib.size), family = "poisson"))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)
  
}


#ZINB distribution
model_Chizinb <- function(x, lib.size){
  library(pscl)
  library(emdbook)
  
  1-pchibarsq(2*(logLik(zeroinfl(x ~ 1 | 1, offset=log(lib.size), dist = "negbin"))-logLik(glm.nb(x ~ 1 + offset(log(lib.size))))),df=1,mix=0.5,lower.tail=TRUE, log.p = FALSE)
}


BM_genes <- Best_model(result_BIC, stim_data_list)


Chipoi_pval <- future_lapply(BM_genes[[1]], model_Chipoi, lib.size <- lib.size, future.seed = 0xBEEF)
Chinb_pval <- future_lapply(BM_genes[[2]], model_Chinb, lib.size <- lib.size, future.seed = 0xBEEF)
Chizip_pval <- future_lapply(BM_genes[[3]], model_Chizip, lib.size <- lib.size, future.seed = 0xBEEF)
Chizinb_pval <- future_lapply(BM_genes[[4]], model_Chizinb, lib.size <- lib.size, future.seed = 0xBEEF)


result_pval <- list(Chipoi_pval, Chinb_pval, Chizip_pval, Chizinb_pval)


save(result_pval, file = "pbmc3k/n3000/result_pval.rds")



