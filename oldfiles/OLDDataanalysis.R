##############################################################
#########                                        #############
######### Data Analysis: Norming change scores   #############
#########                                        #############
#########                                        #############
#########                                        #############
# Note: the descriptive analysis is at the bottom of this file 
##############################################################


####################### PART I: Calculating IPR and reorganizing the data for further analysis  ####################################

library(Kendall)
library(foreach)
library(doSNOW)
library(doRNG)
library(psych)
 
IPR_reg <- matrix(NA, 270, 9) 
IPR_Tscore <- matrix(NA, 270, 9)
Rankcorrelation_mean <- matrix(NA, 270, 2)
Rankcorrelation_sd <- matrix(NA, 270, 2)
Par_changescoreMean <- matrix(NA, 270, 5)  #"change_rel", "var_pre", "var_post", "cor_prepost", "cor_preD"
Par_changescoreSD <- matrix(NA, 270, 5)
colnames(Par_changescoreMean) <- c("change_rel", "var_pre", "var_post", "cor_prepost", "cor_preD")
colnames(Par_changescoreSD) <- c("change_rel", "var_pre", "var_post", "cor_prepost", "cor_preD")
reg_coef_mean <- matrix(NA, 270, 5)   # regression based: coefficent X1, and X2; Tscore: pre, X1, X2
reg_coef_max <- matrix(NA, 270, 5)
reg_coef_min <- matrix(NA, 270, 5)
colnames(reg_coef_mean) <- c("reg-based: X1", "reg-based: X2", "Tscore: pre", "Tscore: X1", "Tscore: X2")
colnames(reg_coef_max) <- c("reg-based: X1", "reg-based: X2", "Tscore: pre", "Tscore: X1", "Tscore: X2")
colnames(reg_coef_min) <- c("reg-based: X1", "reg-based: X2", "Tscore: pre", "Tscore: X1", "Tscore: X2")
####### two functions ###############
changescore <- function(Data){
  change <- Data[, 2] - Data[, 1]
  return(change)
}

# interprecentile range
calculate_IPR <- function(DATA){
  result <- quantile(DATA, c(.025, .975))
  IPR <- result[2] - result[1]
}

########################################

num_test <- 1
while(num_test <= 270){
  print(num_test)
  
  # 1. reorganize the data.
  
  filename <- paste("results_", num_test, ".RData", sep = "")  
     
  load(filename)

  ## 2. standardize scores. 
  #mydata <- scale(sim_result, center = TRUE, scale = TRUE)  #Rescale pretest and posttest raw scores. note that sim_result contains estimated reliability, but scaling does not affect reliability.

  # 3. reorganize the data 
  datalist <- list()
  rellist <- matrix(NA, 1000, 5)
  
  for(i in 0:999){
    j <- i+1
    k <- 2*i+1
    l <- 2*i+2
    datalist[[j]] <- sim_result[[k]]  #contains the simulated pretest and posttest scores
    rellist[j, ] <- sim_result[[l]]   #contains the change-score reliability, variance of pretest scores, variance of posttest, correlation betwen pretest posttest, correlation between pretest and change.
  }
    
  Par_changescoreMean[num_test,] <- colMeans(rellist)
  Par_changescoreSD[num_test,] <- apply(rellist, 2, sd)
  
  
  ##############################################################
  ###### norming methods
  ##############################################################
  # 3. The regression-based change approach and T-score

  changescores <- lapply(datalist, changescore)  
  #identical(changescores[[15]], datalist[[15]][,2] - datalist[[15]][,1])

  # Parallel computing
  cl <- makeCluster(10)
  registerDoSNOW(cl)

  set.seed(112)  # set seed, gonna use parallel computing
  ZT_result <- foreach(i = 1:1000, .combine='cbind') %dorng% {
    
    # regression-based change approach
    fit <- lm(changescores[[i]] ~ X1 + X2)
    coef_regbased <- summary(fit)$coefficients[c(2,3), 1]
    Escore <- changescores[[i]] - predict(fit) #residual
    SD_e <- sqrt(sum(Escore^2)/(length(Escore) - 2))
    Zscore <- Escore/SD_e
    qZ <- quantile(Zscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    rank_cor_Z <- Kendall::Kendall(theta_D, Zscore)$tau
    
    # Tscore
    fit2 <- lm(datalist[[i]][, 2] ~ datalist[[i]][, 1] + X1 + X2)
    coef_Tscore <- summary(fit2)$coefficients[c(2, 3, 4), 1]
    Escore2 <-  datalist[[i]][, 2] - predict(fit2)
    SD_e2 <- sqrt(sum(Escore2^2)/(length(Escore2) - 2))
    Tscore <- Escore2/SD_e2
    qT <- quantile(Tscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    rank_cor_T <- Kendall::Kendall(res_trueTheta, Tscore)$tau  # note: there the res_trueTheta is the residual true change score (theta)

    perc <- list()
    perc[[1]] <- rbind(qZ, qT)
    perc[[2]] <- c(rank_cor_Z, rank_cor_T)
    perc[[3]] <- c(coef_regbased, coef_Tscore)
    return(perc)
  }
  stopCluster(cl)

  qZmatrix <- matrix(NA, 1000, 9)
  qTmatrix <- matrix(NA, 1000, 9)
  rank_corMat <- matrix(NA, 1000, 2)
  coeff_ZT <- matrix(NA, 1000, 5)  #regression coefficients
  for(i in 0:999){
    
    j <- 3*i+1
    k <- 3*i+2
    m <- 3*i+3
    l <- i + 1
    qZmatrix[l, ] <- ZT_result[[j]][1,]
    qTmatrix[l, ] <- ZT_result[[j]][2,]
    rank_corMat[l, ] <- ZT_result[[k]]
    coeff_ZT[l, ] <- ZT_result[[m]]
    
  }
 
  
  IPR_reg[num_test, ] <- apply(qZmatrix, 2, calculate_IPR)
  IPR_Tscore[num_test, ] <- apply(qTmatrix, 2, calculate_IPR)
  Rankcorrelation_mean[num_test, ] <- apply(rank_corMat, 2, mean)
  Rankcorrelation_sd[num_test, ] <- apply(rank_corMat, 2, sd)
  reg_coef_mean[num_test, ] <- apply(coeff_ZT, 2, mean) # regression based: coefficent X1, and X2; Tscore: pre, X1, X2
  reg_coef_max[num_test, ] <- apply(coeff_ZT, 2, max)
  reg_coef_min[num_test, ] <- apply(coeff_ZT, 2, min)
  
  
  num_test <- num_test + 1 
}
colnames(IPR_reg) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")
colnames(IPR_Tscore) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")
colnames(Rankcorrelation_mean) <- c("regression based", "T score")
colnames(Rankcorrelation_sd) <- c("regression based", "T score")
save(IPR_reg, IPR_Tscore, Rankcorrelation_mean, Rankcorrelation_sd, Par_changescoreMean, Par_changescoreSD, reg_coef_mean, reg_coef_max, reg_coef_min, file = "??.RData")  #save data


########################### PART II: Comparing rank correlations and ANOVAs ###########################################################



# 1. get the design factors.
# Note: The following code is from simulationsNEW.R, where we know which design factors are in which cell.
sample_sizeV = c(100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 9000, 10000)  #sample size 
propEV <- c("R2=.065", "R2=.13", "R2=.26")  # proportion of explained by each predictor. Note that I relabelled them (compared to the simulation)
polytomousV <- c(TRUE, FALSE)  # if true, simulate polytomous response data 
num_itemsV <- c(10, 20, 40) # 10, 20, and 40
num_conditions <- length(sample_sizeV) * length(propEV) * length(polytomousV) * length(num_itemsV)
conditions <- list()
p <- 1
for(i in 1:length(sample_sizeV)){
  for(j in 1:length(propEV)){
    for(k in 1:length(polytomousV)){
      for(l in 1:length(num_itemsV)){
        conditions[[p]] <- c(sample_sizeV[i], propEV[j], polytomousV[k], num_itemsV[l])
        p <- p+1
      }
    }
  }
}

df <- data.frame(matrix(unlist(conditions), nrow=num_conditions, byrow = T))
colnames(df) <- c("sample_size", "proportionExplained", "polytomous", "num_items")

# 2. Combine results with design factors. Note that there are in total 6 data files from simulationsNEW.R

#-- file 1
load(file ="varD_14_rho_0.RData") #Load simulation results

IPR_regData <- cbind(IPR_reg, df, "regression-based", "VarD=.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData <- cbind(IPR_Tscore, df, "TScore", "VarD=.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                        "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                        "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                        "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                        "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                        "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                        "Kendall (mean)", "Kendall (sd)")
coef_data1 <- cbind(df, reg_coef_mean)
write.csv(coef_data1, file = "estimatedCoef_varD_14_rho_0.csv")  # these coefficients will be submitted as supplementary material

#-- file 2
load(file ="varD_14_rho_1.RData") #Load simulation results
IPR_regData_temp1 <- cbind(IPR_reg, df, "regression-based", "VarD=.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp1 <- cbind(IPR_Tscore, df, "TScore", "VarD=.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp1)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                        "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                        "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                        "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp1)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                      "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                      "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                      "Kendall (mean)", "Kendall (sd)")
coef_data2 <- cbind(df, reg_coef_mean)
write.csv(coef_data2, file = "estimatedCoef_varD_14_rho_1.csv") # these coefficients will be submitted as supplementary material

#-- file 3
load(file ="varD_14_rho_neg1.RData")
IPR_regData_temp2 <- cbind(IPR_reg, df, "regression-based", "VarD=.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp2 <- cbind(IPR_Tscore, df, "TScore", "VarD=.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp2)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp2)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")
coef_data3 <- cbind(df, reg_coef_mean)
write.csv(coef_data3, file = "estimatedCoef_varD_14_rho_neg1.csv") # these coefficients will be submitted as supplementary material

#-- file 4
load(file ="varD_114_rho_0.RData")
IPR_regData_temp3 <- cbind(IPR_reg, df, "regression-based", "VarD=1.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp3 <- cbind(IPR_Tscore, df, "TScore", "VarD=1.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp3)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp3)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")
coef_data4 <- cbind(df, reg_coef_mean)
write.csv(coef_data4, file = "estimatedCoef_varD_114_rho_0.csv") # these coefficients will be submitted as supplementary material

#-- file 5
load(file ="varD_114_rho_1.RData")
IPR_regData_temp4 <- cbind(IPR_reg, df, "regression-based", "VarD=1.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp4 <- cbind(IPR_Tscore, df, "TScore", "VarD=1.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp4)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp4)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")
coef_data5 <- cbind(df, reg_coef_mean)
write.csv(coef_data5, file = "estimatedCoef_varD_114_rho_1.csv") # these coefficients will be submitted as supplementary material

#-- file 6
load(file ="varD_114_rho_neg1.RData")
IPR_regData_temp5 <- cbind(IPR_reg, df, "regression-based", "VarD=1.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp5 <- cbind(IPR_Tscore, df, "TScore", "VarD=1.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp5)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp5)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (sd)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")
write.csv(coef_data6, file = "estimatedCoef_varD_114_rho_neg1.csv") # these coefficients will be submitted as supplementary material


IPR_regF <- rbind(IPR_regData, IPR_regData_temp1, IPR_regData_temp2, IPR_regData_temp3, IPR_regData_temp4, IPR_regData_temp5)
IPR_TDataF <- rbind(IPR_TData, IPR_TData_temp1, IPR_TData_temp2, IPR_TData_temp3, IPR_TData_temp4, IPR_TData_temp5)

IPR_regF <- cbind(seq(1620), IPR_regF)
IPR_TDataF <- cbind(seq(1620), IPR_TDataF)

IPR_FINALData <- rbind(IPR_regF, IPR_TDataF)  #important! 
colnames(IPR_FINALData)[1] <- "Cell No."

save(IPR_FINALData, file = "IPR_FINALData.RData")

###########################################
#### 4. Comparing rank correlations
mean(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based", 28])  #=0.3576953
mean(IPR_FINALData[IPR_FINALData$`norming method`=="TScore", 28])  #=0.4294376


### 4.1 a few descriptives 
# 4.1.1 var(D)=.14, 0 correlation between pretest and change (theta level)
summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel

load(file ="varD_14_rho_0.RData")
apply(reg_coef_mean, 2, mean)  # here we check the estimated regression coefficents, and summerized them in a table.

# 4.1.2 var(D)=.14, .1 correlation between pretest and change
summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel

# 4.1.3 var(D)=.14, -0.1 correlation between pretest and change
summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel

# 4.1.4 var(D)=1.14, 0 correlation between pretest and change
summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel

# 4.1.5 var(D)=1.14, 0.1 correlation between pretest and change 
summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel

# 4.1.6 var(D)=1.14, -0.1 correlation between pretest and change
summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 10, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 20, ])  #note that the results are copied to excel

summary(IPR_FINALData[IPR_FINALData$`norming method`=="regression-based" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel
summary(IPR_FINALData[IPR_FINALData$`norming method`=="TScore" & IPR_FINALData$`var true change`=="VarD=1.14" & IPR_FINALData$Rho_preD=="Rho_preD=-0.1" & IPR_FINALData$num_items == 40, ])  #note that the results are copied to excel


# 5. mixed anova 
library(DescTools)  #to calculate eta squared

summary(IPR_FINALData)  #design factors need to be coded as factors (categorical)
IPR_FINALData$num_items <- factor(IPR_FINALData$num_items, levels = c(10, 20, 40))
IPR_FINALData$polytomous <- factor(IPR_FINALData$polytomous, levels = c("TRUE", "FALSE"))
IPR_FINALData$proportionExplained <- factor(IPR_FINALData$proportionExplained, levels = c("R2=.065", "R2=.13", "R2=.26"))
IPR_FINALData$`Cell No.` <- factor(IPR_FINALData$`Cell No.`, levels = seq(1:1620))
IPR_FINALData$`norming method` <- factor(IPR_FINALData$`norming method`, levels = c("regression-based", "TScore"))
IPR_FINALData$`var true change` <- factor(IPR_FINALData$`var true change`, levels = c("VarD=.14", "VarD=1.14"))
IPR_FINALData$Rho_preD <- factor(IPR_FINALData$Rho_preD, levels = c("Rho_preD=0", "Rho_preD=0.1", "Rho_preD=-0.1"))
summary(IPR_FINALData) #check again

fit1 <- aov(`1%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
                  sample_size : proportionExplained + 
                  sample_size : polytomous +
                  sample_size : num_items +
                  sample_size : `norming method` +
                  sample_size : `var true change` +
                  sample_size : Rho_preD +
                  Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit1)

eta1 <- EtaSq(fit1, type = 1)


fit5 <- aov(`5%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
              sample_size : proportionExplained + 
              sample_size : polytomous +
              sample_size : num_items +
              sample_size : `norming method` +
              sample_size : `var true change` +
              sample_size : Rho_preD +
              Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit5)
eta5 <- EtaSq(fit5, type = 1)

fit10 <- aov(`10%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
              sample_size : proportionExplained + 
              sample_size : polytomous +
              sample_size : num_items +
              sample_size : `norming method` + 
              sample_size : `var true change` +
              sample_size : Rho_preD +
              Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit10)
eta10 <- EtaSq(fit10, type = 1)

fit25 <- aov(`25%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : `norming method` +
               sample_size : `var true change` +
               sample_size : Rho_preD +
               Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit25)
eta25 <- EtaSq(fit25, type = 1)

fit50 <- aov(`50%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : `norming method` +
               sample_size : `var true change` +
               sample_size : Rho_preD +
               Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit50)
eta50 <- EtaSq(fit50, type = 1)


fit75 <- aov(`75%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : `norming method` +
               sample_size : `var true change` +
               sample_size : Rho_preD +
               Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit75)
eta75 <- EtaSq(fit75, type = 1)

fit90 <- aov(`90%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : `norming method` +
               sample_size : `var true change` +
               sample_size : Rho_preD +
               Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit90)
eta90 <- EtaSq(fit90, type = 1)

fit95 <- aov(`95%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : `norming method` +
               sample_size : `var true change` +
               sample_size : Rho_preD +
               Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit95)
eta95 <- EtaSq(fit95, type = 1)

fit99 <- aov(`99%` ~ sample_size + proportionExplained + polytomous + num_items + `norming method` + `var true change` + Rho_preD +
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : `norming method` +
               sample_size : `var true change` +
               sample_size : Rho_preD +
               Error(`Cell No.`/`norming method`), data=IPR_FINALData)
summary(fit99)
eta99 <- EtaSq(fit99, type = 1)

etamatrix <- cbind(eta1[,2], eta5[, 2], eta10[, 2], eta25[, 2], eta50[, 2], eta75[, 2], eta90[, 2], eta95[, 2], eta99[, 2])
colnames(etamatrix) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")
write.table(etamatrix, file = 'etamatrix.txt', sep = ',')



# 6. relationship between sample size and IPRs
# var(theta_D) =  .14
IPR_smallvar <- IPR_FINALData[IPR_FINALData$`var true change`=="VarD=.14", ]

IPR_smallvar <- aggregate(IPR_smallvar[ , 2:10], by = list(IPR_smallvar$'sample_size', IPR_smallvar$`norming method`), FUN = mean) 
IPR_smallvar[, 1] <- as.numeric(levels(IPR_smallvar[, 1]))[IPR_smallvar[, 1]]
IPR_smallvar  <- IPR_smallvar[order(IPR_smallvar[, 1]), ]


#plot
m <- matrix(c(1,2,3,4,5,6,7, 8, 9, 10, 10, 10),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(.3, .3, .3, .1))


maintext = paste(c("1th", "5th", "10th", "25th", "50th", "75th", "90th", "95th", "99th"), "percentile")

for(i in 1:9){
  
  #filename <- paste("D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\0submissionAssessment\\1ReviseResubmit\\IPR_", i, ".png", sep = "")
  #png(file=filename, width = 1200, height = 1200, units = "px")
  
  #layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
  par(mar=c(3,3,3,3))
  plot(IPR_smallvar[IPR_smallvar[, 2] == "regression-based", i+2], type = "b", xlab = "Sample size", ylab = "IPR", xaxt="n", main = maintext[i], ylim = c(0, 1), pch = 15)
  axis(1, at=1:15,labels=IPR_smallvar[IPR_smallvar[, 2] == "regression-based", 1], las=2)
  lines(IPR_smallvar[IPR_smallvar[, 2] == "TScore", i+2], pch = 18, type = "b", lty = 2)
  
  #par(mar=c(0,0,0,0))
 # plot.new()
  #legend("center", "groups",
  #       c("regression-based change approach", "Tscore method"),
   #      pch=c(15, 18),
  #       lty = c(1, 2),
  #       ncol=2, bty = "n")
  #dev.off()
}

par(mar=c(0,0,0,0))
plot.new()
legend("center", "groups",
        c("regression-based change approach", "Tscore method"),
        pch=c(15, 18),
        lty = c(1, 2),
        ncol=2, bty = "n", 
        cex = 1.75)


# var(theta_D) =  1.14
IPR_largevar <- IPR_FINALData[IPR_FINALData$`var true change`=="VarD=1.14", ]

IPR_largevar <- aggregate(IPR_largevar[ , 2:10], by = list(IPR_largevar$'sample_size', IPR_largevar$`norming method`), FUN = mean) 
IPR_largevar[, 1] <- as.numeric(levels(IPR_largevar[, 1]))[IPR_largevar[, 1]]
IPR_largevar  <- IPR_largevar[order(IPR_largevar[, 1]), ]


#plot
m <- matrix(c(1,2,3,4,5,6,7, 8, 9, 10, 10, 10),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(.3, .3, .3, .1))


maintext = paste(c("1th", "5th", "10th", "25th", "50th", "75th", "90th", "95th", "99th"), "percentile")

for(i in 1:9){
  
  #filename <- paste("D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\0submissionAssessment\\1ReviseResubmit\\IPR_", i, ".png", sep = "")
  #png(file=filename, width = 1200, height = 1200, units = "px")
  
  #layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
  par(mar=c(3,3,3,3))
  plot(IPR_largevar[IPR_largevar[, 2] == "regression-based", i+2], type = "b", xlab = "Sample size", ylab = "IPR", xaxt="n", main = maintext[i], ylim = c(0, 1), pch = 15)
  axis(1, at=1:15,labels=IPR_largevar[IPR_largevar[, 2] == "regression-based", 1], las=2)
  lines(IPR_largevar[IPR_largevar[, 2] == "TScore", i+2], pch = 18, type = "b", lty = 2)
  
  #par(mar=c(0,0,0,0))
  # plot.new()
  #legend("center", "groups",
  #       c("regression-based change approach", "Tscore method"),
  #      pch=c(15, 18),
  #       lty = c(1, 2),
  #       ncol=2, bty = "n")
  #dev.off()
}

par(mar=c(0,0,0,0))
plot.new()
legend("center", "groups",
       c("regression-based change approach", "Tscore method"),
       pch=c(15, 18),
       lty = c(1, 2),
       ncol=2, bty = "n", 
       cex = 1.75)


################  # ANOVAs for IPR generated by Equation (4a) and T-scores  
Data_4ANOVA <- FINALmat[, c(1:6, 12:20)]

Data_4ANOVA$var_D <- factor(Data_4ANOVA$var_D, levels = c(0.14, 1.14))
Data_4ANOVA$rho_preD <- factor(Data_4ANOVA$rho_preD, levels = c(-0.1, 0, 0.1))
#Data_4ANOVA$sample_size <- factor(Data_4ANOVA$sample_size, levels = seq(100, 1500, by = 100))
Data_4ANOVA$proportionExplained <- factor(Data_4ANOVA$proportionExplained, levels = c(0.065, 0.13, 0.26))
Data_4ANOVA$polytomous <- factor(Data_4ANOVA$polytomous, levels = c(0, 1))
Data_4ANOVA$num_items <- factor(Data_4ANOVA$num_items, levels = c(10, 20, 40))


model1 = ols(log(`qTZ1%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result1 <- robcov(model1)
qqnorm(Data_4ANOVA$`qTZ1%` - result1$residuals)  #plot is quite ok

model5 = ols(log(`qTZ5%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result5 <- robcov(model5)
qqnorm(Data_4ANOVA$`qTZ5%` - result5$residuals)  #normality seems to be ok

model10 = ols(log(`qTZ10%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result10 <- robcov(model10)
qqnorm(Data_4ANOVA$`qTZ10%` - result10$residuals) # plot is quite ok

model25 = ols(log(`qTZ25%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result25 <- robcov(model25)
qqnorm(Data_4ANOVA$`qTZ25%` - result25$residuals) # plot is quite ok

model50 = ols(log(`qTZ50%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result50 <- robcov(model50)
qqnorm(Data_4ANOVA$`qTZ50%` - result50$residuals) # plot is quite good

model75 = ols(log(`qTZ75%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result75 <- robcov(model75)
qqnorm(Data_4ANOVA$`qTZ75%` - result75$residuals) # plot is quite good

model90 = ols(log(`qTZ90%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result90 <- robcov(model90)
qqnorm(Data_4ANOVA$`qTZ90%` - result90$residuals) # plot is quite good

model95 = ols(log(`qTZ95%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result95 <- robcov(model95)
qqnorm(Data_4ANOVA$`qTZ95%` - result95$residuals) # plot is quite good

model99 = ols(log(`qTZ99%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result99 <- robcov(model99)
qqnorm(Data_4ANOVA$`qTZ99%` - result99$residuals) # plot is quite good

coef_table <- round(cbind(result1$coefficients, result5$coefficients, result10$coefficients, result25$coefficients,
                          result50$coefficients, result75$coefficients, result90$coefficients, result95$coefficients,
                          result99$coefficients), digits = 2)
