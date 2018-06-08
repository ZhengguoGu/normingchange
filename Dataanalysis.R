##############################################################
#########                                        #############
######### Data Analysis: Norming change scores   #############
#########                                        #############
#########      Zhengguo Gu, Tilburg University   #############
#########      Last update: 1 June, 2018         #############
#########                                        #############
# Note: the descriptive analysis is at the bottom of this file 
##############################################################


####################### PART I: Calculating IPR ####################################

####################################################################################
###### NOTE: PART I was used to generate the dataset.                         ######
###### 20171011IPR.RData.RData                                                ######
####################################################################################
library(Kendall)
library(foreach)
library(doSNOW)
library(doRNG)
 
IPR_reg <- matrix(NA, 270, 9) 
IPR_Tscore <- matrix(NA, 270, 9)
Rankcorrelation_mean <- matrix(NA, 270, 2)
Rankcorrelation_sd <- matrix(NA, 270, 2)
Par_changescoreMean <- matrix(NA, 270, 2)  #first column: averaged reliability, second column: sd
Par_changescoreSD <- matrix(NA, 270, 2)


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
    datalist[[j]] <- sim_result[[k]]
    rellist[j, ] <- sim_result[[l]]
  }
    
  Par_changescoreMean[num_test,] <- colMeans(rellist)
  Par_changescoreSD[num_test,] <- apply(rellist, 2, sd)
  
  colnames(Par_changescoreMean) <- c("change_rel", "var_pre", "var_post", "cor_prepost", "cor_preD")
  colnames(Par_changescoreSD) <- c("change_rel", "var_pre", "var_post", "cor_prepost", "cor_preD")
  ##############################################################
  ###### norming methods
  ##############################################################
  # 3. The regression-based change approach and T-score

  changescores <- lapply(datalist, changescore)  
  #identical(changescores[[15]], datalist[[15]][,2] - datalist[[15]][,1])

  # Parallel computing
  cl <- makeCluster(4)
  registerDoSNOW(cl)

  set.seed(112)  # set seed, gonna use parallel computing
  ZT_result <- foreach(i = 1:1000, .combine='cbind') %dorng% {
    
    # regression-based change approach
    fit <- lm(changescores[[i]] ~ X1 + X2)
    Escore <- changescores[[i]] - predict(fit) #residual
    SD_e <- sqrt(sum(Escore^2)/(length(Escore) - 2))
    Zscore <- Escore/SD_e
    qZ <- quantile(Zscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    rank_cor_Z <- Kendall::Kendall(theta_D, Zscore)$tau
    
    # Tscore
    fit2 <- lm(datalist[[i]][, 2] ~ datalist[[i]][, 1] + X1 + X2)
    Escore2 <-  datalist[[i]][, 2] - predict(fit2)
    SD_e2 <- sqrt(sum(Escore2^2)/(length(Escore2) - 2))
    Tscore <- Escore2/SD_e2
    qT <- quantile(Tscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    rank_cor_T <- Kendall::Kendall(theta_D, Tscore)$tau

    perc <- list()
    perc[[1]] <- rbind(qZ, qT)
    perc[[2]] <- c(rank_cor_Z, rank_cor_T)
    
    return(perc)
  }
  stopCluster(cl)

  qZmatrix <- matrix(NA, 1000, 9)
  qTmatrix <- matrix(NA, 1000, 9)
  rank_corMat <- matrix(NA, 1000, 2)
  
  for(i in 0:999){
    
    j <- 2*i+1
    k <- 2*i+2
    l <- i + 1
    qZmatrix[l, ] <- ZT_result[[j]][1,]
    qTmatrix[l, ] <- ZT_result[[j]][2,]
    rank_corMat[l, ] <- ZT_result[[k]]
    
  }
 
  
  IPR_reg[num_test, ] <- apply(qZmatrix, 2, calculate_IPR)
  IPR_Tscore[num_test, ] <- apply(qTmatrix, 2, calculate_IPR)
  Rankcorrelation_mean[num_test, ] <- apply(rank_corMat, 2, mean)
  Rankcorrelation_sd[num_test, ] <- apply(rank_corMat, 2, sd)
  
  num_test <- num_test + 1 
}
colnames(IPR_reg) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")
colnames(IPR_Tscore) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")
colnames(Rankcorrelation_mean) <- c("regression based", "T score")
colnames(Rankcorrelation_sd) <- c("regression based", "T score")
save(IPR_reg, IPR_Tscore, Rankcorrelation_mean, Rankcorrelation_sd, Par_changescoreMean, Par_changescoreSD, file = "D:/ZG/20171011IPR.RData")


########################### PART II: ANOVA ###########################################################

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

# 2. Combine results with design factors.
load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/0submissionAssessment/1ReviseResubmit/20180604Simulation/smallVarRho0_IPR.RData") #Load simulation results

IPR_regData <- cbind(IPR_reg, df, "regression-based", "VarD=.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData <- cbind(IPR_Tscore, df, "TScore", "VarD=.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                        "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                        "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                        "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                        "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                        "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                        "Kendall (mean)", "Kendall (sd)")

load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/0submissionAssessment/1ReviseResubmit/20180604Simulation/smallVarRho1_IPR.RData") #Load simulation results
IPR_regData_temp1 <- cbind(IPR_reg, df, "regression-based", "VarD=.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp1 <- cbind(IPR_Tscore, df, "TScore", "VarD=.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp1)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                        "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                        "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                        "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp1)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                      "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                      "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                      "Kendall (mean)", "Kendall (sd)")


load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/0submissionAssessment/1ReviseResubmit/20180604Simulation/smallVarRho1Neg_IPR.RData")
IPR_regData_temp2 <- cbind(IPR_reg, df, "regression-based", "VarD=.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp2 <- cbind(IPR_Tscore, df, "TScore", "VarD=.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp2)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp2)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")


load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/0submissionAssessment/1ReviseResubmit/20180604Simulation/largeVarRho0_IPR.RData")
IPR_regData_temp3 <- cbind(IPR_reg, df, "regression-based", "VarD=1.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp3 <- cbind(IPR_Tscore, df, "TScore", "VarD=1.14", "Rho_preD=0", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp3)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp3)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")

load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/0submissionAssessment/1ReviseResubmit/20180604Simulation/largeVarRho1_IPR.RData")
IPR_regData_temp4 <- cbind(IPR_reg, df, "regression-based", "VarD=1.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp4 <- cbind(IPR_Tscore, df, "TScore", "VarD=1.14", "Rho_preD=0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp4)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp4)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")

load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/0submissionAssessment/1ReviseResubmit/20180604Simulation/largeVarRho1Neg_IPR.RData")
IPR_regData_temp5 <- cbind(IPR_reg, df, "regression-based", "VarD=1.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,1], Rankcorrelation_sd[,1])
IPR_TData_temp5 <- cbind(IPR_Tscore, df, "TScore", "VarD=1.14", "Rho_preD=-0.1", Par_changescoreMean, Par_changescoreSD, Rankcorrelation_mean[,2], Rankcorrelation_sd[,2])
colnames(IPR_regData_temp5)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                              "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                              "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                              "Kendall (mean)", "Kendall (sd)")
colnames(IPR_TData_temp5)[seq(14, 28)] <- c("norming method", "var true change", "Rho_preD",
                                            "Rel_change (mean)", "var_pre (mean)", "var_post (mean)", "cor_prepost (mean)", "cor_preD (mean)",
                                            "Rel_change (mean)", "var_pre (sd)", "var_post (sd)", "cor_prepost (sd)", "cor_preD (sd)",
                                            "Kendall (mean)", "Kendall (sd)")

IPR_regF <- rbind(IPR_regData, IPR_regData_temp1, IPR_regData_temp2, IPR_regData_temp3, IPR_regData_temp4, IPR_regData_temp5)
IPR_TDataF <- rbind(IPR_TData, IPR_TData_temp1, IPR_TData_temp2, IPR_TData_temp3, IPR_TData_temp4, IPR_TData_temp5)

IPR_regF <- cbind(seq(1620), IPR_regF)
IPR_TDataF <- cbind(seq(1620), IPR_TDataF)

IPR_FINALData <- rbind(IPR_regF, IPR_TDataF)  #important! 
colnames(IPR_FINALData)[1] <- "Cell No."

save(IPR_FINALData, file = "D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/0submissionAssessment/1ReviseResubmit/20180604Simulation/IPR_FINALData.RData")

###########################################
# 3. mixed anova 
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

fit1 <- aov(`1%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
                  sample_size : proportionExplained + 
                  sample_size : polytomous +
                  sample_size : num_items +
                  sample_size : NormingMethod +
                  Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit1)

eta1 <- EtaSq(fit1, type = 1)


fit5 <- aov(`5%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
              sample_size : proportionExplained + 
              sample_size : polytomous +
              sample_size : num_items +
              sample_size : NormingMethod +
              Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit5)
eta5 <- EtaSq(fit5, type = 1)

fit10 <- aov(`10%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
              sample_size : proportionExplained + 
              sample_size : polytomous +
              sample_size : num_items +
              sample_size : NormingMethod +
               Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit10)
eta10 <- EtaSq(fit10, type = 1)

fit25 <- aov(`25%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : NormingMethod +
               Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit25)
eta25 <- EtaSq(fit25, type = 1)

fit50 <- aov(`50%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : NormingMethod +
               Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit50)
eta50 <- EtaSq(fit50, type = 1)


fit75 <- aov(`75%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : NormingMethod +
               Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit75)
eta75 <- EtaSq(fit75, type = 1)

fit90 <- aov(`90%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : NormingMethod +
               Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit90)
eta90 <- EtaSq(fit90, type = 1)

fit95 <- aov(`95%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : NormingMethod +
               Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit95)
eta95 <- EtaSq(fit95, type = 1)

fit99 <- aov(`99%` ~ sample_size + proportionExplained + polytomous + num_items + NormingMethod+
               sample_size : proportionExplained + 
               sample_size : polytomous +
               sample_size : num_items +
               sample_size : NormingMethod +
               Error(simulationcell/NormingMethod), data=IPR_Data)
summary(fit99)
eta99 <- EtaSq(fit99, type = 1)

etamatrix <- cbind(eta1[,2], eta5[, 2], eta10[, 2], eta25[, 2], eta50[, 2], eta75[, 2], eta90[, 2], eta95[, 2], eta99[, 2])
colnames(etamatrix) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")
write.table(etamatrix, file = 'D:/ZG/etamatrix.txt', sep = ',')  #save at TiSEM blade server
write.table(etamatrix, file = 'D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\20171010 newdata/etamatrix.txt', sep = ',')  #save at home computer


# 4. relationship between sample size and IPRs
IPR_Data # make use of this datamatrix (obtained from part II, see above)

aggreIPR <- aggregate(IPR_Data[, 1:9], by = list(IPR_Data$'sample_size', IPR_Data$'NormingMethod'), FUN = mean) 

#plot

maintext = paste(c("1th", "5th", "10th", "25th", "50th", "75th", "90th", "95th", "99th"), "percentile")

for(i in 1:9){
  
  filename <- paste("D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\20171010 newdata/IPR_", i, ".png", sep = "")
  png(file=filename, width = 1200, height = 1200, units = "px")
  
  layout(rbind(1,2), heights=c(10,1))# put legend on bottom 1/10th of the chart (note, this is from http://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot)
  par(mar=c(5,4,4,5)+.1)
  plot(aggreIPR[1:15, i+2], type = "b", xlab = "Sample size", ylab = "IPR", xaxt="n", main = maintext[i], ylim = c(0, 1), pch = 15)
  axis(1, at=1:15,labels=aggreIPR[1:15, 1], las=2)
  lines(aggreIPR[16:30, i+2], pch = 18, type = "b", lty = 2)
  
  par(mar=c(0,0,0,0))
  plot.new()
  legend("center", "groups",
         c("regression-based change approach", "Tscore method"),
         pch=c(15, 18),
         lty = c(1, 2),
         ncol=2, bty = "n")
  dev.off()
}


################ PART III: Descriptives ###################
# Note: Ideally descriptive table should be provided before the analysis
# but I forgot to provide the table until Wilco read my draft and asked about it.
# Here I include the descriptive table. 
# Because sample sizes differ, I choose sample size = 1000. 


sample_sizeV = c(100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 9000, 10000)  #sample size 
propEV <- c(0, .065/2, .13/2, .26/2)  # proportion of explained by each predictor. 
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

TESTSCORE <- matrix(NA, 24, 6)
colnames(TESTSCORE) <- c("mean_pretest", "sd_pretest", "mean_posttest", "sd_posttest", "mean_change", "sd_change")
num_test <- 49
while(num_test <= 72){
  
  # 1. reorganize the data.

  #filename <- paste("D:/ZhengguoProj2Blad analysis/results_", num_test, ".RData", sep = "")
  filename <- paste("D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\20171010 newdata/results_", num_test, ".RData", sep = "")
  
  # note that "D:/ZhengguoProj2Blad analysis/results..." means that the analysis is done on Blade Server
  load(filename)
  
  #dim(sim_result)  #1000 (persons) x 3000 matrix
  
  # 2. reorganize the data 
  pretest <- array()
  posttest <- array()
  changescore <- array()
  j <- 1
  for(i in seq(from = 1, to = 3000, by = 3)){
    pretest[j] <- mean(sim_result[, i])
    posttest[j] <- mean(sim_result[, i+1])
    changescore[j] <- posttest[j] - pretest[j]
    j <- j + 1
  }
  
  TESTSCORE[num_test-48,1] <- mean(pretest) 
  TESTSCORE[num_test-48,2] <- sd(pretest)
  TESTSCORE[num_test-48,3] <- mean(posttest)
  TESTSCORE[num_test-48,4] <- sd(posttest)
  TESTSCORE[num_test-48,5] <- mean(changescore)
  TESTSCORE[num_test-48,6] <- sd(changescore)

  num_test <- num_test + 1
}

  # 3. Combine with the data we have in PART III
load(file ="D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\20171010 newdata/20171011IPR.RData")
DescripTable <- cbind(df[49:72, ], TESTSCORE, REL_changescore[49:72,])
colnames(DescripTable)[c(1:4, 11:12)] <- c("sample_size", "proportionExplained", "polytomous", "num_items", "change_reliability", "sd_reliability")

save(DescripTable, file ="D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\20171010 newdata/DescriTable.RData")
write.table(DescripTable, file = 'D:\\Dropbox\\Tilburg office\\Research Individual change\\Project 2 - norming change\\20171010 newdata/DescriTable.txt', sep = ',')
