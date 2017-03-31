 ##############################################################
#########                                        #############
######### Data Analysis: Norming change scores   #############
#########                                        #############
#########      Zhengguo Gu, Tilburg University   #############
#########      Last update: 31/03/2017           #############
##############################################################

library(foreach)
library(doSNOW)
library(doRNG)
 
IPR_reg <- matrix(NA, 360, 9) 
IPR_Tscore <- matrix(NA, 360, 9)
num_test <- 1
while(num_test <= 360){
  
  # 1. reorganize the data.
  # To save space, for each test condition (in total 360), pretest and posttest
  # scores of 1000 datasets are 'cbind'ed. For example, the results of the first 
  # test condition:
  
  filename <- paste("D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170330 simresult/phdpro2BladZhengguo/results_", num_test, ".RData", sep = "")
  load(filename)

  #dim(sim_result)  #100 x 2000 matrix. 100: 100 persons. 2000: 1000 datasets, each contain a pretest and posttest
                 #score, and thus 1000x2=2000

  # 2. standardize scores. 
  mydata <- scale(sim_result, center = TRUE, scale = TRUE)

  # identical((sim_result[, 22] - mean(sim_result[, 22]))/sd(sim_result[, 22]), mydata[,22])

  # 3. reorganize the data 
  datalist <- list()
  j <- 1
  for(i in seq(from = 1, to = 2000, by = 2)){
    datalist[[j]] <- mydata[, c(i, i+1)]
    j <- j + 1
  }

  ##############################################################
  ###### norming methods
  ##############################################################
  # 3. The regression-based change approach and T-score

  changescore <- function(Data){
    change <- Data[, 2] - Data[, 1]
    return(change)
  }
  changescores <- lapply(datalist, changescore)
  # identical(trial[[15]], datalist[[15]][,2] - datalist[[15]][,1])


  # Parallel computing
  cl <- makeCluster(12)
  registerDoSNOW(cl)

  set.seed(112)  # set seed, gonna use parallel computing
  ZT_result <- foreach(i = 1:1000, .combine='rbind') %dorng% {
    
    # regression-based change approach
    fit <- lm(changescores[[i]] ~ X1 + X2)
    Escore <- predict(fit) - changescores[[i]]
    SD_e <- sqrt(sum(Escore^2)/(length(Escore) -2))
    Zscore <- Escore/SD_e
    qZ <- quantile(Zscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    # Tscore
    fit2 <- lm(datalist[[i]][, 2] ~ datalist[[i]][, 1] + X1 + X2)
    Escore2 <- predict(fit2) - datalist[[i]][, 2]
    SD_e2 <- sqrt(sum(Escore2^2)/(length(Escore2) -2))
    Tscore <- Escore2/SD_e2
    qT <- quantile(Tscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    
    perc <- rbind(qZ, qT)

    return(perc)
  }
  stopCluster(cl)


  qZmatrix <- ZT_result[seq(from = 1, to = dim(ZT_result)[1], by = 2), ] #1000 rows (i.e., 1000 samples), each row contains the percentiles based on that sample
  qTmatrix <- ZT_result[seq(from = 2, to = dim(ZT_result)[1], by = 2), ]
  
  # interprecentile range
  calculate_IPR <- function(DATA){
    result <- quantile(DATA, c(.025, .975))
    IPR <- result[2] - result[1]
  }

  IPR_reg[num_test, ] <- apply(qZmatrix, 2, calculate_IPR)
  IPR_Tscore[num_test, ] <- apply(qTmatrix, 2, calculate_IPR)

  num_test <- num_test + 1
}
#
