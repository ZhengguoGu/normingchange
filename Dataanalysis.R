##############################################################
#########                                        #############
######### Data Analysis: Norming change scores   #############
#########                                        #############
#########      Zhengguo Gu, Tilburg University   #############
#########      Last update: 05/05/2017           #############
##############################################################


####################### PART I: Calculating IPR ####################################

####################################################################################
###### NOTE: PART I was used to generate the dataset.                         ######
###### 20170505IPR.RData.RData                                                ######
####################################################################################

library(foreach)
library(doSNOW)
library(doRNG)
 
IPR_reg <- matrix(NA, 360, 9) 
IPR_Tscore <- matrix(NA, 360, 9)
REL_changescore <- matrix(NA, 360, 2)  #first column: averaged reliability, second column: sd
num_test <- 1
while(num_test <= 360){
  
  # 1. reorganize the data.
  # To save space, for each test condition (in total 360), pretest and posttest
  # scores of 1000 datasets are 'cbind'ed. For example, the results of the first 
  # test condition:
  
  filename <- paste("D:/ZhengguoProj2Blad analysis/results_", num_test, ".RData", sep = "")
  load(filename)

  #dim(sim_result)  #100 x 2000 matrix. 100: 100 persons. 2000: 1000 datasets, each contain a pretest and posttest
                 #score, and thus 1000x2=2000

  # 2. standardize scores. 
  mydata <- scale(sim_result, center = TRUE, scale = TRUE)  #note that sim_result contains estimated reliability. But recale does not change the results of reliability (because the 
                                                            # entire column has the same reliability) But somehow the reliabilities in mydata matrix have been replaced by NaN.

  # identical((sim_result[, 22] - mean(sim_result[, 22]))/sd(sim_result[, 22]), mydata[,22])

  # 3. reorganize the data 
  datalist <- list()
  rellist <- array()
  j <- 1
  for(i in seq(from = 1, to = 3000, by = 3)){
    datalist[[j]] <- mydata[, c(i, i+1)]
    rellist[j] <- sim_result[1, i+2]   #this is reliability, we pick the first row of the column (or whichever row we like), because the entire column has the same value.
    j <- j + 1
  }
  REL_changescore[num_test, 1] <- mean(rellist)
  REL_changescore[num_test, 2] <- sd(rellist)
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
colnames(IPR_reg) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")
colnames(IPR_Tscore) <- c("1%", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "99%")

save(IPR_reg, IPR_Tscore, REL_changescore, file = "20170505IPR.RData")


########################### PART II: ANOVA ###########################################################

load(file ="D:/ZhengguoProj2Blad analysis/20170505IPR.RData")

# 1. get the design factors.
# Note: The following code is from simulationsNEW.R, where we know which design factors are in which cell.
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

# 2. Combine results with design factors.

IPR_regData <- cbind(IPR_reg, df, "regression-based", seq(1:360), REL_changescore)
colnames(IPR_regData)[c(14, 15, 16, 17)] <- c("NormingMethod", "simulationcell", "ChangeRel", "ChangeRel1SD")

IPR_TData <- cbind(IPR_Tscore, df, "TScore", seq(1:360), REL_changescore)
colnames(IPR_TData)[c(14, 15, 16, 17)] <- c("NormingMethod", "simulationcell", "ChangeRel", "ChangeRel1SD")

IPR_Data <- rbind(IPR_regData, IPR_TData)

# 3. mixed anova 
library(DescTools)  #to calculate eta squared

summary(IPR_Data)  #design factors need to be coded as factors (categorical)
IPR_Data$num_items <- factor(IPR_Data$num_items, levels = c(10, 20, 40))
IPR_Data$polytomous <- factor(IPR_Data$polytomous, levels = c(0, 1))
IPR_Data$proportionExplained <- factor(IPR_Data$proportionExplained, levels = c(0, .065/2, .13/2, .26/2))
IPR_Data$simulationcell <- factor(IPR_Data$simulationcell, levels = seq(1:360))
summary(IPR_Data) #check again

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
write.table(etamatrix, file = 'D:/ZhengguoProj2Blad analysis/etamatrix.txt', sep = ',')


# 4. relationship between sample size and IPRs
IPR_Data # make use of this datamatrix

aggreIPR <- aggregate(IPR_Data[, 1:9], by = list(IPR_Data$'sample_size', IPR_Data$'NormingMethod'), FUN = mean) 

#plot

maintext = paste(c("1th", "5th", "10th", "25th", "50th", "75th", "90th", "95th", "99th"), "percentile")

for(i in 1:9){
  
  filename <- paste("D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170419/IPR_", i, ".png", sep = "")
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
