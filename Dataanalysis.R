##############################################################
#########                                        #############
######### Data Analysis: Norming change scores   #############
#########                                        #############
#########       Last update: 12/02/2019          #############
#########                                        #############
##############################################################

library(ggplot2)
library(rms)
############ 0. Some preparations ############################
# a) set the working directory to the folder containing the .RData files
filenames <- list.files(path = ".", pattern = "var_D")
# b) combine all the dataset and their conditions. 
finalmat <- list()
for(i in 1:length(filenames)){
  
  filename <- filenames[i]
  split_name <- strsplit(filename, split = "_")[[1]]
  var_D <- as.numeric(split_name[3])
  rho_preD <- as.numeric(substr(split_name[6], 1, nchar(split_name[6])-6))
  
  sample_sizeV = seq(100, 1500, by = 100)  #sample size 
  beta_pre <- rho_preD*sqrt(var_D)    #see Equation (19) in the article
  propEV <- c((.065-rho_preD^2)*var_D/2, (.13-rho_preD^2)*var_D/2, (.26-rho_preD^2)*var_D/2)  # proportion of explained by each X1 and X2 (excluding theta_pre).
                                                                                              # R^2=.065: small effect; = .13: medium effect; =.26: large effect.
                                                                                              # check Appendix C and Equation (19)
  polytomousV <- c(1, 0)  # if 1, simulate polytomous response data 
  num_itemsV <- c(10, 20, 40) 
  
  df <- expand.grid(sample_sizeV, propEV, polytomousV, num_itemsV)
  colnames(df) <- c("sample_size", "proportionExplained", "polytomous", "num_items")
  
  load(file = filename)

  finalmat[[i]] <- cbind(var_D, rho_preD, df, fresultMat)
   
}

FINALmat <- rbind(finalmat[[1]], finalmat[[2]], finalmat[[3]], finalmat[[4]], finalmat[[5]], finalmat[[6]])


summary(FINALmat)

FINALmat$proportionExplained <- FINALmat$proportionExplained * 2 /FINALmat$var_D + (FINALmat$rho_preD)^2  #transform proportionExplained back to effect sizes


############ 1. Tables: Rank correlations, using N=1000 as an example ####################
colnames(FINALmat)
FINAL_rank_1000 <- FINALmat[FINALmat$sample_size==1000, c(1:11, 21, 31:32)]

#test length 
round(colMeans(FINAL_rank_1000[FINAL_rank_1000$num_items==10,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==10,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$num_items==20,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==20,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$num_items==40,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==40,], MARGIN = 2, sd), digits = 2)

# #item scores
round(colMeans(FINAL_rank_1000[FINAL_rank_1000$polytomous==0,]), digits = 2) #dichotomous
round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==0,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$polytomous==1,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==1,], MARGIN = 2, sd), digits = 2)

#effect size of covariates
round(colMeans(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .065,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .065,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .13,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .13,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .26,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .26,], MARGIN = 2, sd), digits = 2)

#Pupulation correlation between pretest and change
round(colMeans(FINAL_rank_1000[FINAL_rank_1000$rho_preD==0,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$rho_preD==0,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$rho_preD==0.1,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$rho_preD==0.1,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$rho_preD==-0.1,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$rho_preD==-0.1,], MARGIN = 2, sd), digits = 2)

#Variance of theta-change
round(colMeans(FINAL_rank_1000[FINAL_rank_1000$var_D ==.14,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==.14,], MARGIN = 2, sd), digits = 2)

round(colMeans(FINAL_rank_1000[FINAL_rank_1000$var_D ==1.14,]), digits = 2)
round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==1.14,], MARGIN = 2, sd), digits = 2)


############### 2. IPR plots: IPR using Equation 4a and Tscores   #############################################################

colnames(FINALmat)
FINAL_IPR <- FINALmat[, c(1:6, 12:20, 22:30)]
FINAL_IPR$sample_size <- as.factor(FINAL_IPR$sample_size)
#regression-based using Equation (4a) and T Scores

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
     x ="N", y = "IPR") +
  ylim(0,1) ; p  #the ylimit of all the plots are set to (0,1), as a result, a few outliers in this plot are not shown (that's why ggplot2 generate warnings)

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x ="N", y = "IPR")+
  ylim(0,1) ; p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p


################  # ANOVAs
Data_4ANOVA <- FINALmat[, c(1:6, 12:20)]

summary(Data_4ANOVA)  #design factors need to be coded as factors (categorical)
Data_4ANOVA$var_D <- factor(Data_4ANOVA$var_D, levels = c(0.14, 1.14))
Data_4ANOVA$rho_preD <- factor(Data_4ANOVA$rho_preD, levels = c(-0.1, 0, 0.1))
#Data_4ANOVA$sample_size <- factor(Data_4ANOVA$sample_size, levels = seq(100, 1500, by = 100))
Data_4ANOVA$proportionExplained <- factor(Data_4ANOVA$proportionExplained, levels = c(0.065, 0.13, 0.26))
Data_4ANOVA$polytomous <- factor(Data_4ANOVA$polytomous, levels = c(0, 1))
Data_4ANOVA$num_items <- factor(Data_4ANOVA$num_items, levels = c(10, 20, 40))

Data_4ANOVA <- cbind(seq(1620), Data_4ANOVA)
#colnames(Data_4ANOVA)[1] <- 'celNo'  #no need for this

model1 = ols(log(`qTZ1%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result <- robcov(model1)
qqnorm(Data_4ANOVA$`qTZ1%` - result1$residuals)  #plot is quite good

model5 = ols(log(`qTZ5%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result5 <- robcov(model5)
qqnorm(Data_4ANOVA$`qTZ5%` - result5$residuals)  #normality seems to be ok

model10 = ols(log(`qTZ10%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result10 <- robcov(model10)
qqnorm(Data_4ANOVA$`qTZ10%` - result10$residuals) # plot is quite good

model25 = ols(log(`qTZ25%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result25 <- robcov(model25)
qqnorm(Data_4ANOVA$`qTZ25%` - result25$residuals) # plot is quite good

model50 = ols(log(`qTZ50%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result50 <- robcov(model50)
qqnorm(Data_4ANOVA$`qTZ50%` - result50$residuals) # plot is quite good

model75 = ols(log(`qTZ75%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result75 <- robcov(model75)
qqnorm(Data_4ANOVA$`qTZ75%` - result75$residuals) # plot seems to be ok

model90 = ols(log(`qTZ90%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result90 <- robcov(model90)
qqnorm(Data_4ANOVA$`qTZ90%` - result90$residuals) # plot seems to be ok

model95 = ols(log(`qTZ95%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result95 <- robcov(model95)
qqnorm(Data_4ANOVA$`qTZ95%` - result95$residuals) # plot seems to be ok

model99 = ols(log(`qTZ99%`) ~ var_D + rho_preD + proportionExplained + polytomous + num_items + sample_size, data = Data_4ANOVA, x=TRUE)    #use the suggestion from https://stats.stackexchange.com/questions/76904/robust-regression-inference-and-sandwich-estimators
result99 <- robcov(model99)
qqnorm(Data_4ANOVA$`qTZ99%` - result99$residuals) # plot does not looks good. But we do not deal with this issue further, because if we change the model for 99% percentile, then we have to also use the same model for the other percentiles. 
