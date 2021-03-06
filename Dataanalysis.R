##############################################################
#########                                        #############
######### Data Analysis: Norming change scores   #############
#########                                        #############
#########       Last update: 03/2019             #############
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


############ 1a. Table: Rank correlations, using N=1000 as an example ####################
colnames(FINALmat)
FINAL_rank_1000 <- FINALmat[FINALmat$sample_size==1000, c(1:11, 21, 31:32)]

#test length 
row1 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==10,], MARGIN = 2, median), digits = 2)
row2 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==20, ], MARGIN = 2, median), digits = 2)
row3 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==40,], MARGIN = 2, median), digits = 2)

# #item scores
row4 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==0,], MARGIN = 2, median), digits = 2)
row5 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==1,], MARGIN = 2, median), digits = 2)

#effect size of covariates
row6 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .065,], MARGIN = 2, median), digits = 2)
row7 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .13 ,], MARGIN = 2, median), digits = 2)
row8 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .26 ,], MARGIN = 2, median), digits = 2)

#Population correlation between theta_pre and theta_D
row9 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$rho_preD ==0 ,], MARGIN = 2, median), digits = 2)
row10 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$rho_preD == -0.1 ,], MARGIN = 2, median), digits = 2)
row11 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$rho_preD == 0.1 ,], MARGIN = 2, median), digits = 2)

#Variance of theta-change
row12 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==.14 ,], MARGIN = 2, median), digits = 2)
row13 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==1.14 ,], MARGIN = 2, median), digits = 2)

table1 <- rbind(row1, row2, row3, row4, row5,
                row6, row7, row8, row9, row10,
                row11, row12, row13)

write.csv(table1, file = "table1_afterdisucssingwithWilco.csv")
#comments: in the csv: rank_cor_ZvecT = r_E,T, rank_cor_ZvecregNoXpre = r_E,D, rank_cor_TandregNoXpre = D,T

############### 1b. plot: rank correlation vs. sample size
# Rank correlation r_zt
FINAL_rank_forplot <- FINALmat[, c(1:6, 21)]
FINAL_rank_forplot$sample_size <- as.factor(FINAL_rank_forplot$sample_size)
p <- ggplot(FINAL_rank_forplot, aes(x=FINAL_rank_forplot$sample_size, y=FINAL_rank_forplot$rank_cor_ZvecT)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), axis.line = element_line(colour = "black")) +
  labs(title="(b) T Scores for Change method", x ="N", y = "Rank Correlation") +
  ylim(0,1) ; 
p 
# Rank correlation r_zd 
FINAL_rank_forplot <- FINALmat[, c(1:6, 31)]
FINAL_rank_forplot$sample_size <- as.factor(FINAL_rank_forplot$sample_size)
p <- ggplot(FINAL_rank_forplot, aes(x=FINAL_rank_forplot$sample_size, y=FINAL_rank_forplot$rank_cor_ZvecregNoXpre)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), axis.line = element_line(colour = "black")) +
  labs(title="(a) Regression-based change approach", x ="N", y = "Rank Correlation") +
  ylim(0,1) ; 
p 

############### 2a. IPR plots: IPR generated by T Score for change method #############################################################

colnames(FINALmat)
FINAL_IPR <- FINALmat[, c(1:6, 12:20)]  
FINAL_IPR$sample_size <- as.factor(FINAL_IPR$sample_size)

library(grid)
library(gridExtra)
#regression-based using Equation (4a) and T Scores

p1 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
     x ="N", y = "IPR") +
  ylim(0,1.2) ; p1 #Note! the y limit is truncated, to be consistent with the plots below 

p2 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2); p2

p3 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2) ; p3

p4 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2); p4

p5 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2); p5

p6 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2); p6

p7 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2); p7

p8 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2); p8

p9 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qTZ99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x ="N", y = "IPR")+
  ylim(0,1.2); p9

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
)

################## 2b. more IPR box plots. because ANOVAs are not suitable (assumptions are violated.)  ###########################################
### 2b.1 variance of theta change #############################################################################
colnames(FINALmat)
FINAL_IPR <- FINALmat[, c(1:6, 12:20)]  
FINAL_IPR$var_D <- as.factor(FINAL_IPR$var_D)

p1 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p1 #note: the reason of warming messages: the ylim is truncated. 

p2 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p2

p3 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p3

p4 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p4

p5 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p5

p6 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p6

p7 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p7

p8 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p8

p9 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$var_D, y=FINAL_IPR$`qTZ99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x =expression(paste("Variance of ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p9

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
)

### 2b.2 test length  #############################################################################
FINAL_IPR$num_items <- as.factor(FINAL_IPR$num_items)

p1 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p1

p2 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p2

p3 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p3

p4 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p4

p5 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p5

p6 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p6

p7 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p7

p8 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p8

p9 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$num_items, y=FINAL_IPR$`qTZ99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x ="Test Length", y = "IPR") +
  ylim(0,1.2) ; p9

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
)

### 2b.3 item scores (dich/poly)  #############################################################################
FINAL_IPR$polytomous <- as.factor(FINAL_IPR$polytomous)
New_lables <- c(2, 5)
p1 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p1

p2 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p2

p3 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p3

p4 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p4

p5 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p5

p6 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p6

p7 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p7

p8 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p8

p9 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$polytomous, y=FINAL_IPR$`qTZ99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x ="Number of Item Scores", y = "IPR") +
  ylim(0,1.2) +
  scale_x_discrete(labels= New_lables); p9

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
)
### 2b.4 Effect size of covariates ###########################
FINAL_IPR$proportionExplained <- as.factor(FINAL_IPR$proportionExplained)

p1 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p1

p2 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p2

p3 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p3

p4 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p4

p5 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p5

p6 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p6

p7 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p7

p8 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p8

p9 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$proportionExplained, y=FINAL_IPR$`qTZ99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x ="Effect Size of Covariates", y = "IPR") +
  ylim(0,1.2) ; p9

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
)

### 2b.5 Correlation between theta_pre and theta_change ###########################
FINAL_IPR$rho_preD <- as.factor(FINAL_IPR$rho_preD)

p1 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p1

p2 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p2

p3 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p3

p4 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p4

p5 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p5

p6 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p6

p7 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p7

p8 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p8

p9 <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$rho_preD, y=FINAL_IPR$`qTZ99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x =expression(paste("Correlation Between ", theta[pre], " and ", theta[D])), y = "IPR") +
  ylim(0,1.2) ; p9

grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  nrow = 3
)

######################### END  ##################################

############### extra! Not in the paper!  #######################################################################
# examine IPR generated by regression-based change approach (i.e., X_pre is excluded from the regression)
colnames(FINALmat)
FINAL_IPR <- FINALmat[, c(1:6, 22:30)]
FINAL_IPR$sample_size <- as.factor(FINAL_IPR$sample_size)

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre1%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="1st percentile",
       x ="N", y = "IPR") +
  ylim(0,1) ; p  #the ylimit of all the plots are set to (0,1), as a result, a few outliers in this plot are not shown (that's why ggplot2 generate warnings)

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre5%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="5th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre10%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="10th percentile",
       x ="N", y = "IPR")+
  ylim(0,1) ; p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre25%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="25th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre50%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="50th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre75%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="75th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre90%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="90th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre95%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="95th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p

p <- ggplot(FINAL_IPR, aes(x=FINAL_IPR$sample_size, y=FINAL_IPR$`qZ_noXpre99%`)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  theme(text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) +
  labs(title="99th percentile",
       x ="N", y = "IPR")+
  ylim(0,1); p


############ Note used. Table: Rank correlations, using N=1000 as an example (rho_preD=0 and !=0) ####################
colnames(FINALmat)
FINAL_rank_1000 <- FINALmat[FINALmat$sample_size==1000, c(1:11, 21, 31:32)]

#when theta_pre is not included in the population model: i.e., rho_preD = 0####
index_rhopre <- FINAL_rank_1000$rho_preD == 0
#test length 
row1 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==10 & index_rhopre,], MARGIN = 2, median), digits = 2)
row2 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==20 & index_rhopre, ], MARGIN = 2, median), digits = 2)
row3 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==40 & index_rhopre,], MARGIN = 2, median), digits = 2)

# #item scores
row4 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==0 & index_rhopre,], MARGIN = 2, median), digits = 2)
row5 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==1 & index_rhopre,], MARGIN = 2, median), digits = 2)

#effect size of covariates
row6 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .065 & index_rhopre,], MARGIN = 2, median), digits = 2)
row7 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .13 & index_rhopre,], MARGIN = 2, median), digits = 2)
row8 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .26 & index_rhopre,], MARGIN = 2, median), digits = 2)

#Variance of theta-change
row9 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==.14 & index_rhopre,], MARGIN = 2, median), digits = 2)
row10 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==1.14 & index_rhopre,], MARGIN = 2, median), digits = 2)

table1 <- rbind(row1, row2, row3, row4, row5,
                row6, row7, row8, row9, row10)

#when theta_pre is not included in the population model: i.e., rho_preD != 0####
index_rhopre <- FINAL_rank_1000$rho_preD != 0
#test length 
row1 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==10 & index_rhopre,], MARGIN = 2, median), digits = 2)
row2 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==20 & index_rhopre, ], MARGIN = 2, median), digits = 2)
row3 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$num_items==40 & index_rhopre,], MARGIN = 2, median), digits = 2)

# #item scores
row4 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==0 & index_rhopre,], MARGIN = 2, median), digits = 2)
row5 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$polytomous==1 & index_rhopre,], MARGIN = 2, median), digits = 2)

#effect size of covariates
row6 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .065 & index_rhopre,], MARGIN = 2, median), digits = 2)
row7 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .13 & index_rhopre,], MARGIN = 2, median), digits = 2)
row8 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$proportionExplained == .26 & index_rhopre,], MARGIN = 2, median), digits = 2)

#Variance of theta-change
row9 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==.14 & index_rhopre,], MARGIN = 2, median), digits = 2)
row10 <- round(apply(FINAL_rank_1000[FINAL_rank_1000$var_D ==1.14 & index_rhopre,], MARGIN = 2, median), digits = 2)

table2 <- rbind(row1, row2, row3, row4, row5,
                row6, row7, row8, row9, row10)
