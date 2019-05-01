###################################################
####### simulation: Norming change scores #########
#######                                   #########
#######       Last update: 25/04/2019     #########
###################################################

library(foreach)
library(psychometric)
library(doSNOW)
library(doRNG)
library(Kendall)

#########NEW Parameters that need to be adjusted manually (for 1st revision at Assessment) ################
var_D <- 0.14 # .14 or 1.14   #variance of change 
rho_preD <- -0.1 # 0, .1, -.1    #correlation between theta_pre and theta_D
#######################################################################################################

tmp=proc.time()


set.seed(10) # set seed
sample_sizeV = seq(100, 1500, by = 100)  #sample size 
beta_pre <- rho_preD*sqrt(var_D)    #see Equation (15) in the article
propEV <- c((.065-rho_preD^2)*var_D/2, (.13-rho_preD^2)*var_D/2, (.26-rho_preD^2)*var_D/2)  # proportion of explained by each X1 and X2 (excluding theta_pre).
                                                                                            # R^2=.065: small effect; = .13: medium effect; =.26: large effect.
                                                                                            # check Appendix C
polytomousV <- c(1, 0)  # if 1, simulate polytomous response data 
num_itemsV <- c(10, 20, 40) 




##################### Function for generating responses ##############################
GRM_sim <- function(ability, itempar){
  
  # descrption:
  #
  # ability = ability parameter
  # itempar = item parameter
    n_sub <- length(ability)
    response <- matrix(NA, n_sub, nrow(itempar))
    true_response <- matrix(NA, n_sub, nrow(itempar))
    
    
    for(i in 1:n_sub){
      if(dim(itempar)[2] > 2){
        
        #polytomous
        numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))  #note: the first column of the itempar matrix contains the discrimination parameters; the rest columsn thresholds.
        P_star <- numeritor/(1+numeritor) # this is the "true response"
        
        response[i, ] <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))  # this is the simulated item response
        true_response[i, ] <- rowSums(P_star)  #sum of the probabilities, this is regarded as the true score
      } else{
        
        #dichotomous 
        numeritor <- exp((ability[i] - matrix(itempar[, -1], ncol = 1)) * matrix(itempar[, 1], ncol = 1))
        P_star <- numeritor/(1+numeritor)
        
        response[i, ] <- rowSums(P_star >= matrix(runif(nrow(itempar), min=0, max=1)))
        true_response[i, ] <- rowSums(P_star)
      }
      
      
    }
 
  return(list(response, true_response))
}

############## interprecentile range #########################
calculate_IPR <- function(vec_data){
  result <- quantile(vec_data, c(.025, .975))
  IPR <- result[2] - result[1]
}
##############################################################


########################################################
#######                                   ##############
#######     Now starts the simulation     ##############
#######                                   ##############
########################################################


df <- expand.grid(sample_sizeV, propEV, polytomousV, num_itemsV)

colnames(df) <- c("sample_size", "proportionExplained", "polytomous", "num_items")

#results_responseD <- list() dont need it anymore
#item_par <- list()  # do not save this anymore
fresultMat <- matrix(NA, nrow(df), 27)
num_test <- 1
while(num_test <= nrow(df)){

  print(num_test)
  
  sample_size <- df[num_test, 1]
  propE <- df[num_test, 2]
  polytomous <- df[num_test, 3]
  num_items <- df[num_test, 4]
  
  #### first draw items from the item population, afterwards the items are fixed.
  
  if(polytomous == 1){
    
    itempar <- matrix(NA,num_items,5)
    itempar[,1] <- runif(num_items,1.5,2.5)   # discrimination
    avg_beta <- runif(num_items, 0, 1.25)
    itempar[,2] <- avg_beta - .75
    itempar[,3] <- avg_beta - .25
    itempar[,4] <- avg_beta + .25
    itempar[,5] <- avg_beta + .75
    
  } else{
    
    itempar <- matrix(NA,num_items,2)
    itempar[,1] <- runif(num_items,1.5,2.5)   # discrimination
    itempar[,2] <- runif(num_items, 0, 1.25)
  }
  
  #item_par[[num_test]] <- itempar #we dont save this anymore
  
  beta1 <- sqrt(propE / .25)  # For X1, dichotomous; see Appendix C
  beta2 <- sqrt(propE * 12 / (12 - 4)^2)    # For X2. continuous, uniform; see Appendix C
  
  # 1. simulate person parameters
  
  # 1.1. covariates X1 and X2
  X1 <- rbinom(sample_size, size = 1, prob = .5)
  X2 <- runif(sample_size, min = 4, max = 12)
  
  # 1.2. expectation of theta_pretest, given X1 and X2, and theta_pre
  
  theta_pre <- rnorm(sample_size, 0, 1)
  
  if(propE== (.065-rho_preD^2)*var_D/2){   #propEV <- c((.065-rho_preD^2)*var_D/2, (.13-rho_preD^2)*var_D/2, (.26-rho_preD^2)*var_D/2)
    var_Z <- var_D * (1 - 0.065)
  }else if(propE == (.13-rho_preD^2)*var_D/2){
    var_Z <- var_D * (1 - 0.13)
  }else if (propE == (.26-rho_preD^2)*var_D/2){
    var_Z <- var_D * (1 - 0.26)
  }
  
  vector_of_Z <- rnorm(sample_size, .75, sqrt(var_Z))  #need to record this
  theta_D <- beta1 * X1 + beta2 * X2 + beta_pre * theta_pre + vector_of_Z # notice that beta_0 = .75 is included when we generate Z
  theta_post <- theta_pre + theta_D

  
  cl <- makeCluster(4)
  registerDoSNOW(cl)

  #note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
  sim_result <- foreach(i = 1:1000, .combine='rbind') %dorng% {
    
    ###################### Simulate response data #######################
    
    X_pre <- GRM_sim(theta_pre, itempar)
    #result <- ltm::grm(X_pre)
    #iParameters_hat <- matrix(unlist(result$coefficients), ncol = 5, byrow = TRUE)
    #plot(iParameters_hat[, 5], itempar[, 1])
    #plot(iParameters_hat[, 1], itempar[, 2])
    #plot(iParameters_hat[, 2], itempar[, 3])
    #plot(iParameters_hat[, 3], itempar[, 4])
    X_post <- GRM_sim(theta_post, itempar)
    t_pre <- X_pre[[2]]
    t_post <- X_post[[2]]
    sum_pre <- rowSums(X_pre[[1]])
    sum_post <- rowSums(X_post[[1]])
    Difference_sumscores <- sum_post - sum_pre
    
    var_1 <- var(sum_pre)
    var_2 <- var(sum_post)
    #cor_12 <- cor(sum_pre, sum_post)
    #change_rel <- (psychometric::alpha(X_pre)*var_1 + psychometric::alpha(X_post)*var_2 - 2*cor_12*sqrt(var_1*var_2))/(var_1 + var_2 - 2*cor_12*sqrt(var_1*var_2))
    #cor_preD <- cor(sum_pre, Difference_sumscores)  #correlation between observed pretest and observed change
      
    #Rel_ect <- c(change_rel, var_1, var_2, cor_12, cor_preD )
    
    
    ####### regression-based change approach, comparable to T-Scores
    #fit <- lm(Difference_sumscores ~ sum_pre + X1 + X2)
    #coef_regbased <- summary(fit)$coefficients[c(2,4), 1]    #no need to record coefficients
    #Escore <- Difference_sumscores - predict(fit) #residual
    #SD_e <- sqrt(sum(Escore^2)/(length(Escore) - 4))  #K+2 = 4
    #Zscore <- Escore/SD_e
    #qZ <- quantile(Zscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    
    ####### T-scores
    #fit2 <- lm(sum_post ~ sum_pre + X1 + X2)
    #coef_Tscore <- summary(fit2)$coefficients[c(2, 3, 4), 1]
    #Escore2 <-  sum_post - predict(fit2)
    #SD_e2 <- sqrt(sum(Escore2^2)/(length(Escore2) - 4))
    #Tscore <- Escore2/SD_e2
    #qTZ <- quantile(Tscore, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))  #note, i call it qTZ because in this case, the quantiles here are equal to the regression-based approach including Xpre
    
    #TZ_same <- mean(round(Zscore, digits = 3) == round(Tscore, digits = 3))   #this is trivial: I expect it to be true, both methods include Xpre. Therefore TZ_same == 1
    ### rank correlations, T-score and vector of Z
    #rank_cor_ZT <- Kendall::Kendall(vector_of_Z, Tscore)$tau[1]
    
    ####### extra, if regression-based change approach does not include X_pre
    #fit3 <- lm(Difference_sumscores ~  X1 + X2)
    #Escore3 <- Difference_sumscores - predict(fit3) #residual
    #SD_e <- sqrt(sum(Escore3^2)/(length(Escore3) - 3))  #X_pre is removed K+2 = 3 now
    #Zscore_noXpre <- Escore3/SD_e
    #qZ_noXpre <- quantile(Zscore_noXpre, c(.01, .05, .1, .25, .50, .75, .90, .95, .99))
    #rank_cor_ZregNoXpre <- Kendall::Kendall(vector_of_Z, Zscore_noXpre)$tau[1]
    #rank_cor_TandregNoXpre <- Kendall::Kendall(Tscore, Zscore_noXpre)$tau[1]
  
    final_sim <- c(var_1, var_2)
    
    return(final_sim)
  }

  stopCluster(cl)
  
  #colnames(sim_result) <- c("change_rel", "var_1", "var_2", "cor_12", "cor_preD", 
                           # "qTZ1%", "qTZ5%", "qTZ10%", "qTZ25%", "qTZ50%", "qTZ75%", "qTZ90%", "qTZ95%", "qTZ99%",
                           # "rank_cor_ZvecT",
                           # "qZ_noXpre1%", "qZ_noXpre5%", "qZ_noXpre10%", "qZ_noXpre25%", "qZ_noXpre50%", "qZ_noXpre75%", "qZ_noXpre90%", "qZ_noXpre95%", "qZ_noXpre99%",
                           # "rank_cor_ZvecregNoXpre", "rank_cor_TandregNoXpre", "TZ_same")
  
  column_meansANDipr <- colMeans(sim_result)
  IPR_qTZ <- apply(sim_result[,6:14], 2, calculate_IPR)
  column_meansANDipr[6:14]<- IPR_qTZ
  IPR_qZ_noXpre <- apply(sim_result[,16:24], 2, calculate_IPR)
  column_meansANDipr[16:24]<- IPR_qZ_noXpre
  
  fresultMat[num_test,] <- column_meansANDipr
  
  num_test <- num_test + 1 
}

colnames(fresultMat) <- c("change_rel", "var_1", "var_2", "cor_12", "cor_preD", 
                          "qTZ1%", "qTZ5%", "qTZ10%", "qTZ25%", "qTZ50%", "qTZ75%", "qTZ90%", "qTZ95%", "qTZ99%",
                          "rank_cor_ZvecT",
                          "qZ_noXpre1%", "qZ_noXpre5%", "qZ_noXpre10%", "qZ_noXpre25%", "qZ_noXpre50%", "qZ_noXpre75%", "qZ_noXpre90%", "qZ_noXpre95%", "qZ_noXpre99%",
                          "rank_cor_ZvecregNoXpre", "rank_cor_TandregNoXpre", "TZ_same")
filename <- paste("var_D_", var_D,"_rho_preD_",rho_preD, ".RData", sep = "")
save(fresultMat, file = filename)

proc.time()-tmp


x = GRM_sim(rep(0, 100000), itempar)[[1]]
apply(x, 2, mean)
