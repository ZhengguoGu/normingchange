###################################################
####### simulation: Norming change scores #########
#######   Zhengguo Gu Tilburg University  #########
#######       Last update: 30/05/2018     #########
###################################################

library(foreach)
library(psychometric)
library(doSNOW)
library(doRNG)
library(Kendall)

#########NEW Parameters that need to be adjusted manually (1st revision at Assessment) ################
var_D <- 1.14 # .14 or 1.14
rho_preD <- .1 # 0, .1, -.1
#######################################################################################################

tmp=proc.time()
sample_sizeV = c(100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 9000, 10000)  #sample size 
beta_pre <- rho_preD*sqrt(var_D) 
propEV <- c((.065-rho_preD^2)*var_D/2, (.13-rho_preD^2)*var_D/2, (.26-rho_preD^2)*var_D/2)  # proportion of explained by each predictor. 
polytomousV <- c(TRUE, FALSE)  # if true, simulate polytomous response data 
num_itemsV <- c(10, 20, 40) # 10, 20, and 40




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
        
        numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))
        P_star <- numeritor/(1+numeritor) # this is the "true response"
        
        response[i, ] <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
        true_response[i, ] <- rowSums(P_star)
      } else{
        
        numeritor <- exp((ability[i] - matrix(itempar[, -1], ncol = 1)) * matrix(itempar[, 1], ncol = 1))
        P_star <- numeritor/(1+numeritor)
        
        response[i, ] <- rowSums(P_star >= matrix(runif(nrow(itempar), min=0, max=1)))
        true_response[i, ] <- rowSums(P_star)
      }
      
      
    }
 
  return(list(response, true_response))
}

########################################################
#######                                   ##############
#######     Now starts the simulation     ##############
#######                                   ##############
########################################################

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

#results_responseD <- list() dont need it anymore
item_par <- list()

num_test <- 1
while(num_test <= nrow(df)){
 
  
  set.seed(112)  # set seed
  
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
  
  item_par[[num_test]] <- itempar
  
  ####################### theta_D ###########################################
  
  beta1 <- sqrt(propE / .25)  # For X1, dichotomous
  beta2 <- sqrt(propE * 12 / (12 - 4)^2)    # For X2. continuous, uniform
  
  # 1. simulate person parameters
  
  # 1.1. covariates X1 and X2
  X1 <- rbinom(sample_size, size = 1, prob = .5)
  X2 <- runif(sample_size, min = 4, max = 12)
  
  # 1.2. expectation of theta_pretest, given X1 and X2
  
  theta_pre <- rnorm(sample_size, 0, 1)
  
  if(propE== (.065-rho_preD^2)*var_D/2){
    var_Z <- var_D * (1 - 0.065)
  }else if(propE == (.13-rho_preD^2)*var_D/2){
    var_Z <- var_D * (1 - 0.13)
  }else if (propE == (.26-rho_preD^2)*var_D/2){
    var_Z <- var_D * (1 - 0.26)
  }
  theta_D <- beta1 * X1 + beta2 * X2 + beta_pre * theta_pre + rnorm(sample_size, .75, sqrt(var_Z))
  theta_post <- theta_pre + theta_D
  
  
  #var(theta_D)
  # 1.3. extra: check whether our setup is correct
  #beta1^2*var(X1)+beta2^2*var(X2) + beta_pre^2 + var_Z #total variance, should be close to var_D
  #var_D
  #fit <- lm(theta_D ~ X1 + X2 + theta_pre)
  #summary(fit)
  #library(lmSupport)
  #modelEffectSizes(fit)
  
  ycl <- makeCluster(12)
  registerDoSNOW(cl)

  #note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
  sim_result <- foreach(i = 1:1000, .combine='cbind') %dorng% {
    
    ###################### Simulate response data #######################
    
    X_pre <- GRM_sim(theta_pre, itempar)[[1]]
    #result <- ltm::grm(X_pre)
    #iParameters_hat <- matrix(unlist(result$coefficients), ncol = 5, byrow = TRUE)
    #plot(iParameters_hat[, 5], itempar[, 1])
    #plot(iParameters_hat[, 1], itempar[, 2])
    #plot(iParameters_hat[, 2], itempar[, 3])
    #plot(iParameters_hat[, 3], itempar[, 4])
    X_post <- GRM_sim(theta_post, itempar)[[1]]
    
    sum_pre <- rowSums(X_pre)
    sum_post <- rowSums(X_post)
    
    Difference_item <- X_post - X_pre
    Differecen_sum <- sum_post - sum_pre
    
    change_rel <- psychometric::alpha(Difference_item)
    var_1 <- var(sum_pre)
    var_2 <- var(sum_post)
  
    list_sum <- cbind(sum_pre, sum_post, change_rel, var_1, var_2)
    
    return(list_sum)
  }

  stopCluster(cl)
  
  
  filename <- paste("results_", num_test, ".RData", sep = "")
  save(theta_D, sim_result, file = filename)
  beta_paramter <- paste("beta_", num_test, ".RData", sep = "")
  save(beta_pre, beta1, beta2, file = beta_paramter)
  num_test <- num_test + 1 
}


proc.time()-tmp
#save(df, results_responseD, )
