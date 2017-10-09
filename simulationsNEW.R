###################################################
####### simulation: Norming change scores #########
#######   Zhengguo Gu Tilburg University  #########
#######       Last update: 01/05/2017     #########
###################################################

library(foreach)
library(psychometric)
library(doSNOW)
library(doRNG)


tmp=proc.time()
sample_sizeV = c(100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 9000, 10000)  #sample size 
propEV <- c(0, .065/2, .13/2, .26/2)  # proportion of explained by each predictor. 
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
#while(num_test <= 2){ 
  
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
  
  beta1 <- sqrt(propE * .14 / .25)  # For X1, dichotomous
  beta2 <- sqrt(propE * .14 * 12 / (12 - 4)^2)             # For X2. continuous, uniform
  
  # 1. simulate person parameters
  
  # 1.1. covariates X1 and X2
  X1 <- rbinom(sample_size, size = 1, prob = .5)
  X2 <- runif(sample_size, min = 4, max = 12)
  
  # 1.2. expectation of theta_pretest, given X1 and X2
  
  theta_D <- beta1 * X1 + beta2 * X2 + rnorm(sample_size, .75, sqrt(.14 - 2 * propE * .14))
  
  #var(theta_D)
  # 1.3. extra: check whether our setup is correct
  #cor(X1, theta_D) ^ 2
  #cor(X2, theta_D) ^ 2
  
  #fit <- lm(theta_D ~ X1 + X2)
  #summary(fit)
  #library(lmSupport)
  #modelEffectSizes(fit)
  
  cl <- makeCluster(12)
  registerDoSNOW(cl)

  #note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
  sim_result <- foreach(i = 1:1000, .combine='cbind') %dorng% {
    
    ###################### Simulate response data #######################
  
    theta_pre <- rnorm(sample_size, 0, 1)
    theta_post <- theta_pre + theta_D
    
    
    X_pre <- GRM_sim(theta_pre, itempar)[[1]]
    X_post <- GRM_sim(theta_post, itempar)[[1]]
    
    sum_pre <- rowSums(X_pre)
    sum_post <- rowSums(X_post)
    
    Difference_item <- X_post - X_pre
    change_rel <- psychometric::alpha(Difference_item)
    list_sum <- cbind(sum_pre, sum_post, change_rel)
    
    return(list_sum)
  }

  stopCluster(cl)
  
  
  filename <- paste("results_", num_test, ".RData", sep = "")
  save(X1, X2, sim_result, file = filename)
  itemparamter <- paste("item_", num_test, ".RData", sep = "")
  save(item_par, file = itemparamter)
  num_test <- num_test + 1 
}


proc.time()-tmp
#save(df, results_responseD, )
