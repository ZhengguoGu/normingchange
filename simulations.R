###################################################
####### simulation: Norming change scores #########
#######   Zhengguo Gu Tilburg University  #########
###################################################


sample_size = 10000000  #sample size = 100, 500, 1000, 1500, 2000, 3000, 5000
propE <- .26 / 2  # proportion of explained by each predictor. 0, .065/2, .13/2, and .26/2
polytomous <- TRUE  # if true, simulate polytomous response data 
num_items <- 10 # 10, 20, and 40

##################### Function for generating responses ##############################
GRM_sim <- function(ability, itempar){
  
  # descrption:
  #
  # ability = ability parameter
  # itempar = item parameter
    n_sub <- length(ability)
    response <- matrix(NA, n_sub, nrow(itempar))
    true_response <- matrix(NA, n_sub, nrow(itempar))
    
    # if true, then unidimensional 
    for(i in 1:n_sub){
      
      numeritor <- exp(sweep((ability[i]-itempar[, -1]), 1, itempar[, 1], "*"))
      P_star <- numeritor/(1+numeritor) # this is the "true response"
      
      response[i, ] <- rowSums(P_star >= runif(nrow(itempar), min=0, max=1))
      true_response[i, ] <- rowSums(P_star)
      
    }
 
  return(list(response, true_response))
}

####################### theta_D ###########################################

beta1 <- sqrt(propE * .14 / .25)  # For X1, dichotomous
beta2 <- sqrt(propE * .14 * 12 / (12 - 4)^2)             # For X2. continuous, uniform

# 1. simulate person parameters

# 1.1. covariates X1 and X2
set.seed(112)
X1 <- rbinom(sample_size, size = 1, prob = .5)
X2 <- runif(sample_size, min = 4, max = 12)

# 1.2. expectation of theta_pretest, given X1 and X2

theta_D <- beta1 * X1 + beta2 * X2 + rnorm(sample_size, .75, sqrt(.14 - 2 * propE * .14))
var(theta_D)
# 1.3. extra: check whether our setup is correct
cor(X1, theta_D) ^ 2
cor(X2, theta_D) ^ 2

#fit <- lm(theta_D ~ X1 + X2)
#summary(fit)
#library(lmSupport)
#modelEffectSizes(fit)

###################### Simulate response data #######################
theta_pre <- rnorm(sample_size, 0, 1)
theta_post <- theta_pre + theta_D

if(polytomous == TRUE){
  
  itempar <- matrix(NA,num_items,5)
  itempar[,1] <- runif(num_items,1.5,2.5)   # discrimination
  avg_beta <- runif(num_items, 0, 1.25)
  itempar[,2] <- avg_beta - 1
  itempar[,3] <- avg_beta - .5
  itempar[,4] <- avg_beta + .5
  itempar[,5] <- avg_beta + 1
  
} else{
  
  itempar <- matrix(NA,num_items,2)
  itempar[,1] <- runif(num_items,1.5,2.5)   # discrimination
  itempar[,2] <- runif(num_items, 0, 1.25)
}



