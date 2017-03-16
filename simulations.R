###################################################
####### simulation: Norming change scores #########
#######   Zhengguo Gu Tilburg University  #########
###################################################


#################################################################################
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

######################################################
sample_size = 100000  #sample size = 100, 500, 1000, 1500, 2000, 3000, 5000

propE <- .065 / 2  # proportion of explained by each predictor. .065/2, .13/2, and .26/2
beta1 <- sqrt(4*propE/(1-propE)) # For X1, dichotomous
beta2 <- sqrt(propE)             # For X2. continuous

# 1. simulate person parameters

# 1.1. covariates X1 and X2
set.seed(112)
X1 <- rbinom(sample_size, size = 1, prob = .5)

set.seed(112) 
X2 <- runif(sample_size, min = 4, max = 12)

# 1.2. expectation of theta_pretest, given X1 and X2

theta <- beta1 * X1 + beta2 * scale(X2, center = T, scale = T) + rnorm(sample_size,0,sqrt(1 - 2 * propE))    #note that X2 has to be standardized so that beta2^2 is the proportion                                                                                                                       of variance explained. Also rnorm() is added so that the var of theta    

# 1.3. extra: check whether our setup is correct
cor(X1, theta) ^ 2
cor(X2, theta) ^ 2
fit <- lm(theta ~ X1 + X2)
summary(fit)
library(lmSupport)
modelEffectSizes(fit)
