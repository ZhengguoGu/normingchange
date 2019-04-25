library(lmSupport)

var_D <- 0.14 # .14 or 1.14   #variance of change 
rho_preD <- 0 # 0, .1, -.1    #correlation between theta_pre and theta_D

sample_sizeV = seq(100, 1500, by = 100)  #sample size 
beta_pre <- rho_preD*sqrt(var_D)    #see Equation (19) in the article
propEV <- c((.065-rho_preD^2)*var_D/2, (.13-rho_preD^2)*var_D/2, (.26-rho_preD^2)*var_D/2)  # proportion of explained by each X1 and X2 (excluding theta_pre).
                                                                                            # R^2=.065: small effect; = .13: medium effect; =.26: large effect.
# check Appendix C and Equation (19)
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
######################################################################

df <- expand.grid(sample_sizeV, propEV, polytomousV, num_itemsV)

colnames(df) <- c("sample_size", "proportionExplained", "polytomous", "num_items")

sample_size <- 10000  # in the simulation, sample size changes, but here we use a large sample size.
vec_adj_R2 <- vector()
num_test <- 1 # the 225th row of the df dataframe. 
while(num_test <= 270){
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

  # check whether our setup is correct
  #beta1^2*var(X1)+beta2^2*var(X2) + beta_pre^2 + var_Z #total variance, should be close to var_D
  #var_D
  fit <- lm(theta_D ~ X1 + X2 + theta_pre)
  vec_adj_R2[num_test] <- summary(fit)$adj.r.squared
  
  print(num_test)
  num_test <- num_test + 1
}


### Now, let's see whether the adj R2 are close to the values (.065, .13, .26) that we used in the simulation 
propEV_new <- c(0.065, 0.13, .26)  
df_new <- expand.grid(sample_sizeV, propEV_new, polytomousV, num_itemsV)
colnames(df_new) <- c("sample_size", "proportionExplained", "polytomous", "num_items")
df_new <- cbind(df_new, vec_adj_R2)
# compare the column "proportionExplained" and the vec_adj_R2. The values should be very close. 