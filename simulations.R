###################################################
####### simulation: Norming change scores #########
#######   Zhengguo Gu Tilburg University  #########
###################################################

sample_size = 100  #sample size = 100, 500, 1000, 1500, 2000, 3000, 5000
covar_effect = .065 #covariate effect = .065, .13, .26
# 1. simulate person parameters

# 1.1. covariates X1 and X2
set.seed(112)
X1 <- rbinom(sample_size, size = 1, prob = .5)

set.seed(112) 
X2 <- runif(sample_size, min = 4, max = 12)

# 1.2. expectation of theta_pretest, given X1 and X2

varX1 <- var(X1)
varX2 <- var(X2)
