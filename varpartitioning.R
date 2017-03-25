g <- c(rep(0, 500000), rep(1, 500000))
x <- rnorm(1000000, 0, 1)

th <- .5 * x + rnorm(1000000, 0, sqrt(.75))

cor(x, th) ^ 2

results <- lm(th ~ x)
summary(results)

th <- sqrt(4 * .25 /.75) * g + rnorm(1000000, 0, sqrt(1))
cor(g, th)^2
results <- lm(th ~ g)
summary(results)

library(lmSupport)
modelEffectSizes(results)

b1 = sqrt(0.2)
b2 = sqrt(4 * .2/(1 - .2))      

n = 10000000
x1 = runif(n, min = 4, max = 12)
g = c(rep(0,n/2),rep(1,n/2))

x1 = rnorm(n)

th = b1*scale(x1, center = T, scale = T) + b2*g + rnorm(n,0,sqrt(.6))
tht = b1*x1 + b2*g
cor(th,tht)^2
cor(th,g)^2
cor(th,x1)^2

results <- lm(th ~ g + x1)
summary(results)


########################  situation 1

x <- rnorm(1000000, 0, 1)
x2 <- rnorm(1000000, 0, 1)
th <- .5 * x  + .5 * x2 + rnorm(1000000, 0, sqrt(.5))

cor(x, th) ^ 2
cor(x2, th) ^ 2

results <- lm(th ~ x + x2)
summary(results)
modelEffectSizes(results)

########################  situation 2

g <- c(rep(0, 500000), rep(1, 500000))

th <- sqrt(4 * .25 /.75) * g + rnorm(1000000, 0, sqrt(1)) #note rnorm(, 0, sqrt(1) not sqrt(.75)) !!
cor(g, th)^2
g <- factor(g, levels = c(0,1))
results <- lm(th ~ g)
summary(results)
modelEffectSizes(results)

#######################  situation 3

x <- rnorm(1000000, 0, 1)
g <- c(rep(0, 500000), rep(1, 500000))

th <-  .5 * x + sqrt(4 * .25 /.75) * g + rnorm(1000000, 0, sqrt(.75)) 
cor(x, th) ^ 2
cor(g, th) ^ 2
cor(g, x)

g <- factor(g, levels = c(0,1))
results <- lm(th ~ g  + x)
summary(results)
modelEffectSizes(results)

library(heplots)
results <- aov(th ~ g  + x)

####################### trial: two stage
n <- 1000000
x <- rnorm(n, 0, 1)
g <- rbinom(n, size = 1, prob = .5)

data <- cbind(x, g)

th1 <- .5 * data[g==0, 1] + rnorm(sum(g==0), 0, sqrt(.75))
th2 <-  .5 * data[g==1, 1] + sqrt(4 * .25 /.75) * data[g==1, 2] + rnorm(sum(g==1), 0, sqrt(.5))

th <- c(th1, th2)
th[g==0] <- th1
th[g==1] <- th2

cor(x, th) ^ 2
cor(g, th) ^ 2
results <- lm(th ~ g  + x)
summary(results)
modelEffectSizes(results)

###################### trial 2

n <- 1000000
x <- rnorm(n, 0, 1)
g <- rbinom(n, size = 1, prob = .5)

th <-  .5 * x + sqrt(4 * .25 /.75) * g + rnorm(1000000, 0, sqrt(.75 )) 
cor(x, th) ^ 2
cor(g, th) ^ 2
results <- lm(th ~ g  + x)
summary(results)
modelEffectSizes(results)


######################  trial 3  CORRECT! 

n <- 10000000
x <- rnorm(n, 0, 1)
g <- rbinom(n, size = 1, prob = .5)

b1 <- sqrt(.1)
b2 <- sqrt(.1/.25)

th <- b1 * x + b2 * g + rnorm(n, 0, sqrt(.8))
cor(x, th) ^ 2
cor(g, th) ^ 2

results <- lm(th ~ g  + x)
summary(results)
library(lmSupport)
modelEffectSizes(results)

#################### trial 3 extra  CORRECT!
#var_th -> .2
#var_z -> .2*.8, thus 80% of variance is from noise
#var_x -> .2*.1, thus 10%
#var_g -> .2*.1, thus 10%

n <- 1000000
x <- rnorm(n, 0, 1)
g <- rbinom(n, size = 1, prob = .5)

b1 <- sqrt(.02)
b2 <- sqrt(.02/.25)

th <- b1 * x + b2 * g + rnorm(n, 0, sqrt(.8*.2))
cor(x, th) ^ 2
cor(g, th) ^ 2

results <- lm(th ~ g  + x)
summary(results)
library(lmSupport)
modelEffectSizes(results)

#################### trial 3 extra final  CORRECT!!!
#var_th -> .2
#var_z -> .2*.8, thus 80% of variance is from noise
#var_x -> .2*.1, thus 10%
#var_g -> .2*.1, thus 10%

n <- 1000000
x <- runif(n, min = 4, max = 12)
g <- rbinom(n, size = 1, prob = .5)

b1 <- sqrt(.02 / (1/12*(12-4)^2) )
b2 <- sqrt(.02/.25)

th <- b1 * x + b2 * g + rnorm(n, 0, sqrt(.8*.2))
cor(x, th) ^ 2
cor(g, th) ^ 2
var(th)
results <- lm(th ~ g  + x)
summary(results)
library(lmSupport)
modelEffectSizes(results)
