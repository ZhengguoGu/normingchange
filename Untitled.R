g <- c(rep(0, 500000), rep(1, 500000))
x <- rnorm(1000000, 0, 1)

th <- .5 * x + rnorm(1000000, 0, sqrt(.75))

cor(x, th)

results <- lm(th ~ x)
summary(results)

th <- sqrt(4 * .25 /.75) * g + rnorm(1000000, 0, sqrt(1))
cor(g, th)^2
results <- lm(th ~ g)
summary(results)

library(lmSupport)
modelEffectSizes(results)

b1 = sqrt(0.1)
b2 = sqrt(.4/.9)      
cvar = 1 -.1 - .1
n = 1000000
x1 = runif(n, min = 4, max = 12)
g = c(rep(0,n/2),rep(1,n/2))

x1 = rnorm(n)

th = b1*scale(x1, center = T, scale = T) + b2*g + rnorm(n,0,sqrt(cvar))
tht = b1*x1 + b2*g
cor(th,tht)^2
cor(th,g)^2
cor(th,x1)^2

results <- lm(th ~ g + x1)
summary(results)
