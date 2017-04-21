##### Toy example to show that RCI is not qualified as a norming method #####
##### Zhengguo Gu Tilburg University 

set.seed(110)
sample_size <- 10000000

T1 <- rnorm(sample_size, 0, sqrt(8.45^2*0.9)) 
E1 <- rnorm(sample_size, 0, sqrt(8.45^2*0.1))

X1 <- T1 + E1 
sd(X1)
rel_x <- var(T1)/var(X1)
SE <- sd(X1)*sqrt(1-rel_x )
sdiff <- sqrt(2*SE^2)
RCInorm <- sdiff*1.96


sample_size <- 10000
rel_sample <- c(.5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1)
false_reliaChange <- array()
rel_samp1 <- array()
rel_samp2 <- array()
for(i in 1:length(rel_sample)){
  samp_T1 <- sample(T1, size=sample_size, replace=FALSE)
  samp_X1 <- samp_T1 + rnorm(sample_size, 0, sqrt(var(samp_T1)/rel_sample[i]-var(samp_T1)))
  rel_samp1[i] <- var(samp_T1)/var(samp_X1)

#suppose no change in True scores.
  samp_X2 <- samp_T1 + rnorm(sample_size, 0, sqrt(var(samp_T1)/rel_sample[i]-var(samp_T1)))
  rel_samp2[i] <- var(samp_T1)/var(samp_X2)
  false_reliaChange[i] <- sum(abs(samp_X2 - samp_X1)>7)/sample_size
  
}

ytick <- paste(seq(1:10)*10, "%", sep = "")

plot(false_reliaChange, type = "b", xlab = "Reliability at Pretest", ylab = "Proportion of persons falsely identified", xaxt="n", yaxt="n", ylim=c(0,1), pch = 15)
axis(1, at=1:11,labels=rel_sample, las=1)
axis(2, at=seq(.1, 1, by = .1),labels = ytick, las=1)

