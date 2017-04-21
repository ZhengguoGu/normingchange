##### Toy example to show that RCI is not qualified as a norming method #####
##### Zhengguo Gu Tilburg University 

set.seed(110)

sample_size <- 10000
rel_sample <- c(.5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1)
false_reliaChange <- array()
rel_samp1 <- array()
rel_samp2 <- array()
for(i in 1:length(rel_sample)){
  samp_T1 <- rnorm(sample_size, 0, sqrt(8.45^2*0.91))
  samp_X1 <- samp_T1 + rnorm(sample_size, 0, sqrt(var(samp_T1)/rel_sample[i]-var(samp_T1)))
  rel_samp1[i] <- var(samp_T1)/var(samp_X1)

#suppose no change in True scores.
  samp_X2 <- samp_T1 + rnorm(sample_size, 0, sqrt(var(samp_T1)/rel_sample[i]-var(samp_T1)))
  rel_samp2[i] <- var(samp_T1)/var(samp_X2)
  false_reliaChange[i] <- sum(abs(samp_X2 - samp_X1)>7)/sample_size
  
}

ytick <- paste(c(0,seq(1:10))*10, "%", sep = "")

plot(false_reliaChange, type = "b", xlab = "Reliability at Pretest", ylab = "Proportion of persons falsely identified", xaxt="n", yaxt="n", ylim=c(0,1), pch = 15)
abline(h=c(.05), lty=2)
axis(1, at=1:11,labels=rel_sample, las=1)
axis(2, at=seq(0, 1, by = .1),labels = ytick, las=1)


