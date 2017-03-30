##############################################################
#########                                        #############
######### Data Analysis: Norming change scores   #############
#########                                        #############
#########      Zhengguo Gu, Tilburg University   #############
#########      Last update: 30/03/2017           #############
##############################################################


# 1. reorganize the data.
# To save space, for each test condition (in total 360), pretest and posttest
# scores of 1000 datasets are 'cbind'ed. For example, the results of the first 
# test condition:
load("D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170330 simresult/phdpro2BladZhengguo/results_1.RData")

dim(sim_result)  #100 x 2000 matrix. 100: 100 persons. 2000: 1000 datasets, each contain a pretest and posttest
                 #score, and thus 1000x2=2000

# 2. standardize scores. 
mydata <- scale(sim_result, center = TRUE, scale = TRUE)

 # identical((sim_result[, 22] - mean(sim_result[, 22]))/sd(sim_result[, 22]), mydata[,22])

##############################################################
###### norming methods
##############################################################
# 3. The regression-based change approach
