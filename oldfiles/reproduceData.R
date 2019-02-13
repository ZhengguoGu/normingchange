##########################################################################################
#### This is to double check the reproducability of the Rcode "Dataanalysis.R" 
####
#### Here, we compared the results from "20170402firstrun.RData" and those from 
#### "20170403secondrun.RData". They MUST be identical. (Note that the reason 
#### to check this is because parallel computing is used in "Dataanalysis.R"
#### and we must make sure that the results after parallel computing are assembled
#### in the same way - thus reproducable. As a side, when writing the code 
#### "Dataanalysis.R", Zhengguo already checked the reproducability of the code,
#### and thus this is just to DOUBLE check to make sure. )
#### 
#### The two datasets were checked on Zhengguo's office computer on 2017 - 04 - 05. 

load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170403 dataanalysis/phdproj2Zhengguo/20170402firstrun.RData")
first <- IPR_reg
load(file ="D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170403 dataanalysis/phdproj2Zhengguo/20170403secondrun.RData")
second <- IPR_reg

identical(first, second)  # the result is TRUE. 

