####################################################################
######  Code for checking reproducability  #########################
####################################################################


### READ ME ########################################################
# This code is for checking whether the code is reproducable.
# When the code for simulation was finished, Zhengguo simulated 
# a few toy examples to see whether the results were exactly the same
# and yes, they were. But to make sure that simulated datasets were
# reproducable, The code for simulating data were executed twice 
# (see the explanation files in the datapackage. )
#
# Thus, this code is for checking whether the simulated data were
# exactly the same.
#
# RESULT: the data are exactly the same (tested on 31/03/2017).
# This code is execuated on Zhengguo's PC at home: Windows 10 Home, 
# processor: Intel Core i7-5820K, 3.30 GHz, RAM: 32GB, 
# 64-bit operating system, x64 based processor.
####################################################################

num_test <- 1
while(num_test <= 360){
  
  filename <- paste("D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170330 simresult/phdpro2BladZhengguo/results_", num_test, ".RData", sep = "")
  load(filename)
  sim_result1st <- sim_result
  
  itemparamter <- paste("D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170330 simresult/phdpro2BladZhengguo/item_", num_test, ".RData", sep = "")
  load(itemparamter)
  itemparamter1st <- item_par
  
  
  filename <- paste("D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170330 simresult 2nd/PhdProj2Zhengguo/results_", num_test, ".RData", sep = "")
  load(filename)
  sim_result2nd <- sim_result
  
  itemparamter <- paste("D:/Dropbox/Tilburg office/Research Individual change/Project 2 - norming change/20170330 simresult 2nd/PhdProj2Zhengguo/item_", num_test, ".RData", sep = "")
  load(itemparamter)
  itemparamter2nd <- item_par
  
  if(identical(sim_result1st, sim_result2nd) == TRUE  & identical(itemparamter1st, itemparamter2nd) == TRUE) {
    rm(sim_result, sim_result1st, sim_result2nd, filename, item_par, itemparamter, itemparamter1st, itemparamter2nd, X1, X2) # clean environment and make room for the next num_test
    print(num_test)
  } else{
    print("HORROR! Data are not the same!")
    stop()
  }
  
  num_test <- num_test + 1 
}
