########################
#Code to run the BDMCMC method using the BDGraph package
#For the Simulation section of my paper
#Bayesian Regression Approach
#Author: Jami Jackson Mulgrave
#
########################

current_dir <- getwd()
setwd(current_dir)  #set the working directory

#clear the workspace
rm(list = ls())

library(R.matlab)
library(BDgraph)
library(Matrix)


source("BDGraph_copula.R") #call the function


reps <- 100

###########Simulation combination: n=500, p=100, sparsity = AR1#######
result_list <- list()

data <- readMat('BayesNonpar_p100_n500_AR1_prior.mat', package = "R.matlab") #Read in all the data


for (iters in 1:reps) {
  cat('Iteration = ', iters)

set.seed(100 + iters)
  
#Run this in a loop for each iteration 
  
  #read in the Sigma_true 
  
  sigma_true <- data$sigma.true
  

  #read in the omega_true
  
  omega_true <- data$omega.true
  
  p = data$p
  
  #read in the x matrix (observed data)
  
xmat <- data$x.matrix.n500.p100[[iters]][[1]]

#nonparanormal truncation with graphical lasso
result_list[[iters]] <- BDGraph_copula(omega_true, sigma_true, xmat,p)


}

###################################################################

#Save the data for latex later

save.image(file = "BDGraph_p100_n500_AR1.rdata")
     
