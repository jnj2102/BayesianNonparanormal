########################
#Code to run the Frequentist Nonparanormal Truncation
#For the Simulation section of my paper
#Cholesky Decomposition
#Author: Jami Jackson Mulgrave
#
########################

current_dir <- getwd()
setwd(current_dir)  #set the working directory

#clear the workspace
rm(list = ls())

library(R.matlab)
library(huge)

source("Frequentist_nonparanormal.R") #call the function


reps <- 100

###########Simulation combination: n=50, p=25, sparsity = circle#######
result_list <- list()

lambda <- seq(from = 0.16, to = 1.2, length.out = 50) #in the Liu et al 2009 paper

data <- readMat('BayesNonpar_p25_n50_circle_prior.mat', package = "R.matlab") #Read in all the data


for (iters in 1:reps) {
  cat('Iteration = ', iters)

set.seed(100 + iters)
  
#Run this in a loop for each iteration 
  
  #read in the Sigma_true 
  
  sigma_true <- data$sigma.true
  

  #read in the omega_true
  
  omega_true <- data$omega.true
  
  #read in the x matrix (observed data)
  
  xmat <- data$x.matrix.n50.p25[[iters]][[1]]
  
  p = ncol(xmat)
  n = nrow(xmat)
  
  
  
  #True edge matrix for the model 
  edge_matrix_true <- matrix(0, nrow = p, ncol = p)
  
  for (j in 1:p) {
    for (k in 1:p) {
      
      if ( omega_true[j,k] == 0) {
        edge_matrix_true[j,k] = 0
      } else {
        edge_matrix_true[j,k] = 1
      }
      
    }
  }
  
  #find the upper indices
  indmx = matrix(1:p^2, nrow = p, ncol = p)
  
  upperind_diag = which(triu(indmx)>0)  #include the diagonal
  
  upperind = which(triu(indmx,1)>0)  #do not include the diagonal
  
  
#nonparanormal truncation with graphical lasso
result_list[[iters]] <- Frequentist_nonparanormal(lambda, omega_true, sigma_true, xmat,n,p,edge_matrix_true,upperind)


}

###################################################################

#Save the data for latex later

save.image(file = "Frequentist_p25_n50_circle.rdata")
     
