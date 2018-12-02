##############################################
#Code to create the graph using the nonparanormal
#and the graphical lasso
#Author: Jami Jackson Mulgrave
##############################################




#clear the workspace
rm(list = ls())

library(R.matlab)
library(huge)
source('Frequentist_nonparanormal.R')


data = read.table('gb-2004-5-11-r92-s1_removedrows.txt', header = TRUE) #read in the data

#I'm not reading in the MAT file because R.matlab doesn't read -v7.3 files.

#data_matrix = readMat('Bsplines_paper_realdata_initialdata.mat', package = "R.matlab")$data.matrix
 #Read in all the data)

data_matrix = data[, c(7:124)]

data_matrix_log = t(log(data_matrix))

#I need to flip the matrix to be nxp which is 118x39

                         
n = dim(data_matrix_log)[1]

p =  dim(data_matrix_log)[2]

lambda <- seq(from = 0.16, to = 1.2, length.out = 50) #in the Liu et al 2009 paper

set.seed(100)

omega_true = diag(p)
sigma_true = omega_true
edge_matrix_true = omega_true

 indmx = matrix(1:p^2, nrow = p, ncol = p)
  
  upperind_diag = which(triu(indmx)>0)  #include the diagonal
  
  upperind = which(triu(indmx,1)>0)  #do not include the diagonal

#Nonparanormal graph
result <- Frequentist_nonparanormal(lambda, omega_true, sigma_true, data_matrix_log,n,p,edge_matrix_true,upperind)


writeMat('RealData_Frequentist.mat', edge_matrix_est_stars = result$edge_matrix_est_stars, 
         edge_matrix_est_ebic=  result$edge_matrix_est_ebic , edge_matrix_est_ric = result$edge_matrix_est_ric)


save.image(file = "RealData_Frequentist.rdata")
