########################
#Code to run the Frequentist Nonparanormal Truncation
#For the Simulation section of my paper
#Author: Jami Jackson Mulgrave
#
########################

current_dir <- getwd()
setwd(current_dir)  #set the working directory

#clear the workspace
rm(list = ls())

library(R.matlab)
library(huge)


# A function to remove the NAs to find the standard error
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


reps <- 100

###########Simulation combination: n=500, p=100, sparsity = circle#######


L1_freq_n500_p100_circle_ebic <- c()
SP_freq_n500_p100_circle_ebic <- c()
SE_freq_n500_p100_circle_ebic <- c()
MCC_freq_n500_p100_circle_ebic <- c()
ans_final_ebic <- list()
TP_ebic <- c()
TN_ebic <- c()
FP_ebic <- c()
FN_ebic <- c()

L1_freq_n500_p100_circle_stars <- c()
SP_freq_n500_p100_circle_stars <- c()
SE_freq_n500_p100_circle_stars <- c()
MCC_freq_n500_p100_circle_stars <- c()
ans_final_stars <- list()
TP_stars <- c()
TN_stars <- c()
FP_stars <- c()
FN_stars <- c()

L1_freq_n500_p100_circle_stars_alternative <- c()
SP_freq_n500_p100_circle_stars_alternative <- c()
SE_freq_n500_p100_circle_stars_alternative <- c()
MCC_freq_n500_p100_circle_stars_alternative <- c()
ans_final_stars_alternative <- list()
TP_stars_alternative <- c()
TN_stars_alternative <- c()
FP_stars_alternative <- c()
FN_stars_alternative <- c()

L1_freq_n500_p100_circle_ric <- c()
SP_freq_n500_p100_circle_ric <- c()
SE_freq_n500_p100_circle_ric <- c()
MCC_freq_n500_p100_circle_ric <- c()
ans_final_ric <- list()
TP_ric <- c()
TN_ric <- c()
FP_ric <- c()
FN_ric <- c()

algorithm_time <- list()
out_npn <- list()
Ymat_npn <- list()


lambda <- seq(from = 0.16, to = 1.2, length.out = 50) #in the Liu et al 2009 paper



for (iters in 1:reps) {
  cat('Iteration = ', iters)

set.seed(iters)
  
#Run this in a loop for each iteration 
  
  #read in the Sigma_true 
  
  sigma_true <- readMat('Sigma_true_n500_p100_circle_1to101.mat',
                        package = "R.matlab")$sigma.true
  

  
  
  #read in the omega_true
  
  omega_true <- readMat('Omega_true_n500_p100_circle_1to101.mat', package = R.matlab)$omega.true
  
  p = dim(omega_true)[1]
  
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
  
  
  #read in the x matrix (observed data)
  
xmat <- readMat(sprintf("Bsplines_Iter_%d_x_matrix_n500_p100_circle_1to100.mat", iters), 
                package = "R.matlab")$x.matrix

ptm <- proc.time()

Ymat_npn[[iters]] = huge.npn(xmat, npn.func="truncation") # Nonparanormal


out_npn[[iters]] = huge(Ymat_npn[[iters]],lambda = lambda, nlambda = length(lambda),
               method = "glasso", cov.output = TRUE) 

algorithm_time[[iters]] <- proc.time() - ptm  #I would use the user value for the algorithm time.

### Use EBIC for Model Selection #####

ans_ebic = huge.select(out_npn[[iters]], criterion = "ebic") #said to use ebic for
#glasso in vignette

ans_final_ebic[[iters]] = ans_ebic

#optimal precision and covariance matrix

opt_precisionMat_ebic = ans_ebic$opt.icov
opt_covarianceMat_ebic = ans_ebic$opt.cov


#Estimated edge matrix for the model 
edge_matrix_est_ebic <- ans_ebic$refit

#Find the L1 loss function using the optimal covariance matrix

#trace is the sum of the diagonal elements.  I am using the covariance
#matrix now.

L1_freq_n500_p100_circle_ebic[iters] = sum(diag(solve(opt_covarianceMat_ebic) %*% sigma_true)) -
  log(det(solve(opt_covarianceMat_ebic) %*% sigma_true)) - p


#Find the TP, TN, FP, and FN
TP_matrix_ebic = edge_matrix_est_ebic[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_ebic[iters] = sum(TP_matrix_ebic) #the colon sums all elements in the matrix


TN_matrix_ebic = edge_matrix_est_ebic[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_ebic[iters] = sum(TN_matrix_ebic) #the colon sums all elements in the matrix


FP_matrix_ebic = edge_matrix_est_ebic[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_ebic[iters] = sum(FP_matrix_ebic) #the colon sums all elements in the matrix


FN_matrix_ebic = edge_matrix_est_ebic[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_ebic[iters] = sum(FN_matrix_ebic) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_n500_p100_circle_ebic[iters] = TN_ebic[iters]/(TN_ebic[iters] + FP_ebic[iters])

SE_freq_n500_p100_circle_ebic[iters] = TP_ebic[iters]/(TP_ebic[iters] + FN_ebic[iters])

temp1_ebic = as.numeric((TP_ebic[iters] + FP_ebic[iters])*(TP_ebic[iters] + FN_ebic[iters]))
temp2_ebic = as.numeric((TN_ebic[iters] + FP_ebic[iters])*(TN_ebic[iters] + FN_ebic[iters]))

MCC_freq_n500_p100_circle_ebic[iters] = ((TP_ebic[iters] * TN_ebic[iters]) -
                                 (FP_ebic[iters] * FN_ebic[iters]))/sqrt(temp1_ebic*temp2_ebic)

### Use STARS for Model Selection using Threshold = .1 (default) #####

ans_stars = huge.select(out_npn[[iters]], criterion = "stars") 

ans_final_stars[[iters]] = ans_stars

#optimal precision and covariance matrix

opt_precisionMat_stars = ans_stars$opt.icov
opt_covarianceMat_stars = ans_stars$opt.cov


#Estimated edge matrix for the model 
edge_matrix_est_stars <- ans_stars$refit

#Find the L1 loss function using the optimal covariance matrix

#trace is the sum of the diagonal elements.  I am using the covariance
#matrix now.

L1_freq_n500_p100_circle_stars[iters] = sum(diag(solve(opt_covarianceMat_stars) %*% sigma_true)) -
  log(det(solve(opt_covarianceMat_stars) %*% sigma_true)) - p


#Find the TP, TN, FP, and FN
TP_matrix_stars = edge_matrix_est_stars[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_stars[iters] = sum(TP_matrix_stars) #the colon sums all elements in the matrix


TN_matrix_stars = edge_matrix_est_stars[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_stars[iters] = sum(TN_matrix_stars) #the colon sums all elements in the matrix


FP_matrix_stars = edge_matrix_est_stars[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_stars[iters] = sum(FP_matrix_stars) #the colon sums all elements in the matrix


FN_matrix_stars = edge_matrix_est_stars[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_stars[iters] = sum(FN_matrix_stars) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_n500_p100_circle_stars[iters] = TN_stars[iters]/(TN_stars[iters] + FP_stars[iters])

SE_freq_n500_p100_circle_stars[iters] = TP_stars[iters]/(TP_stars[iters] + FN_stars[iters])

temp1_stars = as.numeric((TP_stars[iters] + FP_stars[iters])*(TP_stars[iters] + FN_stars[iters]))
temp2_stars = as.numeric((TN_stars[iters] + FP_stars[iters])*(TN_stars[iters] + FN_stars[iters]))

MCC_freq_n500_p100_circle_stars[iters] = ((TP_stars[iters] * TN_stars[iters]) -
                                              (FP_stars[iters] * FN_stars[iters]))/sqrt(temp1_stars*temp2_stars)

### Use STARS for Model Selection using Threshold = .05 (used in STARS paper) #####

ans_stars_alternative = huge.select(out_npn[[iters]], criterion = "stars", stars.thresh = .05) 

ans_final_stars_alternative[[iters]] = ans_stars_alternative

#optimal precision and covariance matrix

opt_precisionMat_stars_alternative = ans_stars_alternative$opt.icov
opt_covarianceMat_stars_alternative = ans_stars_alternative$opt.cov


#Estimated edge matrix for the model 
edge_matrix_est_stars_alternative <- ans_stars_alternative$refit

#Find the L1 loss function using the optimal covariance matrix

#trace is the sum of the diagonal elements.  I am using the covariance
#matrix now.

L1_freq_n500_p100_circle_stars_alternative[iters] = sum(diag(solve(opt_covarianceMat_stars_alternative) %*% sigma_true)) -
  log(det(solve(opt_covarianceMat_stars_alternative) %*% sigma_true)) - p


#Find the TP, TN, FP, and FN
TP_matrix_stars_alternative = edge_matrix_est_stars_alternative[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_stars_alternative[iters] = sum(TP_matrix_stars_alternative) #the colon sums all elements in the matrix


TN_matrix_stars_alternative = edge_matrix_est_stars_alternative[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_stars_alternative[iters] = sum(TN_matrix_stars_alternative) #the colon sums all elements in the matrix


FP_matrix_stars_alternative = edge_matrix_est_stars_alternative[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_stars_alternative[iters] = sum(FP_matrix_stars_alternative) #the colon sums all elements in the matrix


FN_matrix_stars_alternative = edge_matrix_est_stars_alternative[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_stars_alternative[iters] = sum(FN_matrix_stars_alternative) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_n500_p100_circle_stars_alternative[iters] = TN_stars_alternative[iters]/(TN_stars_alternative[iters] + FP_stars_alternative[iters])

SE_freq_n500_p100_circle_stars_alternative[iters] = TP_stars_alternative[iters]/(TP_stars_alternative[iters] + FN_stars_alternative[iters])

temp1_stars_alternative = as.numeric((TP_stars_alternative[iters] + FP_stars_alternative[iters])*(TP_stars_alternative[iters] + FN_stars_alternative[iters]))
temp2_stars_alternative = as.numeric((TN_stars_alternative[iters] + FP_stars_alternative[iters])*(TN_stars_alternative[iters] + FN_stars_alternative[iters]))

MCC_freq_n500_p100_circle_stars_alternative[iters] = ((TP_stars_alternative[iters] * TN_stars_alternative[iters]) -
                                                          (FP_stars_alternative[iters] * FN_stars_alternative[iters]))/sqrt(temp1_stars_alternative*temp2_stars_alternative)



### Use STARS for Model Selection using RIC #####

ans_ric = huge.select(out_npn[[iters]], criterion = "ric") 

ans_final_ric[[iters]] = ans_ric

#optimal precision and covariance matrix

opt_precisionMat_ric = ans_ric$opt.icov
opt_covarianceMat_ric = ans_ric$opt.cov


#Estimated edge matrix for the model 
edge_matrix_est_ric <- ans_ric$refit

#Find the L1 loss function using the optimal covariance matrix

#trace is the sum of the diagonal elements.  I am using the covariance
#matrix now.

L1_freq_n500_p100_circle_ric[iters] = sum(diag(solve(opt_covarianceMat_ric) %*% sigma_true)) -
  log(det(solve(opt_covarianceMat_ric) %*% sigma_true)) - p


#Find the TP, TN, FP, and FN
TP_matrix_ric = edge_matrix_est_ric[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_ric[iters] = sum(TP_matrix_ric) #the colon sums all elements in the matrix


TN_matrix_ric = edge_matrix_est_ric[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_ric[iters] = sum(TN_matrix_ric) #the colon sums all elements in the matrix


FP_matrix_ric = edge_matrix_est_ric[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_ric[iters] = sum(FP_matrix_ric) #the colon sums all elements in the matrix


FN_matrix_ric = edge_matrix_est_ric[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_ric[iters] = sum(FN_matrix_ric) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_n500_p100_circle_ric[iters] = TN_ric[iters]/(TN_ric[iters] + FP_ric[iters])

SE_freq_n500_p100_circle_ric[iters] = TP_ric[iters]/(TP_ric[iters] + FN_ric[iters])

temp1_ric = as.numeric((TP_ric[iters] + FP_ric[iters])*(TP_ric[iters] + FN_ric[iters]))
temp2_ric = as.numeric((TN_ric[iters] + FP_ric[iters])*(TN_ric[iters] + FN_ric[iters]))

MCC_freq_n500_p100_circle_ric[iters] = ((TP_ric[iters] * TN_ric[iters]) -
                                            (FP_ric[iters] * FN_ric[iters]))/sqrt(temp1_ric*temp2_ric)
}

###  For EBIC Model Selection ####

#Find the means 
mean_L1_freq_n500_p100_circle_ebic = mean(L1_freq_n500_p100_circle_ebic)
mean_SP_freq_n500_p100_circle_ebic = mean(SP_freq_n500_p100_circle_ebic)
mean_SE_freq_n500_p100_circle_ebic = mean(SE_freq_n500_p100_circle_ebic)

#I removed the NAs for MCC

mean_MCC_freq_n500_p100_circle_ebic = mean(MCC_freq_n500_p100_circle_ebic, na.rm = TRUE)  

#Find the standard errors

stderr_L1_freq_n500_p100_circle_ebic = stderr(L1_freq_n500_p100_circle_ebic)
stderr_SP_freq_n500_p100_circle_ebic = stderr(SP_freq_n500_p100_circle_ebic)
stderr_SE_freq_n500_p100_circle_ebic = stderr(SE_freq_n500_p100_circle_ebic)
stderr_MCC_freq_n500_p100_circle_ebic = stderr(MCC_freq_n500_p100_circle_ebic)

###  For STARS Model Selection using Threshold = .1####

#Find the means 
mean_L1_freq_n500_p100_circle_stars = mean(L1_freq_n500_p100_circle_stars)
mean_SP_freq_n500_p100_circle_stars = mean(SP_freq_n500_p100_circle_stars)
mean_SE_freq_n500_p100_circle_stars = mean(SE_freq_n500_p100_circle_stars)

#I removed the NAs for MCC

mean_MCC_freq_n500_p100_circle_stars = mean(MCC_freq_n500_p100_circle_stars, na.rm = TRUE)  

#Find the standard errors

stderr_L1_freq_n500_p100_circle_stars = stderr(L1_freq_n500_p100_circle_stars)
stderr_SP_freq_n500_p100_circle_stars = stderr(SP_freq_n500_p100_circle_stars)
stderr_SE_freq_n500_p100_circle_stars = stderr(SE_freq_n500_p100_circle_stars)
stderr_MCC_freq_n500_p100_circle_stars = stderr(MCC_freq_n500_p100_circle_stars)


###  For STARS Model Selection using Threshold = .05####

#Find the means 
mean_L1_freq_n500_p100_circle_stars_alternative = mean(L1_freq_n500_p100_circle_stars_alternative)
mean_SP_freq_n500_p100_circle_stars_alternative = mean(SP_freq_n500_p100_circle_stars_alternative)
mean_SE_freq_n500_p100_circle_stars_alternative = mean(SE_freq_n500_p100_circle_stars_alternative)

#I removed the NAs for MCC

mean_MCC_freq_n500_p100_circle_stars_alternative = mean(MCC_freq_n500_p100_circle_stars_alternative, na.rm = TRUE)  

#Find the standard errors

stderr_L1_freq_n500_p100_circle_stars_alternative = stderr(L1_freq_n500_p100_circle_stars_alternative)
stderr_SP_freq_n500_p100_circle_stars_alternative = stderr(SP_freq_n500_p100_circle_stars_alternative)
stderr_SE_freq_n500_p100_circle_stars_alternative = stderr(SE_freq_n500_p100_circle_stars_alternative)
stderr_MCC_freq_n500_p100_circle_stars_alternative = stderr(MCC_freq_n500_p100_circle_stars_alternative)

###  For RIC Model Selection ####

#Find the means 
mean_L1_freq_n500_p100_circle_ric = mean(L1_freq_n500_p100_circle_ric)
mean_SP_freq_n500_p100_circle_ric = mean(SP_freq_n500_p100_circle_ric)
mean_SE_freq_n500_p100_circle_ric = mean(SE_freq_n500_p100_circle_ric)

#I removed the NAs for MCC

mean_MCC_freq_n500_p100_circle_ric = mean(MCC_freq_n500_p100_circle_ric, na.rm = TRUE)  

#Find the standard errors

stderr_L1_freq_n500_p100_circle_ric = stderr(L1_freq_n500_p100_circle_ric)
stderr_SP_freq_n500_p100_circle_ric = stderr(SP_freq_n500_p100_circle_ric)
stderr_SE_freq_n500_p100_circle_ric = stderr(SE_freq_n500_p100_circle_ric)
stderr_MCC_freq_n500_p100_circle_ric = stderr(MCC_freq_n500_p100_circle_ric)

###################################################################

#Save the data for latex later

save.image(file = "Bsplines_paper_frequentist_p100_n500_circle.rdata")
     
