########################
#Code to run the Frequentist Nonparanormal Truncation
#For the Simulation section of my paper
#Cholesky Decomposition
#Author: Jami Jackson Mulgrave
#
########################

Frequentist_nonparanormal <- function(lambda, omega_true, sigma_true, xmat,n,p,edge_matrix_true,upperind) {
  
ptm <- proc.time()

Ymat_npn = huge.npn(xmat, npn.func="truncation") # Nonparanormal


out_npn = huge(Ymat_npn,lambda = lambda, nlambda = length(lambda),
               method = "glasso", cov.output = TRUE) 

algorithm_time <- proc.time() - ptm  #I would use the user value for the algorithm time.

### Use EBIC for Model Selection #####

ans_ebic = huge.select(out_npn, criterion = "ebic") #said to use ebic for
#glasso in vignette

#ans_final_ebic = ans_ebic

#optimal precision and covariance matrix

opt_precisionMat_ebic = ans_ebic$opt.icov
opt_covarianceMat_ebic = ans_ebic$opt.cov


#Estimated edge matrix for the model 
edge_matrix_est_ebic <- ans_ebic$refit

# #trace is the sum of the diagonal elements. 
# #Use the optimal precision matrix for the entropy loss

entropy_loss_ebic  = sum(diag(opt_precisionMat_ebic%*% sigma_true)) -
  log(det(opt_precisionMat_ebic%*%sigma_true)) - p;


#Frobenius loss
FrobeniusLoss_precision_ebic = (sum(diag(crossprod(opt_precisionMat_ebic - omega_true)))) 


FrobeniusLoss_covariance_ebic = (sum(diag(crossprod(opt_covarianceMat_ebic - sigma_true)))) 

#Bounded loss
bounded_loss_ebic = 1/(p^2) * sum(sum(abs(opt_precisionMat_ebic - omega_true)))


#Find the TP, TN, FP, and FN
TP_matrix_ebic = edge_matrix_est_ebic[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_ebic = sum(TP_matrix_ebic) #the colon sums all elements in the matrix


TN_matrix_ebic = edge_matrix_est_ebic[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_ebic = sum(TN_matrix_ebic) #the colon sums all elements in the matrix


FP_matrix_ebic = edge_matrix_est_ebic[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_ebic = sum(FP_matrix_ebic) #the colon sums all elements in the matrix


FN_matrix_ebic = edge_matrix_est_ebic[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_ebic = sum(FN_matrix_ebic) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_ebic = TN_ebic/(TN_ebic + FP_ebic)

SE_freq_ebic= TP_ebic/(TP_ebic + FN_ebic)

temp1_ebic = as.numeric((TP_ebic + FP_ebic)*(TP_ebic + FN_ebic))
temp2_ebic = as.numeric((TN_ebic + FP_ebic)*(TN_ebic + FN_ebic))

MCC_freq_ebic = as.numeric(((TP_ebic * TN_ebic) - (FP_ebic * FN_ebic))/sqrt(temp1_ebic*temp2_ebic))

### Use STARS for Model Selection using Threshold = .1 (default) #####

ans_stars = huge.select(out_npn, criterion = "stars") 


#optimal precision and covariance matrix

opt_precisionMat_stars = ans_stars$opt.icov
opt_covarianceMat_stars = ans_stars$opt.cov


#Estimated edge matrix for the model 
edge_matrix_est_stars <- ans_stars$refit


#trace is the sum of the diagonal elements. 


#Use the optimal precision matrix for the entropy loss

entropy_loss_stars  = sum(diag(opt_precisionMat_stars%*% sigma_true)) -
  log(det(opt_precisionMat_stars%*%sigma_true)) - p;


#Frobenius loss
FrobeniusLoss_precision_stars = (sum(diag(crossprod(opt_precisionMat_stars - omega_true)))) 


FrobeniusLoss_covariance_stars = (sum(diag(crossprod(opt_covarianceMat_stars - sigma_true)))) 

#Bounded loss
bounded_loss_stars = 1/(p^2) * sum(sum(abs(opt_precisionMat_stars - omega_true)))

#Find the TP, TN, FP, and FN
TP_matrix_stars = edge_matrix_est_stars[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_stars = sum(TP_matrix_stars) #the colon sums all elements in the matrix


TN_matrix_stars = edge_matrix_est_stars[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_stars = sum(TN_matrix_stars) #the colon sums all elements in the matrix


FP_matrix_stars = edge_matrix_est_stars[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_stars = sum(FP_matrix_stars) #the colon sums all elements in the matrix


FN_matrix_stars = edge_matrix_est_stars[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_stars = sum(FN_matrix_stars) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_stars = TN_stars/(TN_stars + FP_stars)

SE_freq_stars = TP_stars/(TP_stars + FN_stars)

temp1_stars = as.numeric((TP_stars + FP_stars)*(TP_stars + FN_stars))
temp2_stars = as.numeric((TN_stars + FP_stars)*(TN_stars + FN_stars))

MCC_freq_stars = as.numeric(((TP_stars * TN_stars) -(FP_stars * FN_stars))/sqrt(temp1_stars*temp2_stars))



### Use RIC for Model Selection #####

ans_ric = huge.select(out_npn, criterion = "ric") 


#optimal precision and covariance matrix

opt_precisionMat_ric = ans_ric$opt.icov
opt_covarianceMat_ric = ans_ric$opt.cov


#Estimated edge matrix for the model 
edge_matrix_est_ric <- ans_ric$refit


#trace is the sum of the diagonal elements.  

#Use the optimal precision matrix for the entropy loss

entropy_loss_ric  = sum(diag(opt_precisionMat_ric%*% sigma_true)) -
  log(det(opt_precisionMat_ric%*%sigma_true)) - p;


#Frobenius loss
FrobeniusLoss_precision_ric = (sum(diag(crossprod(opt_precisionMat_ric - omega_true)))) 


FrobeniusLoss_covariance_ric = (sum(diag(crossprod(opt_covarianceMat_ric - sigma_true)))) 

#Bounded loss
bounded_loss_ric = 1/(p^2) * sum(sum(abs(opt_precisionMat_ric - omega_true)))

#Find the TP, TN, FP, and FN
TP_matrix_ric = edge_matrix_est_ric[upperind] == 1 & edge_matrix_true[upperind] == 1

TP_ric = sum(TP_matrix_ric) #the colon sums all elements in the matrix


TN_matrix_ric = edge_matrix_est_ric[upperind] == 0 & edge_matrix_true[upperind] == 0

TN_ric = sum(TN_matrix_ric) #the colon sums all elements in the matrix


FP_matrix_ric = edge_matrix_est_ric[upperind] == 1 & edge_matrix_true[upperind] == 0
FP_ric = sum(FP_matrix_ric) #the colon sums all elements in the matrix


FN_matrix_ric = edge_matrix_est_ric[upperind] == 0 & edge_matrix_true[upperind] == 1
FN_ric = sum(FN_matrix_ric) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP_freq_ric = TN_ric/(TN_ric + FP_ric)

SE_freq_ric = TP_ric/(TP_ric + FN_ric)

temp1_ric = as.numeric((TP_ric + FP_ric)*(TP_ric + FN_ric))
temp2_ric = as.numeric((TN_ric + FP_ric)*(TN_ric + FN_ric))

MCC_freq_ric = as.numeric(((TP_ric * TN_ric) - (FP_ric * FN_ric))/sqrt(temp1_ric*temp2_ric))



return(list(algorithm_time = algorithm_time,ans_ebic = ans_ebic,
            TP_ebic=TP_ebic,
            TN_ebic = TN_ebic, FP_ebic = FP_ebic, FN_ebic = FN_ebic,SP_freq_ebic = SP_freq_ebic,
            SE_freq_ebic = SE_freq_ebic, SE_freq_ebic= SE_freq_ebic,
            MCC_freq_ebic = MCC_freq_ebic,
            ans_stars = ans_stars,
            TP_stars=TP_stars,
            TN_stars = TN_stars, FP_stars = FP_stars, FN_stars = FN_stars,SP_freq_stars = SP_freq_stars,
            SE_freq_stars = SE_freq_stars, SE_freq_stars= SE_freq_stars,
            MCC_freq_stars = MCC_freq_stars,
            ans_ric = ans_ric,
            TP_ric=TP_ric,
            TN_ric = TN_ric, FP_ric = FP_ric, FN_ric = FN_ric,SP_freq_ric = SP_freq_ric,
            SE_freq_ric = SE_freq_ric, SE_freq_ric= SE_freq_ric,
            MCC_freq_ric = MCC_freq_ric,
			entropy_loss_stars = entropy_loss_stars,
FrobeniusLoss_precision_stars = FrobeniusLoss_precision_stars,
FrobeniusLoss_covariance_stars = FrobeniusLoss_covariance_stars,
bounded_loss_stars = bounded_loss_stars,
entropy_loss_ebic = entropy_loss_ebic,
FrobeniusLoss_precision_ebic = FrobeniusLoss_precision_ebic,
FrobeniusLoss_covariance_ebic = FrobeniusLoss_covariance_ebic,
bounded_loss_ebic = bounded_loss_ebic,
entropy_loss_ric = entropy_loss_ric,
FrobeniusLoss_precision_ric  = FrobeniusLoss_precision_ric,
FrobeniusLoss_covariance_ric = FrobeniusLoss_covariance_ric,
bounded_loss_ric  = bounded_loss_ric,
opt_precisionMat_stars = opt_precisionMat_stars,
opt_covarianceMat_stars = opt_covarianceMat_stars,
edge_matrix_est_stars = edge_matrix_est_stars,
opt_precisionMat_ric = opt_precisionMat_ric,
opt_covarianceMat_ric = opt_covarianceMat_ric,
edge_matrix_est_ric = edge_matrix_est_ric,
opt_precisionMat_ebic = opt_precisionMat_ebic,
opt_covarianceMat_ebic = opt_covarianceMat_ebic,
edge_matrix_est_ebic = edge_matrix_est_ebic))
}
