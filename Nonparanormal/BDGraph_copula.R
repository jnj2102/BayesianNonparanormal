########################
#Code to run the Frequentist Nonparanormal Truncation
#For the Simulation section of my paper
#Cholesky Decomposition
#Author: Jami Jackson Mulgrave
#
########################

BDGraph_copula <- function(omega_true, sigma_true, xmat,p) {
  
  
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
  

ptm <- proc.time()

bdgraph.obj <- bdgraph(xmat, method = "gcgm", iter = 15000, burnin = 5000, multi.update = 1,
                       save.all = FALSE) #Gaussian copula graphical models
#10000 MCMCs after a burnin of 5000 using  multi.update = 1

algorithm_time <- proc.time() - ptm  #I would use the user value for the algorithm time.

#posterior estimate of the precision 

precisionMat_estimate = bdgraph.obj$K_hat


#Estimated edge matrix for the model 
edge_matrix_selected  <- select(bdgraph.obj) #selected graph using Bayesian model averaging

#true correlation matrix and true inverse correlation matrix
sigma_true_correlation = cov2cor(sigma_true)
inverse_sigma_true_correlation = solve(sigma_true_correlation)


#trace is the sum of the diagonal elements. 

#Use the optimal precision matrix for the entropy loss

entropy_loss_correlation = sum(diag(precisionMat_estimate%*% sigma_true_correlation)) -
  log(det(precisionMat_estimate%*%sigma_true_correlation)) - p;



#Frobenius loss
FrobeniusLoss_precision = (sum(diag(crossprod(precisionMat_estimate - inverse_sigma_true_correlation)))) 

#Bounded loss
bounded_loss = 1/(p^2) * sum(sum(abs(precisionMat_estimate - inverse_sigma_true_correlation)))



#Find the TP, TN, FP, and FN
TP_matrix = edge_matrix_selected[upperind] == 1 & edge_matrix_true[upperind] == 1

TP = sum(TP_matrix) #the colon sums all elements in the matrix


TN_matrix = edge_matrix_selected[upperind] == 0 & edge_matrix_true[upperind] == 0

TN = sum(TN_matrix) #the colon sums all elements in the matrix


FP_matrix = edge_matrix_selected[upperind] == 1 & edge_matrix_true[upperind] == 0
FP = sum(FP_matrix) #the colon sums all elements in the matrix


FN_matrix = edge_matrix_selected[upperind] == 0 & edge_matrix_true[upperind] == 1
FN = sum(FN_matrix) #the colon sums all elements in the matrix

#Find Specificity, Sensitivity, Matthews Correlation Coefficient

SP = TN/(TN + FP)

SE= TP/(TP + FN)

temp1 = as.numeric((TP + FP)*(TP + FN))
temp2 = as.numeric((TN + FP)*(TN + FN))

MCC = as.numeric(((TP * TN) - (FP * FN))/sqrt(temp1*temp2))

return(list(algorithm_time = algorithm_time,
            TP=TP,
            TN = TN, FP = FP, FN = FN,SP = SP,
            SE = SE, SE= SE,
            MCC = MCC,
            bdgraph.obj = bdgraph.obj,
            entropy_loss_correlation = entropy_loss_correlation,
      			 edge_matrix_selected = edge_matrix_selected,
      			 precisionMat_estimate = precisionMat_estimate,
      			 FrobeniusLoss_precision = FrobeniusLoss_precision,
      			 bounded_loss = bounded_loss))
}
