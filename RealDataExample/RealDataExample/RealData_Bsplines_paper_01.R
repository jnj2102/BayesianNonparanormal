##############################################
#Code to create the graph using the nonparanormal
#and the graphical lasso
#Author: Jami Jackson Mulgrave
##############################################




#clear the workspace
rm(list = ls())

library(R.matlab)
library(huge)


data = read.table('gb-2004-5-11-r92-s1_removedrows.txt', header = TRUE) #read in the data


data_matrix = data[, c(7:124)]

data_matrix_log = t(log(data_matrix))

#I need to flip the matrix to be nxp which is 118x39

#data_matrix_std = scale(t(data_matrix_log), center = TRUE, scale = TRUE)

data_matrix_std = (data_matrix_log - min(data_matrix_log))/(max(data_matrix_log) - min(data_matrix_log))

                         
n = dim(data_matrix_std)[1]

p =  dim(data_matrix_std)[2]

lambda <- seq(from = 0.16, to = 1.2, length.out = 50) #in the Liu et al 2009 paper

#Nonparanormal graph

Ymat.npn = huge.npn(data_matrix_std, npn.func="truncation") # Nonparanormal

#modify this length
out.npn = huge(Ymat.npn,lambda = lambda, nlambda = length(lambda),
               method = "glasso", cov.output = TRUE) #using defaults for glasso
#for nlambda and lambda.min.ratio explained in vignette huge paper

#Consider using the screening to speed up the algorithm.

ans = huge.select(out.npn, criterion = "stars")
#Estimated edge matrix for the model 
edge_matrix_nonparanormal_stars <- ans$refit

ans_stars_alternative = huge.select(out.npn, criterion = "stars", stars.thresh = .05) #alternative threshold that was used in the paper.
#Estimated edge matrix for the model 
edge_matrix_nonparanormal_stars_alternative <- ans$refit

ans_ric = huge.select(out.npn, criterion = "ric") 
edge_matrix_nonparanormal_ric <- ans_ric$refit

ans_ebic = huge.select(out.npn, criterion = "ebic")  #said to use ebic for
#glasso in vignette, but ebic you have to tune per huge package.  I am using stars or ric.  RIC is definitely more
#efficient than STARS.
edge_matrix_nonparanormal_ebic <- ans_ebic$refit

#ans$opt.icov #optimal precision matrix.  It does estimate zeros so I 
#create the edge matrix from the precision matrix to get
#the graph

opt_precisionMat = ans$opt.icov
opt_covarianceMat = ans$opt.cov


#Create an adjacency matrix



# ###################################################################
# #Glasso
# 
# out = huge(data_matrix_std,method = "glasso", lambda = lambda, nlambda = length(lambda), cov.output = TRUE)
# 
# 
# ans_glasso = huge.select(out, criterion = "stars") #said to use ebic for
# #glasso in vignette
# 
# 
# opt_precisionMat_glasso = ans_glasso$opt.icov
# opt_covarianceMat_glasso = ans_glasso$opt.cov
# 
# 
# #Create an adjacency matrix
# 
# #Estimated edge matrix for the model 
# edge_matrix_glasso <- ans_glasso$refit

###################################################################

#Save the data for making a graph


writeMat('RealData_Bsplines_paper_nonparanormal_stars_ric_ebic_01.mat', edge_matrix_nonparanormal_stars = 
           edge_matrix_nonparanormal_stars,
         edge_matrix_nonparanormal_ric =  edge_matrix_nonparanormal_ric,  
         edge_matrix_nonparanormal_ebic= edge_matrix_nonparanormal_ebic,edge_matrix_nonparanormal_stars_alternative = 
           edge_matrix_nonparanormal_stars_alternative)


save.image(file = "RealData_Bsplines_paper_newparams_01.rdata")
