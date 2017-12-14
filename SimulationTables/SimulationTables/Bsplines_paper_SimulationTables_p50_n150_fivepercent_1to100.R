########################################################
#Code to put the data from simulation study into a table
#
#Author: Jami Jackson Mulgrave
#
########################################################

#clear the workspace
rm(list = ls())

library(R.matlab)
library(xtable)


## Frequentist p=50 n=150 sparsity=fivepercent

load("Bsplines_paper_frequentist_p50_n150_fivepercent.rdata")

#####EBIC

#find SP and SE
means_SP_p50_n150_fivepercent_ebic <- c(
  round(mean_SP_freq_n150_p50_fivepercent_ebic,2))


 stderr_SP_p50_n150_fivepercent_ebic <- c(
   round(stderr_SP_freq_n150_p50_fivepercent_ebic,2))
 
 stderr_SP_p50_n150_fivepercent_ebic <- paste("(", 
                                formatC(as.vector(stderr_SP_p50_n150_fivepercent_ebic), format = "f", digits =2),
                                ")", sep="")
 
 means_SE_p50_n150_fivepercent_ebic <- c(
   round(mean_SE_freq_n150_p50_fivepercent_ebic,2))
 
 
 stderr_SE_p50_n150_fivepercent_ebic <- c(
   round(stderr_SE_freq_n150_p50_fivepercent_ebic,2))
 
 stderr_SE_p50_n150_fivepercent_ebic <- paste("(", 
                                            formatC(as.vector(stderr_SE_p50_n150_fivepercent_ebic), format = "f", digits =2),
                                            ")", sep="")

 
 #Find L1
 
 means_L1_p50_n150_fivepercent_ebic <- c(
   round(mean_L1_freq_n150_p50_fivepercent_ebic,2))
 
 
 stderr_L1_p50_n150_fivepercent_ebic <- c(
   round(stderr_L1_freq_n150_p50_fivepercent_ebic,2))
 
 
 stderr_L1_p50_n150_fivepercent_ebic <- paste("(", 
                                         formatC(as.vector(stderr_L1_p50_n150_fivepercent_ebic), format = "f", digits =2),
                                         ")", sep="")
 
 #Find MCC
 means_MCC_p50_n150_fivepercent_ebic <- c(
   round(mean_MCC_freq_n150_p50_fivepercent_ebic,2))
 
 
 stderr_MCC_p50_n150_fivepercent_ebic <- c(
   round(stderr_MCC_freq_n150_p50_fivepercent_ebic,2))
 
 stderr_MCC_p50_n150_fivepercent_ebic <- paste("(", 
                                         formatC(as.vector(stderr_MCC_p50_n150_fivepercent_ebic), format = "f", digits =2),
                                         ")", sep="")
 
 
 
 #####STARS with threshold = .1
 
 #find SP and SE

 means_SP_p50_n150_fivepercent_stars <- c(
   round(mean_SP_freq_n150_p50_fivepercent_stars,2))
 
 
 stderr_SP_p50_n150_fivepercent_stars <- c(
   round(stderr_SP_freq_n150_p50_fivepercent_stars,2))
 
 stderr_SP_p50_n150_fivepercent_stars <- paste("(", 
                                            formatC(as.vector(stderr_SP_p50_n150_fivepercent_stars), format = "f", digits =2),
                                            ")", sep="")
 
 means_SE_p50_n150_fivepercent_stars <- c(
   round(mean_SE_freq_n150_p50_fivepercent_stars,2))
 
 
 stderr_SE_p50_n150_fivepercent_stars <- c(
   round(stderr_SE_freq_n150_p50_fivepercent_stars,2))
 
 stderr_SE_p50_n150_fivepercent_stars <- paste("(", 
                                            formatC(as.vector(stderr_SE_p50_n150_fivepercent_stars), format = "f", digits =2),
                                            ")", sep="")
 
 #Find L1
 
 means_L1_p50_n150_fivepercent_stars <- c(
   round(mean_L1_freq_n150_p50_fivepercent_stars,2))
 
 
 stderr_L1_p50_n150_fivepercent_stars <- c(
   round(stderr_L1_freq_n150_p50_fivepercent_stars,2))
 
 
 stderr_L1_p50_n150_fivepercent_stars <- paste("(", 
                                             formatC(as.vector(stderr_L1_p50_n150_fivepercent_stars), format = "f", digits =2),
                                             ")", sep="")
 
 #Find MCC
 means_MCC_p50_n150_fivepercent_stars <- c(
   round(mean_MCC_freq_n150_p50_fivepercent_stars,2))
 
 
 stderr_MCC_p50_n150_fivepercent_stars <- c(
   round(stderr_MCC_freq_n150_p50_fivepercent_stars,2))
 
 stderr_MCC_p50_n150_fivepercent_stars <- paste("(", 
                                              formatC(as.vector(stderr_MCC_p50_n150_fivepercent_stars), format = "f", digits =2),
                                              ")", sep="")
 
 
 
 
 #####STARS with threshold = .05
 
 #find SP and SE
 
 means_SP_p50_n150_fivepercent_stars_alternative <- c(
   round(mean_SP_freq_n150_p50_fivepercent_stars_alternative,2))
 
 
 stderr_SP_p50_n150_fivepercent_stars_alternative <- c(
   round(stderr_SP_freq_n150_p50_fivepercent_stars_alternative,2))
 
 stderr_SP_p50_n150_fivepercent_stars_alternative <- paste("(", 
                                                         formatC(as.vector(stderr_SP_p50_n150_fivepercent_stars_alternative), format = "f", digits =2),
                                                         ")", sep="")
 
 means_SE_p50_n150_fivepercent_stars_alternative <- c(
   round(mean_SE_freq_n150_p50_fivepercent_stars_alternative,2))
 
 
 stderr_SE_p50_n150_fivepercent_stars_alternative <- c(
   round(stderr_SE_freq_n150_p50_fivepercent_stars_alternative,2))
 
 stderr_SE_p50_n150_fivepercent_stars_alternative <- paste("(", 
                                                         formatC(as.vector(stderr_SE_p50_n150_fivepercent_stars_alternative), format = "f", digits =2),
                                                         ")", sep="")
 
 #Find L1
 
 means_L1_p50_n150_fivepercent_stars_alternative <- c(
   round(mean_L1_freq_n150_p50_fivepercent_stars_alternative,2))
 
 
 stderr_L1_p50_n150_fivepercent_stars_alternative <- c(
   round(stderr_L1_freq_n150_p50_fivepercent_stars_alternative,2))
 
 
 stderr_L1_p50_n150_fivepercent_stars_alternative <- paste("(", 
                                                         formatC(as.vector(stderr_L1_p50_n150_fivepercent_stars_alternative), format = "f", digits =2),
                                                         ")", sep="")
 
 #Find MCC
 means_MCC_p50_n150_fivepercent_stars_alternative <- c(
   round(mean_MCC_freq_n150_p50_fivepercent_stars_alternative,2))
 
 
 stderr_MCC_p50_n150_fivepercent_stars_alternative <- c(
   round(stderr_MCC_freq_n150_p50_fivepercent_stars_alternative,2))
 
 stderr_MCC_p50_n150_fivepercent_stars_alternative <- paste("(", 
                                                          formatC(as.vector(stderr_MCC_p50_n150_fivepercent_stars_alternative), format = "f", digits =2),
                                                          ")", sep="")
 
 
 
 #####RIC
 
 #find SP and SE
 
 means_SP_p50_n150_fivepercent_ric <- c(
   round(mean_SP_freq_n150_p50_fivepercent_ric,2))
 
 
 stderr_SP_p50_n150_fivepercent_ric <- c(
   round(stderr_SP_freq_n150_p50_fivepercent_ric,2))
 
 stderr_SP_p50_n150_fivepercent_ric <- paste("(", 
                                           formatC(as.vector(stderr_SP_p50_n150_fivepercent_ric), format = "f", digits =2),
                                           ")", sep="")
 
 means_SE_p50_n150_fivepercent_ric <- c(
   round(mean_SE_freq_n150_p50_fivepercent_ric,2))
 
 
 stderr_SE_p50_n150_fivepercent_ric <- c(
   round(stderr_SE_freq_n150_p50_fivepercent_ric,2))
 
 stderr_SE_p50_n150_fivepercent_ric <- paste("(", 
                                           formatC(as.vector(stderr_SE_p50_n150_fivepercent_ric), format = "f", digits =2),
                                           ")", sep="")
 
 
 #Find L1
 
 means_L1_p50_n150_fivepercent_ric <- c(
   round(mean_L1_freq_n150_p50_fivepercent_ric,2))
 
 
 stderr_L1_p50_n150_fivepercent_ric <- c(
   round(stderr_L1_freq_n150_p50_fivepercent_ric,2))
 
 
 stderr_L1_p50_n150_fivepercent_ric <- paste("(", 
                                           formatC(as.vector(stderr_L1_p50_n150_fivepercent_ric), format = "f", digits =2),
                                           ")", sep="")
 
 #Find MCC
 means_MCC_p50_n150_fivepercent_ric <- c(
   round(mean_MCC_freq_n150_p50_fivepercent_ric,2))
 
 
 stderr_MCC_p50_n150_fivepercent_ric <- c(
   round(stderr_MCC_freq_n150_p50_fivepercent_ric,2))
 
 stderr_MCC_p50_n150_fivepercent_ric <- paste("(", 
                                            formatC(as.vector(stderr_MCC_p50_n150_fivepercent_ric), format = "f", digits =2),
                                            ")", sep="") 

 
 ###Bayesian nonparanormal ###
 
 
 #Read in the Matlab data and put it in a latex table for Bayes method
 
 
 data_summary <- readMat('Bsplines_paper_Calcs_p50_n150_fivepercent_SpikeSlab_1to100.mat',
                         package = "R.matlab")
 #L1
 Bayes_L1_mean <- data_summary$L1.n150.p50.fivepercent.mean.rounded
 
 Bayes_L1_stderr <- data_summary$L1.n150.p50.fivepercent.SE.rounded
 
 Bayes_L1_stderr <- paste("(", formatC(as.vector(Bayes_L1_stderr), format = "f", digits =2),
                          ")", sep="")
 
 
 #SP
 Bayes_SP_mean <- data_summary$SP.n150.p50.fivepercent.mean.rounded
 
 Bayes_SP_stderr <- data_summary$SP.n150.p50.fivepercent.SE.rounded
 
 Bayes_SP_stderr <- paste("(", formatC(as.vector(Bayes_SP_stderr), format = "f", digits =2),
                          ")", sep="")
 
 #SE
 Bayes_SE_mean <- data_summary$SE.n150.p50.fivepercent.mean.rounded
 
 Bayes_SE_stderr <- data_summary$SE.n150.p50.fivepercent.SE.rounded
 
 Bayes_SE_stderr <- paste("(", formatC(as.vector(Bayes_SE_stderr), format = "f", digits =2),
                          ")", sep="")
 
 #MCC
 Bayes_MCC_mean <- data_summary$MCC.n150.p50.fivepercent.mean.rounded
 
 Bayes_MCC_stderr <- data_summary$MCC.n150.p50.fivepercent.SE.rounded
 
 Bayes_MCC_stderr <- paste("(", formatC(as.vector(Bayes_MCC_stderr), format = "f", digits =2),
                           ")", sep="")
 
###

#L1 Matrix
summary_mat_L1 <- cbind(
                     rbind(formatC(means_L1_p50_n150_fivepercent_ebic, format = "f", digits=2),
                     stderr_L1_p50_n150_fivepercent_ebic),
                     
                     rbind(formatC(means_L1_p50_n150_fivepercent_stars, format = "f", digits=2),
                     stderr_L1_p50_n150_fivepercent_stars),
                     
                     rbind(formatC(means_L1_p50_n150_fivepercent_stars_alternative, format = "f", digits=2),
                     stderr_L1_p50_n150_fivepercent_stars_alternative),
                     
                     rbind(formatC(means_L1_p50_n150_fivepercent_ric, format = "f", digits=2),
                     stderr_L1_p50_n150_fivepercent_ric),
                     
                     rbind(formatC(Bayes_L1_mean, format = "f", digits=2),Bayes_L1_stderr ))


rownames(summary_mat_L1) <- c("fivepercentMean", "fivepercentSE")




colnames(summary_mat_L1) <- c("EBICL1", "STARSL1", "STARALTL1", "RICL1", "Bayes")


print(xtable(summary_mat_L1),type="latex")   


#SP Matrix
summary_mat_SP <- cbind(
  rbind(formatC(means_SP_p50_n150_fivepercent_ebic, format = "f", digits=2),
        stderr_SP_p50_n150_fivepercent_ebic),
  
  rbind(formatC(means_SP_p50_n150_fivepercent_stars, format = "f", digits=2),
        stderr_SP_p50_n150_fivepercent_stars),
  
  rbind(formatC(means_SP_p50_n150_fivepercent_stars_alternative, format = "f", digits=2),
        stderr_SP_p50_n150_fivepercent_stars_alternative),
  
  rbind(formatC(means_SP_p50_n150_fivepercent_ric, format = "f", digits=2),
        stderr_SP_p50_n150_fivepercent_ric),
  
  rbind(formatC(Bayes_SP_mean, format = "f", digits=2), Bayes_SP_stderr))


rownames(summary_mat_SP) <- c("fivepercentMean", "fivepercentSE")




colnames(summary_mat_SP) <- c("EBICSP", "STARSSP", "STARALTSP", "RICSP", "Bayes")


print(xtable(summary_mat_SP),type="latex")   


#SE Matrix
summary_mat_SE <- cbind(
  rbind(formatC(means_SE_p50_n150_fivepercent_ebic, format = "f", digits=2),
        stderr_SE_p50_n150_fivepercent_ebic),
  
  rbind(formatC(means_SE_p50_n150_fivepercent_stars, format = "f", digits=2),
        stderr_SE_p50_n150_fivepercent_stars),
  
  rbind(formatC(means_SE_p50_n150_fivepercent_stars_alternative, format = "f", digits=2),
        stderr_SE_p50_n150_fivepercent_stars_alternative),
  
  rbind(formatC(means_SE_p50_n150_fivepercent_ric, format = "f", digits=2),
        stderr_SE_p50_n150_fivepercent_ric),
  
  rbind(formatC(Bayes_SE_mean, format = "f", digits=2), Bayes_SE_stderr))



rownames(summary_mat_SE) <- c("fivepercentMean", "fivepercentSE")




colnames(summary_mat_SE) <- c("EBICSE", "STARSSE", "STARALTSE", "RICSE", "Bayes")


print(xtable(summary_mat_SE),type="latex")   


#MCC Matrix
summary_mat_MCC <- cbind(
  rbind(formatC(means_MCC_p50_n150_fivepercent_ebic, format = "f", digits=2),
        stderr_MCC_p50_n150_fivepercent_ebic),
  
  rbind(formatC(means_MCC_p50_n150_fivepercent_stars, format = "f", digits=2),
        stderr_MCC_p50_n150_fivepercent_stars),
  
  rbind(formatC(means_MCC_p50_n150_fivepercent_stars_alternative, format = "f", digits=2),
        stderr_MCC_p50_n150_fivepercent_stars_alternative),
  
  rbind(formatC(means_MCC_p50_n150_fivepercent_ric, format = "f", digits=2),
        stderr_MCC_p50_n150_fivepercent_ric),
  
  rbind(formatC(Bayes_MCC_mean, format = "f", digits=2), Bayes_MCC_stderr))



rownames(summary_mat_MCC) <- c("fivepercentMean", "fivepercentMCC")




colnames(summary_mat_MCC) <- c("EBICMCC", "STARSMCC", "STARALTMCC", "RICMCC", "Bayes")


print(xtable(summary_mat_MCC),type="latex")   



#Time it takes to run

# A function to remove the NAs to find the standard error
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

#Make the algorithm_time into a vector

algorithm_time_vector <- c()

for (iters in 1:reps) {
  algorithm_time_vector    <-  c(algorithm_time_vector, algorithm_time[[iters]][1])
}


mean_algorithm_time_freq_n150_p50_fivepercent = round(mean(algorithm_time_vector),2)

stderr_algorithm_time_freq_n150_p50_fivepercent = stderr(algorithm_time_vector)


stderr_algorithm_time_freq_n150_p50_fivepercent <- paste("(", 
                                                   formatC(as.vector(stderr_algorithm_time_freq_n150_p50_fivepercent), format = "f", digits =2),
                                                   ")", sep="")

#Time it takes to run the Bayesian nonparanormal method

Bayes_total_time_mean <- data_summary$total.time.n150.p50.fivepercent.mean.rounded

Bayes_total_time_stderr <- data_summary$total.time.n150.p50.fivepercent.SE.rounded

Bayes_total_time_stderr <- paste("(", formatC(as.vector(Bayes_total_time_stderr), format = "f", digits =2),
                                 ")", sep="")

