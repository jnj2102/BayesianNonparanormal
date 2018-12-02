%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=50, p=25, sparsity = tenpercent

%Horseshoe
load('BayesNonpar_p25_n50_tenpercent_SpikeSlab.mat');


for iters = 1:reps
	
    rng(iters,'twister'); %set the seed for each replication for reproducibility

 
   fprintf('Iterations = %d', iters);

    %Now do the B-splines estimation method    
    
    for c_index = 1:num_elements1
        for fg_index = 1:num_elements2
               
        Omega_Bayes_est =  Omega_Bayes_est_n50_p25_tenpercent{c_index, fg_index, iters};
        final_edge_matrix = edge_matrix_n50_p25_tenpercent{c_index, fg_index, iters};
        
       mean_Z_Bayes_est = mean_Z_Bayes_est_n50_p25_tenpercent{c_index, fg_index, iters};
       
       %%%Find the MLE using a convex optimization problem for eBIC

     %%% MLE with Bglasso
     
[BIC]= FindBIC(Omega_Bayes_est,final_edge_matrix,mean_Z_Bayes_est,n,p);

     BIC_matrix(c_index, fg_index, iters) = BIC;
   
        end
    end

end


    
    %Find the minimum BIC and keep the data that pertains to it for further
    %analysis 
for iters = 1:reps
[minBIC, ~ ] = min(min(BIC_matrix(:,:,iters)));

[rowMinBIC, colMinBIC] = find(BIC_matrix(:,:,iters) == minBIC);

%save the data with the minimum BIC for further analysis
 
 SP_matrix_finalanalysis(iters)  =  SP_matrix(rowMinBIC, colMinBIC,  iters);
 SE_matrix_finalanalysis(iters) =  SE_matrix(rowMinBIC, colMinBIC,  iters);
 MCC_matrix_finalanalysis(iters) =  MCC_matrix(rowMinBIC, colMinBIC,  iters);
 
BIC_matrix_finalanalysis(iters) =   BIC_matrix(rowMinBIC, colMinBIC,  iters);

 total_time_finalanalysis(iters) =   total_time_n50_p25_tenpercent(rowMinBIC, colMinBIC,  iters);
entropy_loss_finalanalysis(iters)  = entropy_loss_matrix(rowMinBIC, colMinBIC, iters) ;
  bounded_loss_finalanalysis(iters)  = bounded_loss_n50_p25_tenpercent(rowMinBIC, colMinBIC, iters) ;
   Frobenius_norm_precision_finalanalysis(iters)  = Frobenius_norm_precision_n50_p25_tenpercent(rowMinBIC, colMinBIC, iters) ;
  Frobenius_norm_covariance_finalanalysis(iters)  = Frobenius_norm_covariance_n50_p25_tenpercent(rowMinBIC, colMinBIC, iters) ;

 Omega_Bayes_est_finalanalysis{iters} =   Omega_Bayes_est_n50_p25_tenpercent{rowMinBIC, colMinBIC,  iters};
Sigma_Bayes_est_finalanalysis{iters} = Sigma_Bayes_est_n50_p25_tenpercent{rowMinBIC, colMinBIC,  iters};
edge_matrix_finalanalysis{iters} =  edge_matrix_n50_p25_tenpercent{rowMinBIC, colMinBIC,  iters};
   TP_finalanalysis(iters) =  TP_SSVS_n50_p25_tenpercent(rowMinBIC, colMinBIC,  iters);
  TN_finalanalysis(iters) =  TN_SSVS_n50_p25_tenpercent(rowMinBIC, colMinBIC,  iters);
  FP_finalanalysis(iters) =  FP_SSVS_n50_p25_tenpercent(rowMinBIC, colMinBIC,  iters);
   FN_finalanalysis(iters) = FN_SSVS_n50_p25_tenpercent(rowMinBIC, colMinBIC,  iters);
   
  mean_Z_Bayes_est_finalanalysis{iters} = mean_Z_Bayes_est_n50_p25_tenpercent{rowMinBIC, colMinBIC,  iters};
  
end



save('BayesNonpar_p25_n50_tenpercent_SpikeSlab_final.mat', '-v7.3');
