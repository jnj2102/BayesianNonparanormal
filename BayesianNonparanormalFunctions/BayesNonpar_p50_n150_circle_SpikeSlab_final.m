%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=500, p=100, sparsity = circle


 ssIters = 1;
 ssIters2 = ssIters + 24;
 
 while ssIters <= 100
     
     
%Horseshoe
load(sprintf('BayesNonpar_p50_n150_circle_SpikeSlab_%dto%d.mat', ssIters,ssIters2) );

for num_iters = ssIters:ssIters2
    rng(num_iters,'twister'); %set the seed for each replication for reproducibility

 
   fprintf('Iterations = %d', num_iters);

    %Now do the B-splines estimation method    
    
    for c_index = 1:num_elements1
        for fg_index = 1:num_elements2
               
        Omega_Bayes_est =  Omega_Bayes_est_n150_p50_circle{c_index, fg_index, num_iters};
        final_edge_matrix = edge_matrix_n150_p50_circle{c_index, fg_index, num_iters};
        
       mean_Z_Bayes_est = mean_Z_Bayes_est_n150_p50_circle{c_index, fg_index, num_iters};
       
       %%%Find the MLE using a convex optimization problem for eBIC

     %%% MLE with Bglasso
     
[BIC]= FindBIC(Omega_Bayes_est,final_edge_matrix,mean_Z_Bayes_est,n,p);

     BIC_matrix(c_index, fg_index, num_iters) = BIC;
   
        end
    end

end


    
    %Find the minimum BIC and keep the data that pertains to it for further
    %analysis 
for num_iters = ssIters:ssIters2
[minBIC, ~ ] = min(min(BIC_matrix(:,:,num_iters)));

[rowMinBIC, colMinBIC] = find(BIC_matrix(:,:,num_iters) == minBIC);

%save the data with the minimum BIC for further analysis
 
 SP_matrix_finalanalysis(num_iters)  =  SP_matrix(rowMinBIC, colMinBIC ,  num_iters);
 SE_matrix_finalanalysis(num_iters) =  SE_matrix(rowMinBIC, colMinBIC ,  num_iters);
 MCC_matrix_finalanalysis(num_iters) =  MCC_matrix(rowMinBIC, colMinBIC ,  num_iters);
 
BIC_matrix_finalanalysis(num_iters) =   BIC_matrix(rowMinBIC, colMinBIC ,  num_iters);

 total_time_finalanalysis(num_iters) =   total_time_n150_p50_circle(rowMinBIC, colMinBIC ,  num_iters);
entropy_loss_finalanalysis(num_iters)  = entropy_loss_matrix(rowMinBIC, colMinBIC, num_iters) ;
  bounded_loss_finalanalysis(num_iters)  = bounded_loss_n150_p50_circle(rowMinBIC, colMinBIC, num_iters) ;
   Frobenius_norm_precision_finalanalysis(num_iters)  = Frobenius_norm_precision_n150_p50_circle(rowMinBIC, colMinBIC, num_iters) ;
  Frobenius_norm_covariance_finalanalysis(num_iters)  = Frobenius_norm_covariance_n150_p50_circle(rowMinBIC, colMinBIC, num_iters) ;

 Omega_Bayes_est_finalanalysis{num_iters} =   Omega_Bayes_est_n150_p50_circle{rowMinBIC, colMinBIC ,  num_iters};
Sigma_Bayes_est_finalanalysis{num_iters} = Sigma_Bayes_est_n150_p50_circle{rowMinBIC, colMinBIC ,  num_iters};
edge_matrix_finalanalysis{num_iters} =  edge_matrix_n150_p50_circle{rowMinBIC, colMinBIC ,  num_iters};
   TP_finalanalysis(num_iters) =  TP_SSVS_n150_p50_circle(rowMinBIC, colMinBIC ,  num_iters);
  TN_finalanalysis(num_iters) =  TN_SSVS_n150_p50_circle(rowMinBIC, colMinBIC ,  num_iters);
  FP_finalanalysis(num_iters) =  FP_SSVS_n150_p50_circle(rowMinBIC, colMinBIC ,  num_iters);
   FN_finalanalysis(num_iters) = FN_SSVS_n150_p50_circle(rowMinBIC, colMinBIC ,  num_iters);
   
  mean_Z_Bayes_est_finalanalysis{num_iters} = mean_Z_Bayes_est_n150_p50_circle{rowMinBIC, colMinBIC ,  num_iters};
  
end

ssIters = ssIters2 +1;
ssIters2 = ssIters + 24;

 end
 
 
save('BayesNonpar_p50_n150_circle_SpikeSlab_final.mat', '-v7.3');
