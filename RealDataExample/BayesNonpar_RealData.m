%%Real data application for B-splines paper
%Author: Jami Jackson Mulgrave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

load('Bsplines_paper_realdata_initialdata.mat');

%log transform the data and standardize the each gene

data_matrix_log = log(data_matrix)';

data_matrix_std = (data_matrix_log - min(data_matrix_log))./(max(data_matrix_log) - min(data_matrix_log));

%Now try my method out - data_matrix_std is my "x_matrix"

[n,p] = size(data_matrix_std);

%these (mu, tau, sigma2) are for my prior for the transformation functions
mu = 1; %I'm just giving it a mean of 1 for now but it can be any constant
tau = 0.5; %I'm giving it a sd of 0.5 for now but it can be any constant

sigma2 = 1; %I'm making the variance 1 for now but it can be any constant

N_points = 20;

c = [0;1]; %vector of linear constraints

  
	
 rng(5000,'twister'); %set the seed for each replication for reproducibility


[optimalJ, minK, Final_AIC,F_mat_cell,g_vec_cell,...
Fconstraint_cell,RHS_gq_cell,A_mat_cell,LHS_gq_cell,W_mat_cell,...
q_vec_cell,knot_vector_cell,index_1_cell,index_2_cell,...
inverse_variance_prior_reduced_cell,mean_prior_reduced_cell,...
Z_red_cell,Z_two_cell,initial_value] = PriorBsplines_Initialvalues_corrected(n,p, data_matrix_std, mu, tau, sigma2,c,...
N_points);

    
  %now run the method
  
     
    %Hyperparameter settings for prior probability.
 
    a = 1; %for Beta
    b = 10; %for Beta
    lambda = 1; %for exponential distribution
	
    
	c0_list = [0.02; 0.005]; %for the spike scale
    fg_cell = {[1,1]; [10,30]}; %for the Inverse Gamma distribution
        
   [num_elements1, ~] = size(c0_list);
 [num_elements2, ~] = size(fg_cell);
 
 iters = 1;
 
 sigma_true = eye(p);
 omega_true = eye(p);
    
for c_index = 1:num_elements1
       for fg_index = 1:num_elements2
            c0 = c0_list(c_index);
            temp = fg_cell{fg_index};
            f = temp(1);
            g = temp(2);
          
             rng(2000,'twister'); %set the seed for each replication for reproducibility

            [~,~, ~, ~, edge_matrix_ssvs,...
     total_time,~, ~, ~,~,Omega_Bayes_est,Sigma_Bayes_est,...
    mean_Z_Bayes_est,~,~,~] = BayesianNonparanormal_StudentTspikeslab(n,p, sigma_true,...
      data_matrix_std, omega_true, F_mat_cell, g_vec_cell,...
      W_mat_cell, q_vec_cell,...
	knot_vector_cell, index_1_cell, index_2_cell,...
    inverse_variance_prior_reduced_cell, mean_prior_reduced_cell,...
    Z_red_cell, Z_two_cell, c0, a,b,f,g, initial_value, lambda);

%save all of the data for each hyperparameter setting.
     
   total_time_matrix(c_index, fg_index,iters) = total_time; %this is the
        %total time it took to run the sampling using q_iters

    edge_matrix_matrix{c_index, fg_index,iters} = edge_matrix_ssvs;

    Sigma_Bayes_est_matrix{c_index, fg_index,iters} = Sigma_Bayes_est;
    Omega_Bayes_est_matrix{c_index, fg_index,iters} = Omega_Bayes_est;
    mean_Z_Bayes_est_matrix{c_index, fg_index,iters} = mean_Z_Bayes_est;
  
       end
	
end
	
    %Find the MLE using a convex optimization problem
    for c_index = 1:num_elements1
        for fg_index = 1:num_elements2
               
        Omega_Bayes_est =  Omega_Bayes_est_matrix{c_index, fg_index, iters};
        final_edge_matrix = edge_matrix_matrix{c_index, fg_index, iters};
        
       mean_Z_Bayes_est = mean_Z_Bayes_est_matrix{c_index, fg_index, iters};
       
       %%%Find the MLE using a convex optimization problem for eBIC

     %%% MLE with Bglasso
     
[BIC]= FindBIC(Omega_Bayes_est,final_edge_matrix,mean_Z_Bayes_est,n,p);

     BIC_matrix(c_index, fg_index, iters) = BIC;
   
        end
    end
	
	
	    %Find the minimum BIC and keep the data that pertains to it for further
    %analysis 


  [minBIC, ~ ] = min(min(BIC_matrix(:,:,iters)));

[rowMinBIC, colMinBIC] = find(BIC_matrix(:,:,iters) == minBIC);

%save the data with the minimum BIC for further analysis
 
BIC_matrix_finalanalysis(iters) =   BIC_matrix(rowMinBIC, colMinBIC,  iters);

 total_time_finalanalysis(iters) =   total_time_matrix(rowMinBIC, colMinBIC,  iters);
 Omega_Bayes_est_finalanalysis{iters} =   Omega_Bayes_est_matrix{rowMinBIC, colMinBIC,  iters};
Sigma_Bayes_est_finalanalysis{iters} = Sigma_Bayes_est_matrix{rowMinBIC, colMinBIC,  iters};
edge_matrix_finalanalysis{iters} =  edge_matrix_matrix{rowMinBIC, colMinBIC,  iters};
   
  mean_Z_Bayes_est_finalanalysis{iters} = mean_Z_Bayes_est_matrix{rowMinBIC, colMinBIC,  iters};
  
  save('BayesNonpar_RealData.mat', '-v7.3');
