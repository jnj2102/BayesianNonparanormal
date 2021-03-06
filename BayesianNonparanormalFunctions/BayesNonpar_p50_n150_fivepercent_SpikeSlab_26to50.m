%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity for spike and slab
%
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=150, p=50, sparsity = fivepercent

load('BayesNonpar_p50_n150_fivepercent_prior.mat');

    c0_list = [0.02; 0.005]; %for the spike scale
    fg_cell = {[1,1]; [10,30]}; %for the Inverse Gamma distribution
    
   [num_elements1, ~] = size(c0_list);
 [num_elements2, ~] = size(fg_cell);

   
	total_time_n150_p50_fivepercent = zeros(num_elements1,num_elements2, reps);
    Sigma_Bayes_est_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
    Omega_Bayes_est_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
    edge_matrix_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
        
      entropy_loss_matrix=  zeros(num_elements1,num_elements2, reps);
    bounded_loss_n150_p50_fivepercent =  zeros(num_elements1,num_elements2, reps);
      Frobenius_norm_precision_n150_p50_fivepercent =  zeros(num_elements1,num_elements2, reps);
   Frobenius_norm_covariance_n150_p50_fivepercent =  zeros(num_elements1,num_elements2, reps);

    SP_matrix=  zeros(num_elements1,num_elements2, reps);
    SE_matrix=  zeros(num_elements1,num_elements2, reps);
    MCC_matrix=  zeros(num_elements1,num_elements2, reps);
     TP_SSVS_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
    TN_SSVS_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
    FP_SSVS_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
    FN_SSVS_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
    mean_Z_Bayes_est_n150_p50_fivepercent = cell([num_elements1,num_elements2, reps]);
        
    %parfor doesn't work with the inner loop indexing.
    
    %Hyperparameter settings for prior probability.
 
    a = 1; %for Beta
    b = 10; %for Beta
    lambda = 1; %for exponential distribution
    
for iters = 26:50
	
    rng(iters,'twister'); %set the seed for each replication for reproducibility
 
   fprintf('Iterations = %d', iters);

    %Now do the B-splines estimation method    
    
    c0_list = [0.02; 0.005]; %for the spike scale
    fg_cell = {[1,1]; [10,30]}; %for the Inverse Gamma distribution
    
    for c_index = 1:num_elements1
        for fg_index = 1:num_elements2
            c0 = c0_list(c_index);
            temp = fg_cell{fg_index};
            f = temp(1);
            g = temp(2);
          
            [entropy_loss,SP_SSVS_total, SE_SSVS_total, MCC_SSVS_total, edge_matrix_ssvs,...
     total_time,TP_SSVS, TN_SSVS, FP_SSVS,FN_SSVS,Omega_Bayes_est,Sigma_Bayes_est,...
    mean_Z_Bayes_est,Frobenius_norm_covariance,bounded_loss,Frobenius_norm_precision] = BayesianNonparanormal_StudentTspikeslab(n,p, sigma_true,...
      x_matrix_n150_p50{iters}, omega_true, F_mat_cell_iters{iters}, g_vec_cell_iters{iters},...
      W_mat_cell_iters{iters}, q_vec_cell_iters{iters},...
	knot_vector_cell_iters{iters}, index_1_cell_iters{iters}, index_2_cell_iters{iters},...
    inverse_variance_prior_reduced_cell_iters{iters}, mean_prior_reduced_cell_iters{iters},...
    Z_red_cell_iters{iters}, Z_two_cell_iters{iters}, c0, a,b,f,g, initial_value_cell_iters{iters}, lambda);



	
%save all of the data for each hyperparameter setting.
    entropy_loss_matrix(c_index, fg_index,iters) = entropy_loss;
    bounded_loss_n150_p50_fivepercent(c_index, fg_index,iters) = bounded_loss; 
   Frobenius_norm_precision_n150_p50_fivepercent(c_index, fg_index,iters) = Frobenius_norm_precision; 
   Frobenius_norm_covariance_n150_p50_fivepercent(c_index, fg_index,iters) = Frobenius_norm_covariance; 

    SP_matrix(c_index, fg_index,iters) = SP_SSVS_total;
    SE_matrix(c_index, fg_index,iters) = SE_SSVS_total;
    MCC_matrix(c_index, fg_index,iters) = MCC_SSVS_total;
   
   total_time_n150_p50_fivepercent(c_index, fg_index,iters) = total_time; %this is the
        %total time it took to run the sampling using q_iters

    edge_matrix_n150_p50_fivepercent{c_index, fg_index,iters} = edge_matrix_ssvs;

    Sigma_Bayes_est_n150_p50_fivepercent{c_index, fg_index,iters} = Sigma_Bayes_est;
    Omega_Bayes_est_n150_p50_fivepercent{c_index, fg_index,iters} = Omega_Bayes_est;
    mean_Z_Bayes_est_n150_p50_fivepercent{c_index, fg_index,iters} = mean_Z_Bayes_est;
    TP_SSVS_n150_p50_fivepercent{c_index, fg_index,iters} = TP_SSVS;
    TN_SSVS_n150_p50_fivepercent{c_index, fg_index,iters} = TN_SSVS;
    FP_SSVS_n150_p50_fivepercent{c_index, fg_index,iters} = FP_SSVS;
    FN_SSVS_n150_p50_fivepercent{c_index, fg_index,iters} = FN_SSVS;
   
        end
   end
        
end


save('BayesNonpar_p50_n150_fivepercent_SpikeSlab_26to50.mat', '-v7.3');


    