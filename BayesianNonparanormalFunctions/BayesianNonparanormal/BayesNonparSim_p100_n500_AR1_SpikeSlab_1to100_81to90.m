%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity for spike and slab
%For first paper on Bsplines
%
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=500, p=100, sparsity = AR1

    load('Bsplines_paper_Sim_n500_p100_AR1_prior_1to100.mat');

    c0_list = [0.02; 0.00025];
    ab_cell = {[5,25]; [10,30]}; 
    
   [num_elements1, ~] = size(c0_list);
 [num_elements2, ~] = size(ab_cell);

   
	total_time_n500_p100_AR1 = zeros(num_elements1,num_elements2, reps);
    L1_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    SP_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    SE_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    MCC_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    %f_sup_n500_p100_AR1  = cell([num_elements1,num_elements2, reps]);
    q_iters_n500_p100_AR1 = zeros(num_elements1,num_elements2, reps);
    mean_est_n500_p100_AR1  = cell([num_elements1,num_elements2, reps]);
    final_thetas_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    final_mean_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    burnin_n500_p100_AR1 = zeros(num_elements1,num_elements2, reps);
    y_transformed_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    All_tmp_iters_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    complete_Sig_total_sims_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    complete_edge_matrix_total_sims_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    complete_Omega_total_sims_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    edge_matrix_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
        
    BIC_matrix=  zeros(num_elements1,num_elements2, reps);
    eBIC_matrix=  zeros(num_elements1,num_elements2, reps);
    eBIC_half_matrix=  zeros(num_elements1,num_elements2, reps);
    L1_matrix=  zeros(num_elements1,num_elements2, reps);
    SP_matrix=  zeros(num_elements1,num_elements2, reps);
    SE_matrix=  zeros(num_elements1,num_elements2, reps);
    MCC_matrix=  zeros(num_elements1,num_elements2, reps);
     TP_SSVS_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    TN_SSVS_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    FP_SSVS_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    FN_SSVS_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    allOmega_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    allSig_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    allWeights_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    monitor_omega_chain_firsthalf_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    monitor_omega_chain_secondhalf_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    allW_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    Omega_L1_SSVS_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
    complete_Z_Bayes_est_sims_n500_p100_AR1 = cell([num_elements1,num_elements2, reps]);
        
    %parfor doesn't work with the inner loop indexing.
    
    %Hyperparameter settings for prior probability.
    f = 1;
    g = 1;

for iters = 81:90
	
    rng(iters,'twister'); %set the seed for each replication for reproducibility
 
   fprintf('Iterations = %d', iters);

    %Now do the B-splines estimation method    
    
    c0_list = [0.02; 0.00025];
    ab_cell = {[5,25]; [10,30]}; 
    
    for c_index = 1:num_elements1
        for ab_index = 1:num_elements2
            c0 = c0_list(c_index);
            temp = ab_cell{ab_index};
            a = temp(1);
            b = temp(2);
          
[L1_SSVS_save,SP_SSVS_total, SE_SSVS_total, MCC_SSVS_total,edge_matrix_ssvs, q_iters,...
    complete_thetas_total_sims, final_burnin , Y_transformed_matrix, All_tmp_iters, total_time_q_iters, ...
    complete_Sig_total_sims, complete_edge_matrix_total_sims, complete_Omega_total_sims, ...
    complete_mean_sims,TP_SSVS, TN_SSVS, FP_SSVS,FN_SSVS,allOmega, allW, allSig, allWeights,...
    monitor_omega_chain_firsthalf, monitor_omega_chain_secondhalf,Omega_L1_SSVS,complete_Z_Bayes_est_sims] = BayesianNonparanormal_joint_clean(n,p,sigma_true,...
     x_matrix_n500_p100{iters}, omega_true, F_mat_cell_iters{iters}, g_vec_cell_iters{iters},...
     W_mat_cell_iters{iters}, q_vec_cell_iters{iters},...
	knot_vector_cell_iters{iters}, index_1_cell_iters{iters}, index_2_cell_iters{iters}, inverse_variance_prior_reduced_cell_iters{iters},...
    mean_prior_reduced_cell_iters{iters},...
    Z_red_cell_iters{iters}, Z_two_cell_iters{iters},...
    c0, a,b,f,g, initial_value_cell_iters{iters});


	
%save all of the data for each hyperparameter setting.
    L1_matrix(c_index, ab_index,iters) = L1_SSVS_save;
    SP_matrix(c_index, ab_index,iters) = SP_SSVS_total;
    SE_matrix(c_index, ab_index,iters) = SE_SSVS_total;
    MCC_matrix(c_index, ab_index,iters) = MCC_SSVS_total;
   
   total_time_n500_p100_AR1(c_index, ab_index,iters) = total_time_q_iters; %this is the
        %total time it took to run the sampling using q_iters

    edge_matrix_n500_p100_AR1{c_index, ab_index,iters} = edge_matrix_ssvs;

    q_iters_n500_p100_AR1(c_index, ab_index,iters) = q_iters;
    final_thetas_n500_p100_AR1{c_index, ab_index,iters} = complete_thetas_total_sims;
    final_mean_n500_p100_AR1{c_index, ab_index,iters} = complete_mean_sims;
    burnin_n500_p100_AR1(c_index, ab_index,iters) = final_burnin;
    y_transformed_n500_p100_AR1{c_index, ab_index,iters} = Y_transformed_matrix;
    All_tmp_iters_n500_p100_AR1{c_index, ab_index,iters} = All_tmp_iters;
    complete_Sig_total_sims_n500_p100_AR1{c_index, ab_index,iters} = complete_Sig_total_sims;
    complete_edge_matrix_total_sims_n500_p100_AR1{c_index, ab_index,iters} = complete_edge_matrix_total_sims;
    complete_Omega_total_sims_n500_p100_AR1{c_index, ab_index,iters} = complete_Omega_total_sims;
    Omega_L1_SSVS_n500_p100_AR1{c_index, ab_index,iters} = Omega_L1_SSVS;
    complete_Z_Bayes_est_sims_n500_p100_AR1{c_index, ab_index,iters} = complete_Z_Bayes_est_sims;
    TP_SSVS_n500_p100_AR1{c_index, ab_index,iters} = TP_SSVS;
    TN_SSVS_n500_p100_AR1{c_index, ab_index,iters} = TN_SSVS;
    FP_SSVS_n500_p100_AR1{c_index, ab_index,iters} = FP_SSVS;
    FN_SSVS_n500_p100_AR1{c_index, ab_index,iters} = FN_SSVS;
    allOmega_n500_p100_AR1{c_index, ab_index,iters} = allOmega;
    allW_n500_p100_AR1{c_index, ab_index,iters} = allW;
    allWeights_n500_p100_AR1{c_index, ab_index,iters} = allWeights;
    allSig_n500_p100_AR1{c_index, ab_index,iters} = allSig;
    monitor_omega_chain_firsthalf_n500_p100_AR1{c_index, ab_index,iters} = monitor_omega_chain_firsthalf;
    monitor_omega_chain_secondhalf_n500_p100_AR1{c_index, ab_index,iters} = monitor_omega_chain_secondhalf;

        end
   end
        
end


save('Bsplines_paper_n500_p100_AR1_SpikeSlab_1to100_81to90.mat', '-v7.3');


    