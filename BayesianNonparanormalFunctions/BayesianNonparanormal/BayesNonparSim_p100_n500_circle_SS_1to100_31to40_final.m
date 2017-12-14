% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n=500, p=100, sparsity = circle

%Load just the variables that are needed because these
%variables take up a lot of memory

load('Bsplines_paper_n500_p100_circle_SpikeSlab_1to100_31to40.mat', 'Omega_L1_SSVS_n500_p100_circle',...
    'edge_matrix_n500_p100_circle', 'complete_Z_Bayes_est_sims_n500_p100_circle',...
    'SP_matrix', 'SE_matrix', 'MCC_matrix', 'total_time_n500_p100_circle','q_iters_n500_p100_circle',...
    'monitor_omega_chain_firsthalf_n500_p100_circle','monitor_omega_chain_secondhalf_n500_p100_circle',...
    'All_tmp_iters_n500_p100_circle','burnin_n500_p100_circle','num_elements1', 'num_elements2',...
    'p','n');



for iters = 31:40
	
    rng(iters,'twister'); %set the seed for each replication for reproducibility

 
   fprintf('Iterations = %d', iters);

    %Now do the B-splines estimation method    
    
    for c_index = 1:num_elements1
        for ab_index = 1:num_elements2
               
        Omega_L1_SSVS = Omega_L1_SSVS_n500_p100_circle{c_index, ab_index,iters};
        
        edge_matrix_ssvs = edge_matrix_n500_p100_circle{c_index, ab_index,iters};
        
       complete_Z_Bayes_est_sims = complete_Z_Bayes_est_sims_n500_p100_circle{c_index, ab_index,iters};
       
       %Find the MLE using a convex optimization problem

%find nonzero elements of the diagonal of the upper triangular part of omega 
Omega_L1_SSVS_triu = triu(Omega_L1_SSVS);

%we are including all of the diagonal entries.
[diag_index, ~] = find(diag(Omega_L1_SSVS_triu) ~= 0);
%find the nonzero elements that relate to the offdiagonal nonzero elements of the strict upper
%triangular part of omega
edge_matrix_ssvs_strict_triu = triu(edge_matrix_ssvs, 1);

[row_index, ~] = find(edge_matrix_ssvs_strict_triu(:) ~=0);

%convert the index to subscripts.  You have to tell it the dimensions of
%the matrix.
[row_omega, col_omega] = ind2sub(size(edge_matrix_ssvs),row_index);

number_nonzero_E_1 = numel([diag_index;row_omega]); %number of nonzero elements = q
number_nonzero_E_2 = numel([diag_index;col_omega]); 

nonzero_E_1 = [diag_index;row_omega];
nonzero_E_2 = [diag_index;col_omega];

%create the E_1 matrix, a p x q matrix

E_1 = zeros([p, number_nonzero_E_1]);
E_2 = zeros([p, number_nonzero_E_2]);

%loop through row_omega because the loop will happen in order and tell me
%the column number.

for q = 1:number_nonzero_E_1
    row = nonzero_E_1(q);
    col = nonzero_E_2(q);
    
    E_1(row, q) = 1;
    E_2(col, q) = 1;
    
end


%find the average sample covariance matrix
mean_Z_Bayes_est = mean(complete_Z_Bayes_est_sims,3);
mean_sample_covariance_bayes_est = 1/n*(mean_Z_Bayes_est'*mean_Z_Bayes_est); 

%Write an anonymous function that calculates the objective.

fun = @(x)-log(det((E_1*diag(x)*E_2' + E_2*diag(x)*E_1')))  + trace((E_1*diag(x)*E_2' + E_2*diag(x)*E_1')*mean_sample_covariance_bayes_est);


%my initial numbers must satisfy det(*) > 0
x_MLE = fminunc(fun,diag(E_1'*E_2),optimset('MaxFunEvals',50000,...
   'MaxIter',50000));  

% MLE for Omega.  Leave the diagonals as they are
Omega_MLE = E_1*diag(x_MLE)*E_2' + E_2*diag(x_MLE)*E_1';

%Now calculate the BIC and eBIC based on the MLE

%Not multiplying by n because using the sample covariance 1/n*Z'*Z

%find the number of nonzero entries in the upper diagonal portion of
%the estimated precision matrix MLE
Omega_MLE_upperdiag = triu(Omega_MLE);

num_upperdiag_omega_MLE = sum(Omega_MLE_upperdiag(:) ~= 0); %this should be the same as before since it was
%the constraint.
BIC = n*(-log(det(Omega_MLE)) +...
    trace(Omega_MLE*(mean_sample_covariance_bayes_est))) +...
    num_upperdiag_omega_MLE*log(n);

eBIC = n*(-log(det(Omega_MLE)) +...
    trace(Omega_MLE*(mean_sample_covariance_bayes_est))) +...
    num_upperdiag_omega_MLE*log(n) + 4*1*num_upperdiag_omega_MLE*log(p);

eBIC_half = n*(-log(det(Omega_MLE)) +...
    trace(Omega_MLE*(mean_sample_covariance_bayes_est))) +...
    num_upperdiag_omega_MLE*log(n) + 4*1/2*num_upperdiag_omega_MLE*log(p);

     BIC_matrix(c_index, ab_index,iters) = BIC;
     eBIC_matrix(c_index, ab_index,iters) = eBIC;
     eBIC_half_matrix (c_index, ab_index,iters) = eBIC_half;

        end
    end
end

    
    %Find the minimum BIC and keep the data that pertains to it for further
    %analysis
for iters = 31:40

[minBIC, ~ ] = min(min(BIC_matrix(:,:, iters)));

[rowMinBIC, colMinBIC] = find(BIC_matrix(:,:,iters) == minBIC);


%save the data with the minimum BIC for further analysis

%Variables in comments are not being called at this time.
 %L1_finalanalysis(iters) = L1_matrix(rowMinBIC, iters);
 SP_matrix_finalanalysis(iters)  =  SP_matrix(rowMinBIC, colMinBIC,iters);
 SE_matrix_finalanalysis(iters) =  SE_matrix(rowMinBIC, colMinBIC,iters);
 MCC_matrix_finalanalysis(iters) =  MCC_matrix(rowMinBIC, colMinBIC,iters);
BIC_matrix_finalanalysis(iters) =   BIC_matrix(rowMinBIC, colMinBIC,iters);
eBIC_matrix_finalanalysis(iters) =     eBIC_matrix(rowMinBIC, colMinBIC,iters);
 eBIC_half_matrix_finalanalysis(iters) =    eBIC_half_matrix(rowMinBIC, colMinBIC,iters);
    
 total_time_n500_p100_circle_finalanalysis(iters) =   total_time_n500_p100_circle(rowMinBIC, colMinBIC,iters);

 edge_matrix_n500_p100_circle_finalanalysis{iters} =   edge_matrix_n500_p100_circle{rowMinBIC, colMinBIC,iters};
Omega_L1_SSVS_n500_p100_circle_finalanalysis{iters} = Omega_L1_SSVS_n500_p100_circle{rowMinBIC, colMinBIC,iters};
   q_iters_n500_p100_circle_finalanalysis(iters) =  q_iters_n500_p100_circle(rowMinBIC, colMinBIC,iters);
%   final_thetas_n500_p100_circle_finalanalysis{iters} =  final_thetas_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%   final_mean_n500_p100_circle_finalanalysis{iters} = final_mean_n500_p100_circle{rowMinBIC, colMinBIC,iters};
    burnin_n500_p100_circle_finalanalysis(iters) = burnin_n500_p100_circle(rowMinBIC, colMinBIC,iters);
 %  y_transformed_n500_p100_circle_finalanalysis{iters} = y_transformed_n500_p100_circle{rowMinBIC, colMinBIC,iters};
   All_tmp_iters_n500_p100_circle_finalanalysis{iters} =  All_tmp_iters_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%  complete_Sig_total_sims_n500_p100_circle_finalanalysis{iters} =  complete_Sig_total_sims_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%   complete_Z_total_sims_n500_p100_circle_finalanalysis{iters} =  complete_edge_matrix_total_sims_n500_p100_circle{rowMinBIC, colMinBIC,iters};
 % complete_Omega_total_sims_n500_p100_circle_finalanalysis{iters} =   complete_Omega_total_sims_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%complete_Z_Bayes_est_sims_n500_p100_circle_finalanalysis{iters} = complete_Z_Bayes_est_sims_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%   TP_SSVS_n500_p100_circle_finalanalysis(iters) =  TP_SSVS_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%  TN_SSVS_n500_p100_circle_finalanalysis(iters) =  TN_SSVS_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%  FP_SSVS_n500_p100_circle_finalanalysis(iters) =  FP_SSVS_n500_p100_circle{rowMinBIC, colMinBIC,iters};
%   FN_SSVS_n500_p100_circle_finalanalysis(iters) = FN_SSVS_n500_p100_circle{rowMinBIC, colMinBIC,iters};
  %allOmega_n500_p100_circle_finalanalysis{iters} =  allOmega_n500_p100_circle{rowMinBIC, colMinBIC,iters};
 % allSig_n500_p100_circle_finalanalysis{iters} =   allSig_n500_p100_circle{rowMinBIC, colMinBIC,iters};
  monitor_omega_chain_firsthalf_n500_p100_circle_finalanalysis{iters} =  monitor_omega_chain_firsthalf_n500_p100_circle{rowMinBIC, colMinBIC,iters};
  monitor_omega_chain_secondhalf_n500_p100_circle_finalanalysis{iters} =  monitor_omega_chain_secondhalf_n500_p100_circle{rowMinBIC, colMinBIC,iters};



end

save('BayesNonparSim_n500_p100_circle_SS_1to100_31to40_final_limited.mat', '-v7.3');
