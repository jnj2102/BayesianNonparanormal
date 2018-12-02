function [entropy_loss,SP_SSVS_total, SE_SSVS_total, MCC_SSVS_total, edge_matrix_ssvs,...
     total_time,TP_SSVS, TN_SSVS, FP_SSVS,FN_SSVS,Omega_Bayes_est,Sigma_Bayes_est,...
    mean_Z_Bayes_est,Frobenius_norm_covariance,bounded_loss,Frobenius_norm_precision] = BayesianNonparanormal_StudentTspikeslab(n,p, sigma_true,...
     x_matrix, omega_true, F_mat_cell, g_vec_cell, W_mat_cell, q_vec_cell,...
	knot_vector_cell, index_1_cell, index_2_cell, inverse_variance_prior_reduced_cell, mean_prior_reduced_cell,...
    Z_red_cell, Z_two_cell, c0, a,b,f,g, initial_value, lambda)

%Function to run the B-splines approach with the SSVS to the Bayesian nonparanormal
%graphical model method.  

%Author: Jami Jackson
%
%Input:
%Output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


indmx = reshape(1:p^2,p,p); 
  upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Evaluate the function 

basis_full = cell([p,1]);
basis_red = cell([p,1]); 
basis_two = cell([p,1]);
  
Bayes_est_cell_initial = zeros([n,p]);


ind_noi_all = zeros(p-1,p);

for i = 1:p
       if i==1  
       ind_noi = (2:p)'; 
      elseif i==p
       ind_noi = (1:p-1)'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       ind_noi_all(:,i) = ind_noi;
       
end
   

%%%%find the basis functions

for d = 1:p 

    basis_full{d} = bspline_basismatrix(4, knot_vector_cell{d}, x_matrix(:,d)); 

%remove the two rows of basis matrix using the two indices that match the
%thetas that we can solve for.

basis_red{d} = basis_full{d};
basis_red{d}(:, [index_1_cell{d},index_2_cell{d}]) = [];  

basis_two{d} = basis_full{d}(:, [index_1_cell{d},index_2_cell{d}]);
               
end


    
    for d = 1:p
%average the thetas across MCMCs to get theta_hat, the estimated thetas for
%each J. 

%so using the estimated thetas and the estfunction, get the Bayes
%estimate for each predictor.  The Bayes estimator (a function) uses the estimated
%thetas evaluated at x.  For the x, using the original x values.
    
%the final sample is already the dth final sample
   [Bayes_est] = estfunction(mean(initial_value{d}, 2),...  
                   basis_red{d}, basis_two{d}, W_mat_cell{d}, q_vec_cell{d});
               
               Bayes_est_cell_initial(:,d)  = Bayes_est;
    end

mean_initial_chain = mean(Bayes_est_cell_initial)';

sigma_initial_chain = cov(Bayes_est_cell_initial);

omega_initial_chain_temp = invChol_mex(sigma_initial_chain); 
                                                
%when p=n or most likely when p>n, the matrix may not be positive definite,
%so use the nearest positive definite matrix since these are just
%initial values anyway for omega only.                                                
if all(eig(omega_initial_chain_temp) > 0) == 1
    
     omega_initial_chain = omega_initial_chain_temp;
else
    omega_initial_chain = nearestSPD(omega_initial_chain_temp);
end                                                
                     

%draw the initial value for tau^2
tau_initial = 1/gamrnd(f,1/g);

%draw an initial pii
pii_initial = betarnd(a,b);

tau = tau_initial*ones(p);
pii = pii_initial; %now just one value
sl_variance = tau; %these are the spike and slab variances for omega

W = zeros(p,p);


burnin = 5000;
        

num_samples = 15000; 


    saved_samples = num_samples - burnin;
    
     complete_Omega_total_sims = zeros([p,p,saved_samples]); 

    complete_Sig_total_sims = zeros([p,p,saved_samples]);
    complete_Z_Bayes_est_sims =  zeros([n,p,saved_samples]);
 
    complete_edge_matrix_total_sims = zeros([p,p,saved_samples]) ;
    
    thetas_samples_cell = initial_value;
    
   Omega_sparse =  sparse(omega_initial_chain);
    Sig_sparse = sparse(sigma_initial_chain);
    
    mean_posterior = mean_initial_chain;
    
    Y_matrix = Bayes_est_cell_initial;

tic

for iters_num = 1:num_samples



%Step 5 Truncated Multivariate Normal sampling of Marginals using HMC


[thetas_samples_cell,Y_matrix] = NonparanormalBsplinesPrior(F_mat_cell,thetas_samples_cell,g_vec_cell,...
    ind_noi_all,W_mat_cell, q_vec_cell, mean_posterior, Omega_sparse,...
    inverse_variance_prior_reduced_cell, mean_prior_reduced_cell, n, Z_red_cell,...
    Z_two_cell, Y_matrix,basis_red, basis_two,p);




[Z,mean_posterior] = samplemunew_Gibbs(Y_matrix, p, n, Sig_sparse);


 [S_std] = SumProductsMatrix(Z,p,n);
 
 Sig = full(Sig_sparse);
  
C_omega = invChol_mex(Sig);  %Following the code of (Wang, 2015, Scaling it Up)


[C_omega,W, pii,Sig] = StudentTSpikeSlab(ind_noi_all,sl_variance,Sig,S_std,n,p,W,...
    C_omega,lambda,c0,pii,tau, a,b,f,g,upperind);




%update sigma_initial 
Omega_sparse = sparse(C_omega);

Omega_sparse = 1/2*(Omega_sparse + Omega_sparse'); %do this for numerical stability.

Sig_sparse = sparse(Sig);
Sig_sparse = 1/2*(Sig_sparse + Sig_sparse');

%Only save what's after burnin.

if iters_num > burnin
    %create the corresponding index
    
   save_index = iters_num - burnin;
     
complete_Z_Bayes_est_sims(:,:, save_index) = Z;

complete_Omega_total_sims(:,:,save_index) = Omega_sparse; 
complete_Sig_total_sims(:,:,save_index) = Sig_sparse;

complete_edge_matrix_total_sims(:,:, save_index) = W;

end


end %end of number of iterations

 total_time=  toc;
   
 
%Summaries
%Bayes estimate of Omega
Omega_Bayes_est = mean(complete_Omega_total_sims,3);

%Bayes estimate of Z transformed variables
mean_Z_Bayes_est = mean(complete_Z_Bayes_est_sims,3);

Sigma_Bayes_est = mean(complete_Sig_total_sims,3);  

%Find the loss


%Entropy loss
 entropy_loss = trace(Omega_Bayes_est*sigma_true) -...
     log(det(Omega_Bayes_est*sigma_true)) - p;
 
 

 
 %Frobenius norm

    %Find the Frobenius Norm 
    Frobenius_norm_precision = trace((Omega_Bayes_est -...
        omega_true)'*(Omega_Bayes_est -...
        omega_true));
    
        Frobenius_norm_covariance = trace((Sigma_Bayes_est -...
        sigma_true)'*(Sigma_Bayes_est -...
        sigma_true));
    
        %Find the bounded loss
    
    bounded_loss = 1/(p^2) * sum(sum(abs(Omega_Bayes_est -...
        omega_true)));
    
    
%Find the True Positives, True Negatives, False Positives, and False Negatives
% for SSVS method for 1 replication
 edge_matrix_ssvs = mean(complete_edge_matrix_total_sims,3)>.5;  %put .5 here because half the time we get 1, have the time 0
                                      %averaging the ones and zeros across the samples

                                      
%True edge matrix for the model 
edge_matrix_true = zeros(p,p); %1 is edge is present 0 is edge is absent.

for j = 1:p
    for k = 1:p
        if omega_true(j,k) == 0
            edge_matrix_true(j,k) = 0;
        else
            edge_matrix_true(j,k) = 1;
        end
    end
end


                                      
TP_matrix_SSVS = edge_matrix_ssvs(upperind) == 1 & edge_matrix_true(upperind) == 1;
TP_SSVS = sum(TP_matrix_SSVS(:)); %the colon sums all elements in the matrix


TN_matrix_SSVS = edge_matrix_ssvs(upperind) == 0 & edge_matrix_true(upperind) == 0;
TN_SSVS = sum(TN_matrix_SSVS(:)); 


FP_matrix_SSVS = edge_matrix_ssvs(upperind) == 1 & edge_matrix_true(upperind) == 0;
FP_SSVS = sum(FP_matrix_SSVS(:));


FN_matrix_SSVS = edge_matrix_ssvs(upperind) == 0 & edge_matrix_true(upperind) == 1;
FN_SSVS = sum(FN_matrix_SSVS(:)); 

%Find Specificity, Sensitivity, Matthews Correlation Coefficient 

SP_SSVS_total = TN_SSVS/(TN_SSVS + FP_SSVS);

SE_SSVS_total = TP_SSVS/(TP_SSVS + FN_SSVS);

MCC_SSVS_total = ((TP_SSVS*TN_SSVS) - (FP_SSVS*FN_SSVS))/sqrt((TP_SSVS+FP_SSVS)*(TP_SSVS+FN_SSVS)*(TN_SSVS+FP_SSVS)*(TN_SSVS+FN_SSVS));
  


end

