function [edge_matrix_ssvs, q_iters, complete_thetas_total_sims, final_burnin , Y_transformed_matrix,...
    All_tmp_iters, total_time_q_iters, ...
    complete_Sig_total_sims, complete_edge_matrix_total_sims, complete_Omega_total_sims, ...
    complete_mean_sims, allOmega, allW, allSig, allWeights,...
    monitor_omega_chain_firsthalf, monitor_omega_chain_secondhalf,Omega_L1_SSVS,complete_Z_Bayes_est_sims] = BayesianNonparanormal_joint_realdata_clean(n,p,...
     x_matrix, F_mat_cell, g_vec_cell, W_mat_cell, q_vec_cell,...
	knot_vector_cell, index_1_cell, index_2_cell, inverse_variance_prior_reduced_cell, mean_prior_reduced_cell,...
    Z_red_cell, Z_two_cell, c0, a,b,f,g, initial_value)

%Function to run the B-splines approach with the SSVS to the Bayesian nonparanormal
%graphical model method.  Does BDA version of Rhat with split chains. Joint
%algorithm version.  For use with real data.

%Author: Jami Jackson
%
%Input:
%Output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


indmx = reshape(1:p^2,p,p); 
  upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
 
%Values for Convergence Diagnostics

chains = 3;
T_iters = 15000;   
iters = 2*T_iters;
a_iters = max(100, floor(iters/100));
Q_iters = iters/a_iters;


%Pre-allocate for the loops
sims_chain{chains} = [];
mean_chain = cell([chains,1]);
sigma_chain = cell([chains,1]);
W_chain_cell = cell([chains,1]);
W_chain{p} = [];
samples_cell_tmp = cell([p,1]);
allOmega = cell([chains,1]);
allW = cell([chains,1]);
allSig = cell([chains,1]);
allWeights = cell([chains,1]);
W_chain_mean = cell([chains,1]);
Y_chain = cell([chains,1]);
Z_chain = cell([chains,1]);
omega_chain = cell([chains,1]);
Z_Bayes_est_matrix_saved_chain = cell([chains,1]);

mean_monitor_omega_chain_firsthalf = cell([chains,1]);
mean_monitor_omega_chain_secondhalf = cell([chains,1]);
within_variance_monitor_omega_chain_firsthalf = cell([chains, 1]);
within_variance_monitor_omega_chain_secondhalf = cell([chains, 1]);
monitor_omega_chain_firsthalf = cell([chains, 1]);
monitor_omega_chain_secondhalf = cell([chains, 1]);


All_tmp_iters = cell(0);

lambda = 1; 
max_while = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Evaluate the function 

basis_full = cell([p,1]);
basis_red = cell([p,1]); 
basis_two = cell([p,1]);
  
Bayes_est_cell_initial = zeros([n,p]);
Bayes_est_cell_initial_chain = cell([chains,1]);
mean_initial_chain = cell([chains,1]);
sigma_initial_chain = cell([chains,1]);
omega_initial_chain = cell([chains,1]);
    

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

for C = 1:chains    
    
    for d = 1:p
%average the thetas across MCMCs to get theta_hat, the estimated thetas for
%each J. 

%so using the estimated thetas and the estfunction, get the Bayes
%estimate for each predictor.  The Bayes estimator (a function) uses the estimated
%thetas evaluated at x.  For the x, using the original x values.
    
%the final sample is already the dth final sample
   [Bayes_est] = estfunction(mean(initial_value{C}{d}, 2),...  
                   basis_red{d}, basis_two{d}, W_mat_cell{d}, q_vec_cell{d});
               
               Bayes_est_cell_initial(:,d)  = Bayes_est;
    end
Bayes_est_cell_initial_chain{C} = Bayes_est_cell_initial;

mean_initial_chain{C} = mean(Bayes_est_cell_initial)';

sigma_initial_chain{C} = cov(Bayes_est_cell_initial);

omega_initial_chain_temp = invChol_mex(sigma_initial_chain{C}); 
                                                
%when p=n or most likely when p>n, the matrix may not be positive definite,
%so use the nearest positive definite matrix since these are just
%initial values anyway for omega only.                                                
if all(eig(omega_initial_chain_temp) > 0) == 1
    
     omega_initial_chain{C} = omega_initial_chain_temp;
else
    omega_initial_chain{C} = nearestSPD(omega_initial_chain_temp);
end                                                
                     
 end %end of chains

samples_cell = cell([p,1]);


%draw the initial value for tau^2
tau_initial = 1/gamrnd(f,1/g);

%draw an initial pii
pii_initial = betarnd(a,b);

tau = tau_initial*ones(p);
pii = pii_initial*ones(p);
sl_variance = tau; %these are the spike and slab variances for omega

W = zeros(p,p);

tic

for q_iters = 1:Q_iters


    burnin = floor(q_iters*a_iters/2);
        
    
    num_samples = q_iters*a_iters ; 

        saved_samples = num_samples-burnin;

    
     samples_cell_iters = cell([num_samples,1]);

     mean_posterior_iters = zeros([p,num_samples]);  

     Omega_sims_chain_tmp = zeros([p,p,num_samples]); 
    W_sims_chain_tmp = zeros([p,p,num_samples]);

    Sig_sims_chain_tmp = zeros([p,p,num_samples]);
    Bayes_est_z_iters =  zeros([n,p,num_samples]);
    Bayes_est_Y_iters = zeros([n,p,num_samples]);

    weights_sims_chain_tmp = zeros([p,p,num_samples]) ;

final_burnin = q_iters*a_iters;
    
for C = 1:chains 
    
    theta_initial = initial_value{C};
    
    omega_initial = omega_initial_chain{C};
    sigma_initial = sigma_initial_chain{C};
    
    
    mean_initial = mean_initial_chain{C};
    
    Y_initial = Bayes_est_cell_initial_chain{C};

%Step 5 Truncated Multivariate Normal sampling of Marginals using HMC


for iters_num = 1:num_samples 

for d= 1:p
    
    %First find the initial values
     
  if  all(F_mat_cell{d} *theta_initial{d} + g_vec_cell{d} < 0) %Make this <0 to say it doesn't satisfy the constraint

        disp('The initial value does not satisfy the constraint.');
 
   end
   
   
   ind_noi_index = ind_noi_all(:,d);
   
   %mean and sigma are not getting updated at this step because they are
   %random. Only Y is getting updated.
   

  
 [M_mat, r_vec] = posterior_Bsplines_Gibbs_joint(d,ind_noi_index, W_mat_cell{d},...
     q_vec_cell{d}, mean_initial, omega_initial, inverse_variance_prior_reduced_cell{d},...
     mean_prior_reduced_cell{d}, n, Z_red_cell{d}, Z_two_cell{d}, Y_initial);
 
 %Now take these precision matrix and linear term and use it for the HMC.

sM = sparse(M_mat);
sr = sparse(r_vec);
sF = sparse(F_mat_cell{d});
sg = sparse(g_vec_cell{d});

%the HMC algorithm as it is written by the author outputs the initial values as
%the first sample.  So it will be ignored.
    [total, ~] = HMC_exact(sF,sg, sM, sr', false, 2, theta_initial{d} );  %get two samples

    samples_iters = total(:, 2); % I only want the second sample, not the first
     
%remember these thetas are different sizes, that's why I need to use a cell
    samples_cell_tmp{d}  = samples_iters;  %save these for the chain for rhat later after burnin

    % Now Y_initial gets updated by calculating a new Y
    %column.  
    
   Y_new_column =  updateY_basis(basis_red{d}, basis_two{d}, W_mat_cell{d}, q_vec_cell{d},...
        samples_iters);
    
    Y_initial(:,d) = Y_new_column;  %now do the next sample with this new Y column


end %end of d=1:p

%I need to update my thetas with all of the last drawn thetas.
theta_initial = samples_cell_tmp;


[Z,mean_posterior] = samplemunew_Gibbs(Y_initial, p, n, sigma_initial);


%update the mean_initial with the mean_posterior
mean_initial = mean_posterior;

samples_cell_iters{iters_num} = samples_cell_tmp;

%Now sample mu 
     mean_posterior_iters(:,iters_num) = mean_posterior;  

     
Bayes_est_z_iters(:,:, iters_num) = Z;

Bayes_est_Y_iters(:,:,iters_num) = Y_initial; 

%create the sum of products matrix 
  
S = Z'*Z;

S_ExpCorrC = zeros(p,p);
  
for j = 1:p
    for k = 1:p
        S_ExpCorrC(j,k) = S(j,k)/(sqrt(S(j,j))*sqrt(S(k,k)));
    end
end


 S_std = S_ExpCorrC*n; %standardize the S

 Sig = sigma_initial;
  
C_omega = invChol_mex(Sig);  %Following the code of (Wang, 2015, Scaling it Up)

L = 1;
    
while (L <= max_while)

for i = 1:p
    
      ind_noi = ind_noi_all(:,i);
 
      sl_temp = sl_variance(ind_noi,i);
      
      Sig11 = Sig(ind_noi,ind_noi);

      Sig12 = Sig(ind_noi,i);

      invC11 = Sig11 - Sig12*Sig12'/Sig(i,i);

      Ci = (S_std(i,i)+lambda)*invC11+diag(1./sl_temp);
      
       Ci = (Ci+Ci')./2;  
       
       Ci_chol = chol(Ci);   
       
      
      mu_i = -Ci_chol\(Ci_chol'\S_std(ind_noi,i));
      
      beta = mu_i+ Ci_chol\randn(p-1,1);
        
        
          %update Omega
        C_omega(ind_noi,i) = beta;
        C_omega(i,ind_noi) = beta;
        
        a_gam = 0.5*n+1;
        b_gam = (S_std(i,i)+lambda)*0.5;
        gam = gamrnd(a_gam,1/b_gam);
      
        
        c_diag = beta'*invC11*beta;
        C_omega(i,i) = gam+c_diag;
      
        invC11beta = invC11*beta;
                
        %update Sigma
        Sig(ind_noi,ind_noi) = invC11+invC11beta*invC11beta'/gam;
        Sig12 = -invC11beta/gam;
        Sig(ind_noi,i) = Sig12;
        Sig(i,ind_noi) = Sig12';
        Sig(i,i) = 1/gam;
        
        
         pii_vector = pii(ind_noi,i);     
         tau_squared = tau(ind_noi,i);
         
        w1 = -0.5*log(c0*tau_squared) -0.5*beta.^2./(c0*tau_squared)+log(1-pii_vector);
        w2 = -0.5*log(1*tau_squared) -0.5*beta.^2./(1*tau_squared)+log(pii_vector);
        
        
         w_max = max([w1,w2],[],2);

            %Find the Pr(Z=1) 
        prob = exp(w2-w_max)./sum(exp([w1,w2]-repmat(w_max,1,2)),2);

        w = (rand(p-1,1)<prob);  %logical for Z = 1 based on weight w

        %convert z which is a logical array to a numeric array using double
        
        w = double(w);
                
         %change all zeros to c0.
        w(w==0) = c0;
                
        
        %update W 
        W(ind_noi,i) = w;
        W(i,ind_noi) = w;
                
        
       %update the pii
       %Have to use logicals because updating the whole column vector
        %if W_vector = 1 then that's a 1 and that makes W_vector = c0, 0
        %automatically
       pii_update = betarnd(a + (w == 1), b + (w == c0),...
           [p-1,1]);
   
        
       %update pii
       pii(ind_noi,i) = pii_update;
       pii(i,ind_noi) = pii_update;
       
       %Lastly tau gets updated
       
       alpha_tau = (f + 1/2)*ones([p-1,1]); %constant vector
       beta_tau = 1./(g + (beta.^2./(2*w)));
       
       tau_update = 1./(gamrnd(alpha_tau, beta_tau));
         
       %update tau
       tau(ind_noi,i) = tau_update;
       tau(i,ind_noi) = tau_update;

       %update the spike and slab variances for omega
       
       sl_variance(ind_noi,i) = w.*tau_update;
       sl_variance(i,ind_noi) = w.*tau_update;
            
       
 end       
        

 test = eig(C_omega);
    
  if (all(test > 0) ==1) 
      break;
  else
       L = L + 1;
    continue;
  end
  
  
  
end %end of while loop

if L == max_while
    disp('Could not find a pd matrix');
end


%update sigma_initial 
omega_initial = C_omega;
sigma_initial = Sig;

Omega_sims_chain_tmp(:,:,iters_num) = C_omega; 
W_sims_chain_tmp(:,:,iters_num) = W;  

Sig_sims_chain_tmp(:,:,iters_num) = Sig;

weights_sims_chain_tmp(:,:,iters_num) = pii;


end %end of number of iterations

%Save these iterations to create trace plots later.
allOmega{C} = Omega_sims_chain_tmp;
allW{C} = W_sims_chain_tmp;
allSig{C} = Sig_sims_chain_tmp;
allWeights{C} = weights_sims_chain_tmp;


for d = 1:p
    
    samples_tmp = [];
    
for l = 1:saved_samples
    
samples_tmp(:,l) = samples_cell_iters{burnin+l}{d};
end

samples_cell{d} = samples_tmp;

end

sims_chain{C} = samples_cell;

%I need W_chain for each predictor. 

for d=1:p
W_chain{d} = width( samples_cell{d}', .2);
end

W_chain_cell{C} = W_chain;

%Next find the within-chain variability after doing burnin 

mean_saved =  mean_posterior_iters(:, burnin+1:end);

mean_chain{C} = mean_saved;


%now find the within variability 

W_chain_mean{C} = width(mean_saved', .2);

%Do burnin for the Y's

Y_saved = Bayes_est_Y_iters(:,:, burnin+1:end);
Y_chain{C} = Y_saved;

%Do the burnin for the Z's.  These Z's are Z ~ N(0, Sigma). Not the edge
%matrix.
Z_Bayes_est_matrix_saved = Bayes_est_z_iters(:,:,burnin+1:end);
Z_Bayes_est_matrix_saved_chain{C} = Z_Bayes_est_matrix_saved;

%Do burnin for omega

omega_saved = Omega_sims_chain_tmp(:,:, burnin+1:end); %these are the omegas I'm testing 
                                    %for Rhat

omega_chain{C} = omega_saved;

sigma_saved = Sig_sims_chain_tmp(:,:,  burnin+1:end);

sigma_chain{C} = sigma_saved;

%Lastly, do this for the W.

W_saved = W_sims_chain_tmp(:,:,burnin+1:end); %I'm keeping it as Z because I'm lazy

Z_chain{C} = W_saved;


%Find the log det(omega) and sum of edges
monitor_omega_firsthalf = zeros([1,saved_samples/2]);


for samps = 1:saved_samples/2

    edge_matrix_temp = W_saved(:,:,samps);
    sum_edges =  sum(edge_matrix_temp(upperind) == 1);
                monitor_omega_firsthalf(:,samps) = sum_edges;

 
end

    [number_elements, ~] = size(monitor_omega_firsthalf);

%Find the log det(omega) and sum of edges
monitor_omega_secondhalf = zeros([number_elements,saved_samples/2]);


for samps = 1:saved_samples/2
    samps_index = (saved_samples/2 + samps);

      edge_matrix_temp = W_saved(:,:,samps_index); % I want the second half
    sum_edges =  sum(edge_matrix_temp(upperind) == 1);
            monitor_omega_secondhalf(:,samps) = sum_edges;

end


%Find the average of each element in monitor_omega

mean_monitor_omega_firsthalf = mean(monitor_omega_firsthalf, 2);
mean_monitor_omega_secondhalf = mean(monitor_omega_secondhalf, 2);


mean_monitor_omega_chain_firsthalf{C} = mean_monitor_omega_firsthalf;
mean_monitor_omega_chain_secondhalf{C} = mean_monitor_omega_secondhalf;

%within variance


within_variance_monitor_omega_firsthalf = (1/((saved_samples/2)-1)*...
    sum((monitor_omega_firsthalf - repmat(mean_monitor_omega_firsthalf,[1,saved_samples/2])).^2,2))';


within_variance_monitor_omega_secondhalf = (1/((saved_samples/2)-1)*...
    sum((monitor_omega_secondhalf - repmat(mean_monitor_omega_secondhalf,[1,saved_samples/2])).^2,2))';

within_variance_monitor_omega_chain_firsthalf{C} = within_variance_monitor_omega_firsthalf;
within_variance_monitor_omega_chain_secondhalf{C} = within_variance_monitor_omega_secondhalf;

%save the values that I am monitoring for each chain as well - for the
%trace plots

monitor_omega_chain_firsthalf{C} = monitor_omega_firsthalf;
monitor_omega_chain_secondhalf{C} = monitor_omega_secondhalf;


end  

%END OF CHAINS%

    %Find RHAT using BDA way with the split chains.
        %stick the two halves together before taking the mean
    overallmean_monitor_omega = mean([cell2mat(mean_monitor_omega_chain_firsthalf'),...
        cell2mat(mean_monitor_omega_chain_secondhalf')] ,2);
    
    temp_sum_pieces_firsthalf = zeros([number_elements,chains]);
temp_sum_pieces_secondhalf = zeros([number_elements,chains]);

for C = 1:chains

   temp_sum_pieces_firsthalf(:, C)= (mean_monitor_omega_chain_firsthalf{C} - overallmean_monitor_omega).^2;
   temp_sum_pieces_secondhalf(:, C)= (mean_monitor_omega_chain_secondhalf{C} - overallmean_monitor_omega).^2;

end

%now average over 2 times the number of chains since putting them together.

%I could also change the orientation of the within_variance... to remove
%these transposes.
Within_monitor_omega = mean([cell2mat(within_variance_monitor_omega_chain_firsthalf);...
    cell2mat(within_variance_monitor_omega_chain_secondhalf)])';

%now we average over 2 times the number of chains and the number of samples
%is divided by 2
Between_monitor_omega = (saved_samples/2)/(2*chains-1)*sum([temp_sum_pieces_firsthalf,temp_sum_pieces_secondhalf],2);


variance_hat_plus = ((saved_samples/2)-1)/(saved_samples/2)*Within_monitor_omega +...
    1/(saved_samples/2)*Between_monitor_omega;

R_hat_BDA = sqrt(variance_hat_plus./Within_monitor_omega); 

%Now compute the variogram V_t at each lag t
%Compute a partial sum, starting from lag 0 and continuing until
%the sum of autocorrelation estimates for two successive lag is negative

rho_hat_t = zeros([number_elements,saved_samples/2]);
temp_inner_pieces_firsthalf = zeros([number_elements,saved_samples/2]);
temp_inner_pieces_secondhalf = zeros([number_elements,saved_samples/2]);

inner_sum_firsthalf = zeros([number_elements,chains]);
inner_sum_secondhalf = zeros([number_elements,chains]);

%I need to do this for each element of the monitor omega

for elem = 1:number_elements
    t = 0;  %Restart at zero every time.
    
while (t <= ((saved_samples/2)-1))
    
    for C = 1:chains
        monitor_omega_tmp_firsthalf = monitor_omega_chain_firsthalf{C};
        monitor_omega_tmp_secondhalf = monitor_omega_chain_secondhalf{C};

    
        for i = (t+1):(saved_samples/2)
        temp_inner_pieces_firsthalf(elem,i) = (monitor_omega_tmp_firsthalf(elem,i) -...
            monitor_omega_tmp_firsthalf(elem,i-t)).^2;
        
        temp_inner_pieces_secondhalf(elem,i) = (monitor_omega_tmp_secondhalf(elem,i) -...
            monitor_omega_tmp_secondhalf(elem,i-t)).^2;
   
        end
        
        inner_sum_firsthalf(elem,C) = sum(temp_inner_pieces_firsthalf(elem,:));
        inner_sum_secondhalf(elem,C) = sum(temp_inner_pieces_secondhalf(elem,:));

   
    end
    %sum across chains now.  The number of samples (n) is saved_samples/2
 
    V_t = 1/((saved_samples/2)-t) *mean([inner_sum_firsthalf(elem,:),inner_sum_secondhalf(elem,:)]);
    
    rho_hat_t(elem,t+1) = 1-(V_t/(2*variance_hat_plus(elem))); %I can't index with 0
    
    temp_index_1 = t;
    temp_index_2 = t+1;
    
    if  1 <= t && t <= ((saved_samples/2) - 1) && (rho_hat_t(elem,temp_index_1) + rho_hat_t(elem,temp_index_2)) < 0
        break;
    else
        t = t+1;
        continue;
    end
    
    
end

 
end


%effective sample size is based on 2*chains and saved_samples/2
effective_sample_size_hat = (2*chains*(saved_samples/2))./(1+2*sum(rho_hat_t,2));

All_tmp = horzcat(R_hat_BDA');
    
    All_tmp_iters{q_iters} = All_tmp; %save this to plot later
    
   disp([All_tmp,effective_sample_size_hat']);
   
   
   if (all(All_tmp < 1.05) && all(effective_sample_size_hat >= 100))
      
    break
    
   else
       
   %Update the initial values for theta with the total_sims
   
 % The code takes care of the chains and all of the mcmcs
   
        %Initialize the mean of the posterior using the samples 
        %averaged across each chain


for C = 1:chains
    
    for d = 1:p
initial_value{C}{d} = mean(sims_chain{C}{d},2);

    end
                        
Bayes_est_cell_initial_chain{C} = mean(Y_chain{C},3);
mean_initial_chain{C} = mean(mean_chain{C},2);
omega_initial_chain{C} = mean(omega_chain{C},3);
sigma_initial_chain{C} = mean(sigma_chain{C},3);

end


 %show the q_iters
   if q_iters == Q_iters
         disp('The algorithm did not converge.');
   end

continue
 
   end  
   

end

   %%END OF q_ITERS%%
   
 total_time_q_iters=  toc;
   
 %%Now instead of using these completes, just pull all of the simulations
 %%together for each chain.

 
%Pull together the thetas by chain to find the supremums

complete_thetas_total_tmp = cell([1, chains]);
complete_thetas_total_sims = cell([p,1]);

for d = 1:p
    
for C = 1:chains
    complete_thetas_total_tmp{C} = sims_chain{C}{d};
end
    complete_thetas_total_sims{d} = cell2mat( complete_thetas_total_tmp);
end


%%Evaluate the function 

%find the Basis_red and Basis_two for each replication
  
    Bayes_est_cell_total{p} = [];
    
for d = 1:p 
 
%average the thetas across MCMCs to get theta_hat, the estimated thetas for
%each J. 
%so using the estimated thetas and the estfunction, get the Bayes
%estimate for each predictor.  The Bayes estimator (a function) uses the estimated
%thetas evaluated at x.  For the x, using the original x matrix.
    
%the final sample is already the dth final sample
   [Bayes_est_complete] = estfunction(mean(complete_thetas_total_sims{d}, 2),...  
                   basis_red{d}, basis_two{d}, W_mat_cell{d}, q_vec_cell{d});
               
               Bayes_est_cell_total{d}  = Bayes_est_complete;
               
 end
       
   
Y_transformed_matrix = cell2mat(Bayes_est_cell_total);


complete_mean_sims =  cat(3, mean_chain{:});  
complete_Sig_total_sims = cat(3, sigma_chain{:});   %make a list to a 3 dimensional matrix (p x p x total samples) using cat
  
complete_edge_matrix_total_sims = cat(3, Z_chain{:});

complete_Omega_total_sims = cat(3, omega_chain{:});

%these are the Z ~ N(0, Sigma)
complete_Z_Bayes_est_sims = cat(3,Z_Bayes_est_matrix_saved_chain{:});

%Now calculate SE, SP, MCC, L1 and give the final estimated edge matrix


%find the loss function using SSVS

    Omega_L1_SSVS = invChol_mex(mean(complete_Sig_total_sims,3));  %Bayes estimator using SSVS

  
%Find the True Positives, True Negatives, False Positives, and False Negatives
% for SSVS method for 1 replication
 edge_matrix_ssvs = mean(complete_edge_matrix_total_sims,3)>.5;  %put .5 here because half the time we get 1, have the time 0
                                      %averaging the ones and zeros across the samples
    %even though there are .02 in the offdiagonals for edge matrix, this
    %mean gets rid of them.
           
end

