%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Bayesian nonparanormal graphical model 
% Looking at different choices of n, p, and sparsity
%
%For first paper on Bsplines
% Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%add the paths to where the files are located
%addpath ./DEMO/  %SSVS function


clear; %clear the workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation combination: n500, p100, sparsity = AR1

%set the dimension
p=100;
%set the sample size
n=500;

%these are all true omegas (omega_true) - I will need to loop through the models
[star_model, AR2_model, AR4_model, circle_model, full_model, AR1_model, band_model,...
    fivepercent_model,tenpercent_model, twentypercent_model, twopercent_model, fifteenpercent_model] = model(p);
    
    
%set the true omega
omega_true = AR1_model;
sigma_true = inv(omega_true);

%create the true means. 
mean_true = transpose( linspace(1,2,p)); 

mu = 1; %I'm just giving it a mean of 1 for now but it can be any constant
tau = 0.5; %I'm giving it a sd of 0.5 for now but it can be any constant

sigma2 = 1; %I'm making the variance 1 for now but it can be any constant

%set the number of points for the Gauss-Hermite quadrature to 20

N_points = 20;

chains = 3;

c = [0;1]; %vector of linear constraints

%save the sigma matrix for R
save('Sigma_true_n500_p100_AR1_1to101.mat','sigma_true');

%save the omega matrix for R
save('Omega_true_n500_p100_AR1_1to101.mat','omega_true');

%set the number of replications for the simulation
reps = 100;

L1_n500_p100_AR1 = cell([reps,1]);
    SP_n500_p100_AR1 = cell([reps,1]);
    SE_n500_p100_AR1 = cell([reps,1]);
    MCC_n500_p100_AR1 = cell([reps,1]);
    f_sup_n500_p100_AR1 = cell([reps,1]);
    q_iters_n500_p100_AR1= cell([reps,1]);
    mean_est_n500_p100_AR1 = cell([reps,1]);
    final_thetas_n500_p100_AR1 = cell([reps,1]);
    burnin_theta_n500_p100_AR1= cell([reps,1]);
    burnin_omega_n500_p100_AR1 = cell([reps,1]);
    y_transformed_n500_p100_AR1 = cell([reps,1]);
    y_true_n500_p100 = cell([reps,1]);

    optimalJ_n500_p100_AR1 = cell([reps,1]);
    minK_n500_p100_AR1 = cell([reps,1]);
    Final_AIC_n500_p100_AR1 = cell([reps,1]);
	total_time_n500_p100_AR1 =  cell([reps,1]);
    

 x_matrix_n500_p100 =  cell([reps,1]);
    true_f_n500_p100 = cell([reps,1]);
    gridpoint_n500_p100 = cell([reps,1]);
    
    All_tmp_iters_n500_p100_AR1 = cell([reps,1]);
    complete_Sig_total_sims_n500_p100_AR1 = cell([reps,1]);
      complete_Z_total_sims_n500_p100_AR1 = cell([reps,1]);
        complete_Omega_total_sims_n500_p100_AR1 = cell([reps,1]);
        
        
F_mat_cell_iters = cell([reps,1]);
g_vec_cell_iters = cell([reps,1]);
A_mat_cell_iters = cell([reps,1]);
LHS_gq_cell_iters = cell([reps,1]);
W_mat_cell_iters = cell([reps,1]);
q_vec_cell_iters = cell([reps,1]);
knot_vector_cell_iters = cell([reps,1]);
index_1_cell_iters = cell([reps,1]);
index_2_cell_iters = cell([reps,1]);
inverse_variance_prior_reduced_cell_iters = cell([reps,1]);
mean_prior_reduced_cell_iters = cell([reps,1]);
Fconstraint_cell_iters = cell([reps,1]);
RHS_gq_cell_iters = cell([reps,1]);
Z_red_cell_iters = cell([reps,1]);
Z_two_cell_iters = cell([reps,1]);

MLE_norm__n500_p100 = cell([reps,1]);
    MLE_log__n500_p100 = cell([reps,1]);
    MLE_ev_n500_p100 = cell([reps,1]);
    MLE_stable_n500_p100 = cell([reps,1]);
    initial_value_cell_iters = cell([reps,1]);
    seed = cell([reps,1]);

F_mat_cell = cell([p,1]);
g_vec_cell = cell([p,1]);
A_mat_cell = cell([p,1]);
LHS_gq_cell = cell([p,1]);
W_mat_cell = cell([p,1]);
q_vec_cell = cell([p,1]);
knot_vector_cell = cell([p,1]);
index_1_cell = cell([p,1]);
index_2_cell = cell([p,1]);
inverse_variance_prior_reduced_cell = cell([p,1]);
mean_prior_reduced_cell = cell([p,1]);
Fconstraint_cell = cell([p,1]);
RHS_gq_cell = cell([p,1]);
Z_red_cell = cell([p,1]);
Z_two_cell = cell([p,1]);

  initial_tmp{p} = [];
  initial_value_cell{chains} = [];

%Run the iterations

for iters = 1:reps
    
    fprintf('iters');
    
    rng(iters,'twister'); %set the seed for each replication for reproducibility
   
    % generate the true y using the AR1 model
    y_true = mvnrnd(mean_true,sigma_true,n); 
    
    
    %Save the true y
    y_true_n500_p100{iters} = y_true;
     
   
    [x_matrix, true_f, gridpoint, MLE_norm, MLE_ev, MLE_log, MLE_stable] = xfunctions(p, y_true, n);
    
    x_matrix_n500_p100{iters} = x_matrix;
    true_f_n500_p100{iters} = true_f;
    gridpoint_n500_p100{iters} = gridpoint;
    MLE_norm__n500_p100{iters} = MLE_norm;
    MLE_log__n500_p100{iters} = MLE_log;
    MLE_ev_n500_p100{iters} = MLE_ev;
    MLE_stable_n500_p100{iters} = MLE_stable;
    

    %save the x matrix as a mat file for R
    
    save(sprintf('Bsplines_Iter_%d_x_matrix_n500_p100_AR1_1to100.mat', iters), 'x_matrix');
       
    %First find the optimal J for the replication
    
  
    [optimalJ, minK, Final_AIC] = AICloop(n,p, x_matrix, mu, tau, sigma2);
    
    
    optimalJ_n500_p100_AR1{iters} = optimalJ;
    minK_n500_p100_AR1{iters} = minK;
    Final_AIC_n500_p100_AR1{iters} = Final_AIC;
    
    %Initialize the cells again for the dth predictor

    
        
for d = 1:p
%First calculate the elements of the truncated normal prior for the B-splines
[F_mat, g_vec, W_mat, q_vec, Fconstraint, A_mat, knot_vector, index_1, index_2,inverse_variance_prior_reduced,...
    mean_prior_reduced, LHS_gq, RHS_gq, Z_red, Z_two] = prior_Bsplines(optimalJ, minK, mu, tau, sigma2,...
  N_points, c, x_matrix, d);

F_mat_cell{d} = F_mat;
g_vec_cell{d} = g_vec;
Fconstraint_cell{d} = Fconstraint;
RHS_gq_cell{d}= RHS_gq;
A_mat_cell{d} = A_mat;
LHS_gq_cell{d} = LHS_gq;
W_mat_cell{d} = W_mat;
q_vec_cell{d} = double(q_vec);
knot_vector_cell{d} = knot_vector;
index_1_cell{d} = index_1;
index_2_cell{d} = index_2;
inverse_variance_prior_reduced_cell{d} = inverse_variance_prior_reduced;
mean_prior_reduced_cell{d} = mean_prior_reduced;
Z_red_cell{d} = Z_red;
Z_two_cell{d}  = Z_two;


end  %end of predictors loop



%Find the initial values for the HMC sampler for each chain

  
 dispersion = randperm(3,3);
 
 
  for C = 1:chains  
  
  for d=1:p
      
 initial_final_tmp = quadprog((dispersion(C)*RHS_gq_cell{d})'*(dispersion(C)*RHS_gq_cell{d}),...
     -LHS_gq_cell{d}*(dispersion(C)*RHS_gq_cell{d}), -Fconstraint_cell{d}, -repmat(1e-4, [optimalJ(d)-1, 1]), A_mat_cell{d}, c);

initial_final = initial_final_tmp;
initial_final([index_1_cell{d},index_2_cell{d}]) = []; 

         if  all(F_mat_cell{d} *initial_final + g_vec_cell{d} <= 0) 
    
         fprintf('The initials for pred %d Chain %d did not converged', d, C);
  
         end
         
      initial_tmp{d} = initial_final;

  end

        initial_value_cell{C} = initial_tmp;
  end
    
initial_value_cell_iters{iters} = initial_value_cell;

%I have these iters to run in my matlab function for each iters
F_mat_cell_iters{iters} = F_mat_cell;
g_vec_cell_iters{iters} = g_vec_cell;
Fconstraint_cell_iters{iters} = Fconstraint_cell;
RHS_gq_cell_iters{iters} = RHS_gq_cell;
A_mat_cell_iters{iters} = A_mat_cell;
LHS_gq_cell_iters{iters} = LHS_gq_cell;
W_mat_cell_iters{iters} = W_mat_cell;
q_vec_cell_iters{iters} = q_vec_cell;
knot_vector_cell_iters{iters} = knot_vector_cell;
index_1_cell_iters{iters} = index_1_cell;
index_2_cell_iters{iters} = index_2_cell;
inverse_variance_prior_reduced_cell_iters{iters} = inverse_variance_prior_reduced_cell;
mean_prior_reduced_cell_iters{iters} = mean_prior_reduced_cell;
Z_red_cell_iters{iters} = Z_red_cell;
Z_two_cell_iters{iters}  = Z_two_cell;



end %end of iters loop


save('Bsplines_paper_Sim_n500_p100_AR1_prior_1to100.mat');

    