function[thetas_samples_cell,Y_initial] = NonparanormalBsplinesPrior(F_mat_cell,theta_initial,g_vec_cell,...
    ind_noi_all,W_mat_cell, q_vec_cell, mean_initial, omega_initial,...
    inverse_variance_prior_reduced_cell, mean_prior_reduced_cell, n, Z_red_cell,...
    Z_two_cell, Y_initial_matrix,basis_red, basis_two,p)
%Function to implement the random series B-splines prior for the
%transformation functions.
%
%Author: Jami Jackson Mulgrave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samples_cell_tmp = cell([p,1]);


for d= 1:p
    
    %First make sure the initial values satisfy the constraints
     
  if  all(F_mat_cell{d} *theta_initial{d} + g_vec_cell{d} < 0) 

        disp('The initial value does not satisfy the constraint.');
 
   end
   
   
   ind_noi_index = ind_noi_all(:,d);
   
   %mean and sigma are not getting updated at this step

  
 [M_mat, r_vec] = posterior_Bsplines_Gibbs_joint(d,ind_noi_index, W_mat_cell{d},...
     q_vec_cell{d}, mean_initial, omega_initial, inverse_variance_prior_reduced_cell{d},...
     mean_prior_reduced_cell{d}, n, Z_red_cell{d}, Z_two_cell{d}, Y_initial_matrix);

 %Now take these precision matrix and linear term and use it for the HMC.

sM = sparse(M_mat);
sr = sparse(r_vec);
sF = sparse(F_mat_cell{d});
sg = sparse(g_vec_cell{d});

%the HMC algorithm as it is written by the author has the initial values as
%the first sample.  So the first sample will be ignored.
    [total, ~] = HMC_exact(sF,sg, sM, sr', false, 2, theta_initial{d} );  %Obtain two samples

    samples_iters = total(:, 2); %Only keep the second sample, not the first
     
%These thetas are different sizes, so using cell arrays
    samples_cell_tmp{d}  = samples_iters;  %save these for the chain for rhat later after burnin


    % Now Y_initial gets updated 
    
   Y_new_column =  updateY_basis(basis_red{d}, basis_two{d}, W_mat_cell{d}, q_vec_cell{d},...
        samples_iters);

  
    Y_initial_matrix(:,d) = Y_new_column;  


end %end of d=1:p


thetas_samples_cell = samples_cell_tmp;

Y_initial = Y_initial_matrix;















end