function [Bayes_est] = updateY_basis(basis_red,basis_two, W_mat_cell, q_vec_cell,...
                            new_sample)

 
%average the thetas across MCMCs to get theta_hat, the estimated thetas for
%each J. 

%so using the estimated thetas and the estfunction, get the Bayes
%estimate for each predictor.  The Bayes estimator (a function) uses the estimated
%thetas evaluated at x.  For the x, using the original x values.
    
   [Bayes_est] = estfunction(new_sample,...  
                   basis_red, basis_two, W_mat_cell, q_vec_cell);
               
             

end