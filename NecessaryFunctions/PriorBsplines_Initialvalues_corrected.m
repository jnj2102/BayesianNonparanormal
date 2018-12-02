function [optimalJ, minK, Final_AIC,F_mat_cell,g_vec_cell,...
Fconstraint_cell,RHS_gq_cell,A_mat_cell,LHS_gq_cell,W_mat_cell,...
q_vec_cell,knot_vector_cell,index_1_cell,index_2_cell,...
inverse_variance_prior_reduced_cell,mean_prior_reduced_cell,...
Z_red_cell,Z_two_cell,initial_value] = PriorBsplines_Initialvalues_corrected(n,p, x_matrix, mu, tau, sigma2,c,...
N_points)

%Find the optimal number of B-spline functions, the prior mean and variance
%and the initial B-spline coefficients for the transformation function
%estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
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

  initial_value{p} = [];
  
  
    [optimalJ, minK, Final_AIC] = AICloop_updated(n,p, x_matrix, mu, tau, sigma2,c);

  
for d = 1:p
%First calculate the elements of the truncated normal prior for the B-splines
[F_mat, g_vec, W_mat, q_vec, Fconstraint, A_mat, knot_vector, index_1, index_2,inverse_variance_prior_reduced,...
    mean_prior_reduced, LHS_gq, RHS_gq, Z_red, Z_two] = prior_Bsplines_corrected(optimalJ, minK, mu, tau, sigma2,...
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

  
 
  for d=1:p
      
 initial_final_tmp = quadprog(RHS_gq_cell{d}'*RHS_gq_cell{d},...
     -LHS_gq_cell{d}*RHS_gq_cell{d}, -Fconstraint_cell{d}, -repmat(1e-4, [optimalJ(d)-1, 1]), A_mat_cell{d}, c);

initial_final = initial_final_tmp;
initial_final([index_1_cell{d},index_2_cell{d}]) = []; 

         if  all(F_mat_cell{d} *initial_final + g_vec_cell{d} <= 0) 
    
         fprintf('The initials for pred %d Chain %d did not converged', d);
  
         end
         
      initial_value{d} = initial_final;

  end


    










end