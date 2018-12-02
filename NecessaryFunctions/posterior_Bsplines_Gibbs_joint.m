function [M_mat, r_vec] = posterior_Bsplines_Gibbs_joint(d,ind_noi_index,W_mat, q_vec, mean_initial, omega_initial,...
    inverse_variance_prior_reduced, mean_prior_reduced, n, Z_red, Z_two, Y_initial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a function to build the posterior for B-splines version of the 
%Bayesian nonparanormal method for each predictor p in order to find the posterior precision
%matrix and the posterior linear term.  It is using conditonal normal
%theory to find the joint likelihood

%Author: Jami Jackson
%
%Input: d_index: the index for the predictor
%       x_matrix: the original x values
%       knot vector: the knot vector for the B-splines
%       index_1: the first index for the first theta that is solved for
%       index_2: the second index for the second theta that is solved for
%
%Output: M_mat: posterior precision matrix for the predictor p
%        r_vec: posterior linear term for the predictor p       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Build the first part of the posterior precision matrix, call it M1

  
   variance_likelihood = 1/omega_initial(d,d);

%M1 = cell([n,1]);

[row, col] = size((Z_red(1, :) + Z_two(1, :)*W_mat)'*(1/variance_likelihood)*(Z_red(1, :) + Z_two(1, :)*W_mat));

M1 = zeros([row,col,n]);

for i = 1:n
        M1(:,:,i) = (Z_red(i, :) + Z_two(i, :)*W_mat)'*(1/variance_likelihood)*(Z_red(i, :) + Z_two(i, :)*W_mat);
end
      

  %B_mat = cat(3, M1{:}); %concatenate the matrices along the 3rd dimension to 
  %maintain the form of the matrices
  
  C_mat = sum(M1,3); %sum along the 3rd dimension to add the n matrices
  
  %Final posterior precision matrix
  
  M_mat = C_mat + inverse_variance_prior_reduced;


%Now find the posterior linear term

r1 = mean_prior_reduced'*inverse_variance_prior_reduced;


%Find the new mean_initial based on the Y matrix now

mean_removed = mean_initial(ind_noi_index);
omega_removed =  omega_initial(d, ind_noi_index);

[row1, col1] = size(r1);

r2 = zeros([row1,col1,n]);

for i = 1:n
         
    %use vector multiplication to get the sum
    
    conditional_mean_tmp = mean_initial(d) +(-variance_likelihood*omega_removed)*...
        (Y_initial(i,ind_noi_index)' - mean_removed);
    
    r2(:,:,i) = (Z_two(i,:)*q_vec - conditional_mean_tmp)*(Z_red(i,:) + Z_two(i,:)*W_mat);
 
end


F_vec = (1/variance_likelihood)*sum(r2, 3);

%final posterior linear term

r_vec = r1 - F_vec;  %there has to be a subtraction here










end