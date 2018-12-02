function [BIC]= FindBIC(Omega_bayes_est,final_edge_matrix,mean_Z_Bayes_est,n,p)

%Function to compute the BIC, EBIC, and EBIC with 1/2 tuning parameter
%Author: Jami Jackson Mulgrave


%find nonzero elements of the diagonal of the upper triangular part of omega 
Omega_bayes_est_triu = triu(Omega_bayes_est);

%we are including all of the diagonal entries.
[diag_index, ~] = find(diag(Omega_bayes_est_triu) ~= 0);
%find the nonzero elements that relate to the offdiagonal nonzero elements of the strict upper
%triangular part of omega.  It is done this way because the edge matrix has
%all zeros down the diagonal elements.
final_edge_matrix_strict_triu = triu(final_edge_matrix, 1);

[row_index, ~] = find(final_edge_matrix_strict_triu(:) ~=0);

%convert the index to subscripts.  You have to tell it the dimensions of
%the matrix.
[row_omega, col_omega] = ind2sub(size(final_edge_matrix),row_index);

number_nonzero_E_1 = numel([diag_index;row_omega]); %number of nonzero elements = q
number_nonzero_E_2 = numel([diag_index;col_omega]); 

nonzero_E_1 = [diag_index;row_omega];
nonzero_E_2 = [diag_index;col_omega];

%create the E_1 and E_2 matrices, which both are p x q matrices, where the
%p is the index from both the row and the column

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


%find the S = ZZ'
S_matrix = mean_Z_Bayes_est'*mean_Z_Bayes_est; 

%Write an anonymous function that calculates the objective.

fun = @(x)-n*log(det((E_1*diag(x)*E_2' + E_2*diag(x)*E_1')))  +...
    trace((E_1*diag(x)*E_2' + E_2*diag(x)*E_1')*S_matrix);


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

%chol is not used in place of log(det(*)) because Omega_MLE may not be
%positive definite since it is sparse
BIC = -2*n*(log(det(Omega_MLE))) +...
    2*trace(Omega_MLE*S_matrix) +...
    num_upperdiag_omega_MLE*log(n);

end