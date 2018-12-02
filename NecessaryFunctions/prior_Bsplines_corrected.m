function [F_mat, g_vec, W_mat, q_vec, Fconstraint, A_mat, knot_vector, index_1, index_2,inverse_variance_prior_reduced,...
    mean_prior_reduced, LHS_gq, RHS_gq, Z_red, Z_two] = prior_Bsplines_corrected(optimalJ, minK, mu, tau, sigma2, N_points, c, x_matrix, d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functions builds the prior for the B-splines for the Bayesian 
%nonparanormal method.  It finds the constraint matrix, mean of the prior
% and inverse variance matrix of the prior.
%Author: Jami Jackson
%
%Input:  J: the optimal number of B-splines, found from AIC calculations
%        mu: User specified hyperparameter for the prior mean
%       tau: User specified hyperparameter for the prior mean
%       sigma2: User specified hyperparameter for the prior variance
%
%Output: F: Monotonicity constraint matrix
%        g: Monotonicity constraint vector
%       A_ls: Matrix used for the initial values of theta
%        b_ls: Vector used for the initial values of theta
%       W_mat: Matrix used to determine the two thetas
%       q_vec: Vector used to determine the two thetas
%       Index_1: Index of the first theta
%       Index_2: Index of the second theta
%       inverse_variance_prior_reduced: inverse variance of the prior for
%         the (J-2) thetas
%       mean_prior_reduced: mean of the prior for the (J-2) thetas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

%Step 1
%First put a prior on the thetas, distributed Multivariate normal (MVN)

K = minK(d);
J = optimalJ(d);
 j = 1:J;
 prob = (j - .375)/(J-.75+1); %vector of probabilities for the quantile fn.

%These are the initial prior mean and variance
mu_star = mu + tau*norminv(prob,0,1); 
sigma2_star = sigma2*eye(J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Step 3
%Create the b-spline matrix for the linear constraints

%The basis function is to be evaluated at .25,.5, and .75 for the 
%linear constraints.  

x_lincon = [.25,.5,.75];

knot_vector = [zeros(4,1).',1/K:1/K:(K-1)/K,ones(4,1).']; 

A_bas = bspline_basismatrix(4,knot_vector,x_lincon);

%Now I will create a matrix for the constraints

A_mat = [A_bas(2,:);A_bas(3,:)-A_bas(1,:)];

%here are the final prior mean and variance-covariance matrices 
%incorporating the linear constraints

%it is better to do the forward slash (or the backward slash) to calculate
%the inverse than to find inv directly.

mean_prior = mu_star.' + transpose(A_mat)/(A_mat * transpose(A_mat))*(c-A_mat*mu_star.');

variance_prior = sigma2_star*(eye(J) - transpose(A_mat)/(A_mat * transpose(A_mat))*A_mat);


%Determine which columns have a nonzero element for each row
[~, col_row1] = find(A_mat(1,:)); %row 1

[~, col_row2] = find(A_mat(2,:)); %row 2

%make the column numbers unique so I can subset with it
%unique_col = unique(col_row1, col_row2);


%Now I am picking the first and second index.  I should be able to pick any
%of these columns
index_1 = col_row1(1);  
index_2_tmp = col_row2(col_row2 ~= index_1); %remove the index for the first
                                            %theta                                 
index_2 = index_2_tmp(1);
 
% I remove the index_1 and index_2 rows and columns of the prior variance matrix
%to make it nonsingular

%I remove the index_1 and index_2 rows of the prior mean vector to match

variance_prior_reduced = variance_prior;

variance_prior_reduced(:, [index_1,index_2]) = [];
variance_prior_reduced([index_1,index_2],:) = [];

inverse_variance_prior_reduced = inv(variance_prior_reduced);

mean_prior_reduced = mean_prior;
mean_prior_reduced([index_1,index_2]) = [];


%Now solve for the two thetas that can be determined from the 2 linear
%constraints using the symbolic toolbox in MATLAB

%create the vector of symbolic variables
sym_vec = sym('a', [J 1]);

%create the symbolic functions using a loop

eqn1_tmp = [];

%makes the symbolic variables into a vector
for j = 1:J
   eqn1_tmp = horzcat(eqn1_tmp, A_mat(1,j)*sym_vec(j)); 
    
end

%sum the symbolic variables in the vector to get my first equation to solve
%for, the first linear constraint

eqn1 = sum(eqn1_tmp) == 0;

%create the second equation to solve for
eqn2_tmp = [];

%makes the symbolic variables into a vector
for j = 1:J
   eqn2_tmp = horzcat(eqn2_tmp, A_mat(2,j)*sym_vec(j)); 
    
end

%sum the symbolic variables in the vector to get my second equation to solve
%for, the second linear constraint

eqn2 = sum(eqn2_tmp) == 1;

%Decide which two thetas to solve for using the index_1 and index_2 and
%applying it to the sym_vec

theta1 = sym_vec(index_1);
theta2 = sym_vec(index_2);

%Now solve for those two thetas using the two linear constraints

[sol1, sol2] = solve([eqn1, eqn2], [theta1, theta2]);

%Now find the monotonicity constraint matrix F and its vector g_vec and the
%matrix to find the two thetas, W_mat and its vector q_vec 

%Find the reduced sym_vec and the sym_vec of the 2 thetas

sym_vec_red = sym_vec;

sym_vec_red([index_1, index_2]) = [];

%sym_vec_two = sym_vec([index_1,index_2]);

%find the coefficients and terms of the first and second solution

[coeff_sol1, terms_sol1] = coeffs(sol1);

[coeff_sol2, terms_sol2] = coeffs(sol2);

%The W matrix only has two rows because we are only solving for two thetas
%because there are only two linear constraints

W_mat= zeros([2, J-2]);

%find the first row of the W matrix

%I am not using a loop for the rows of the W matrix because the solution 1
%and solution 2 are already separated
indx1 = ismember(sym_vec_red, terms_sol1); %index of the first row of W matrix

%find the constant term and save it for the constant vector and remove it
%to find the matching coefficients of the reduced theta vector

coeff_sol1_red = coeff_sol1;
coeff_sol1_red(terms_sol1 == 1) = [];
        
W_mat(1,indx1) = coeff_sol1_red;


%Find the second row of the W matrix

indx2 = ismember(sym_vec_red, terms_sol2); %index of the first row of W matrix

%find the constant term and save it for the constant vector and remove it
%to find the matching coefficients of the reduced theta vector

coeff_sol2_red = coeff_sol2;
coeff_sol2_red(terms_sol2 == 1) = [];
W_mat(2, indx2) = coeff_sol2_red;

%Create the constant vector q_vec

q1 = coeff_sol1(terms_sol1 == 1); %the constant term is represented as 1 
             %in Matlab so I am indexing the logical to get the coefficient
             %of the constant term.  I can have an empty vector like
             %this because I am not indexing q1 because I am not looping
             
%there may not be constant terms so I am changing empty vectors due to no
%constant terms to zero

if isempty(q1) == 1
    q1 = 0;
end

q2 = coeff_sol2(terms_sol2 == 1); 

if isempty(q2) == 1
    q2 = 0;
end
        
q_vec = double([q1; q2]);  %I need to change this class of symbolic to double


%now simplify the expressions to get each row of the monotonicity
%constraint matrix

sym_vec_subs = sym_vec;

%substitute the symbolic vector with the solutions for the two thetas
sym_vec_subs(index_1) = sol1;
sym_vec_subs(index_2) = sol2;

%I need an inner loop since I don't know ahead of time the length of the
%terms.  And I need another loop for each row

F_mat = zeros([J-1, J-2]);
g_vec = zeros([J-1,1]);
    
for j = 1:J-1 %loop through the simplified expressions for the F matrix
        F_tmp = simplify(sym_vec_subs(j+1)-sym_vec_subs(j));

        [F_coeffs_tmp, F_terms_tmp] = coeffs(F_tmp);

        index = ismember(sym_vec_red, F_terms_tmp);
    
        %put the respective coefficients in the columns of the F matrix using
         %the index
    
         F_coeffs_tmp_red = F_coeffs_tmp;
         F_coeffs_tmp_red(F_terms_tmp == 1) = [];
         F_mat(j, index) = F_coeffs_tmp_red;

        %Create the constant vector g_vec
        
        if isempty(F_coeffs_tmp(F_terms_tmp == 1)) == 1
            g_vec(j) = 0 ; %the constant term is represented as 1 
             %in Matlab so I am indexing the logical to get the coefficient
             %of the constant term
             
              %there may not be constant terms so I am changing empty vectors due to no
              %constant terms to zero
        else
            g_vec(j) =  F_coeffs_tmp(F_terms_tmp == 1);
        end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 3

%Now find the initial values by estimating the integral of the 
    %function from -5 to 5

%this is the LHS of the equation

% 
% v = linspace(-5,5); 
% 
% basis_normcdf =  bspline_basismatrix(4,knot_vector,normcdf(v));
% 
% 
% %loop through the k's
% 
% %this is not matrix multiplication, it's element-wise multiplication
% 
% fv{J} = [];
% estimate = zeros([J,1]);
% 
% for k=1:J
%     fv{k} = transpose(v).*basis_normcdf(:,k).*transpose(normpdf(v,0,1));
%     estimate(k) = simps(v,fv{k});
% end
% 
% 
% %this is the RHS of the equation
% 
% % %loop through j=1:J
% % %create a sequence from 0 to 1
% x_seq = linspace(0.001,.99);
% fns{J}{J} = [];
% estimate_matrix = zeros(J,J);
% j_tmp = bspline_basismatrix(4,knot_vector,x_seq);
%         
% for k=1:J
%     for j = 1:J
%         fns{k}{j} = j_tmp(:,j).*j_tmp(:,k);
%         estimate_matrix(k,j) = simps(x_seq, fns{k}{j});
%     end
% end
% 
% 
% %remove the 2 rows of the estimate_matrix (index_1 and index_2)
% 
% estimate_matrix_red = estimate_matrix;
% estimate_matrix_red(:, [index_1,index_2]) = [];
% %estimate_matrix_red = estimate_matrix(:,[1:3,6:8]);
% %estimate_matrix_two = estimate_matrix(:,4:5);
% 
% estimate_matrix_two = estimate_matrix(:, [index_1,index_2]);
% 
% 
% b_ls = double(estimate - estimate_matrix_two*q_vec);
% A_ls = double(estimate_matrix_red + estimate_matrix_two*W_mat);
% 
% 
% %Create the full basis matrix for the dth predictor
Z_mat = bspline_basismatrix(4,knot_vector,x_matrix(:,d));

%Reduce the basis matrix by the index of the two thetas that we solved for

Z_red = Z_mat;
Z_red(:,[index_1,index_2]) = []; 

%Find the basis matrix for the two thetas that we solved for

Z_two = Z_mat(:,[index_1,index_2]);
% 
% U_mean = repmat(mean(Z_red + Z_two*W_mat),[n,1]);
% 
% U_mat = (Z_red + Z_two*W_mat) - U_mean;




diag_tmp = diag(ones([J-2,1]), 1);  %there is one less 1 once it's in the 
%off diagonal

Fconstraint_tmp = diag(repmat(-1,[J-1,1])) + diag_tmp;

%need to add one more column with the last 1
Fconstraint = horzcat(Fconstraint_tmp, [zeros([1,J-2]), 1]' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try to use Gauss Quadrature to generate the initial values

%Try 20 points over interval [-5,5] to cover the range for N(0,1)
% 
% N_points = 20;
% a_interval = -5;
% b_interval = 5;
% 
% [nodes, weights] = lgwt(N_points,a_interval,b_interval);
% 
% %nodes weighted by interval [a,b]
% nodes_LHS = nodes*(b_interval-a_interval)/2 + (b_interval - a_interval)/2;
% 
% %basis function cdf
% bspline_matrix_cdf =  bspline_basismatrix(4,knot_vector,normcdf(nodes_LHS));
% 
% %find the LHS function
% fn_LHS = bspline_matrix_cdf .* nodes_LHS;
% 
% %weight the LHS function
% fn_LHS_weighted = weights .* fn_LHS;
% 
% %sum the weighted LHS function.  Sum the m points (the 20 points), leaving
% %J
% 
% sum_fn_LHS_weighted = sum(fn_LHS_weighted);
% 
% %final LHS Gauss -Quad approx
% 
% LHS_gq = transpose((b_interval - a_interval)/2 * sum_fn_LHS_weighted);
% 
% %Find the RHS
% 
% fn_RHS_weighted = zeros(J,J);
% fn_RHS_weighted_cell = cell([N_points,1]);
% 
% for pts = 1: N_points
%     for j = 1:J
%     
%          for k = 1:J
%                 %Make the function for each N.  Then sum the N
% fn_RHS_weighted(j,k) = weights(pts) * bspline_matrix_cdf(pts,j) .* bspline_matrix_cdf(pts,k);
% 
%           end
%     
%     end
%     
%     %multiply each of the functions by the weight
%     fn_RHS_weighted_cell{pts} = fn_RHS_weighted;
% end
% 
% %sum the weighted functions
% 
% sum_fn_RHS_weighted = sum( cat(3,fn_RHS_weighted_cell{:}),3);
% 
% %Ok now find the final Right hand side (RHS)
% 
% RHS_gq = (b_interval - a_interval)/2 * sum_fn_RHS_weighted;



%%% Use Gauss Hermite for the initial values

%N_points = 20;

[nodes, weights] = GaussHermite_2(N_points);

   
%basis function cdf
bspline_matrix_cdf =  bspline_basismatrix(4,knot_vector,normcdf(nodes));

%find the LHS function
fn_LHS = bspline_matrix_cdf .* nodes;

%weight the LHS function
fn_LHS_weighted = weights .* fn_LHS;

%sum the weighted LHS function.  Sum the m points (the 20 points), leaving
%J

%final LHS Gauss -Quad approx
LHS_gq = sum(fn_LHS_weighted);

%Find the RHS

fn_RHS_weighted = zeros(J,J);
fn_RHS_weighted_cell = cell([N_points,1]);


for pts = 1: N_points
    for j = 1:J
    
         for k = 1:J
                %Make the function for each N.  Then sum the N
fn_RHS_weighted(j,k) = weights(pts) * bspline_matrix_cdf(pts,j) .* bspline_matrix_cdf(pts,k);

          end
    
    end
    
    %multiply each of the functions by the weight
    fn_RHS_weighted_cell{pts} = fn_RHS_weighted;
end

%sum the weighted functions
%Ok now find the final Right hand side (RHS)

RHS_gq = sum( cat(3,fn_RHS_weighted_cell{:}),3);
     
 
    










end



