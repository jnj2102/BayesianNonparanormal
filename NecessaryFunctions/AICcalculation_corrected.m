function [AIC_pred] = AICcalculation_corrected(n, x_matrix, K, Q, d, mu, tau, sigma2,c)
%Function that calculates the AIC for the dth predictor

J = K + Q -1;

j = 1:J; %this actually gets used for prob calculation


%generate the initial prior means of theta using the number of knots

 prob = (j - .375)/(J-.75+1); %vector of probabilities for the quantile fn.

 

%These are the initial prior mean and variance
mu_star = mu + tau*norminv(prob,0,1); %Matlab doesn't like dots
sigma2_star = sigma2*eye(J);


x.lincon = [.25,.5,.75];

knot_vector = [zeros(4,1).',1/K:1/K:(K-1)/K,ones(4,1).']; 

A_bas = bspline_basismatrix(4,knot_vector,x.lincon);


%Now I will create a matrix for the constraints.  This is used in the
%quadprog calculation

A_mat = [A_bas(2,:);A_bas(3,:)-A_bas(1,:)];


%Monotonicity constraint that is not reduced by the two thetas that we can
%solve

diag_tmp = diag(ones([J-2,1]), 1);  %there is one less 1 once it's in the 
%off diagonal

Fconstraint_tmp = diag(repmat(-1,[J-1,1])) + diag_tmp;

%need to add one more column with the last 1
Fconstraint = horzcat(Fconstraint_tmp, [zeros([1,J-2]), 1]' );


gconstraint = repmat(1e-4,[J-1,1]);  %optimization doesn't like strict inequalities 

%find the basis function evaluated using the original x values

    basis_full_matrix = bspline_basismatrix(4, knot_vector, x_matrix(:,d));

%Minimize the thetas wrt the monotonicity constraint and the 2 linear
%constraints


%Find the mean of the basis matrix for each predictor p (each row is the
%predictor
    
basis_mean =  mean(basis_full_matrix);

Z_mat =  basis_full_matrix - repmat(basis_mean, [n,1]);

%Now try the minimization using CVX

%Yay!  This works!!

Z_tt = transpose(Z_mat)*Z_mat;

%Now minimize the function with respect to the thetas

theta = quadprog(Z_tt, zeros([1,J]), -Fconstraint, gconstraint, A_mat,c);

%Find the AIC for the dth predictor

AIC_pred = n * log( transpose(theta)*Z_tt*theta ) + 2*J ; 



end