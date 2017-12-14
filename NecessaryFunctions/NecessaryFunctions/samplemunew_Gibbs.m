function [Z,mean_posterior] = samplemunew_Gibbs(Y_matrix, p, n, sigma_initial)
%Function to sample mu


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Find the posterior mean for the transformed observations
%The posterior mean for the transformed observations is N(y_bar, 1/n*Sigma)

sigma_posterior = (1/n) * sigma_initial;

 mean_tmp = mean(Y_matrix)';
 
   mean_posterior =  mvnrnd(mean_tmp, sigma_posterior)'; %these are random draws
   
   %Now subtract these means from the transformed variables to get Z ~ N(0, Omega^{-1}).  
   %These Z's are used to estimate the precision matrix
      
    Z = zeros(n,p);

   for d = 1:p
       column = Y_matrix(:,d) - repmat(mean_posterior(d),[n,1]);
        Z(:,d) = column;
    
   end
   

end