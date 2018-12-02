function [optimalJ, minK, Final_AIC] = AICloop_updated(n,p, x_matrix,  mu, tau, sigma2,c)
%
%This is a function that calculates the AIC for multiple choices of J.
%It starts with J=4 and calculates the AIC for J=4 and for each single
%increase of J.  It monitors a decrease in AIC and stops when it has 
%identified an increase in the AIC.  It takes the J that is at the lowest
%AIC.

%Input:
%sample size: n
%number of predictors: p
%iters: the current iteration number
%sparsity_string: the sparsity that is being used in string form
%
%Output:
%OptimalJ: the optimal J selected
%minK: the optimal K selected
%Final_AIC: the AIC values that were found
%
%Author: Jami Jackson Mulgrave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%specify the order of the B-spline 
Q = 4;
 
%specify the number of knot intervals with K. The J are the number of
%Basis functions.  

 %The J increments from K=1 until it reaches a K that makes the AIC
 %increase after having decreased.

 
%initialize the loop with the first AIC calculation for J=4 (K=1)

 %Horzcat each of these AIC_pred for each J.


minK = zeros([p,1]);

optimalJ = zeros([p,1]);
    


Final_AIC = cell([p,1]);


for d = 1:p
    
    AIC_matrix = [];
    
    K = 1; % for each predictor, 
            %we will start with K=1 and incrementally add to the K in the loop.

    AIC_matrix = AICcalculation_corrected(n, x_matrix, K, Q, d,  mu, tau, sigma2,c); %initial calculation for 
                                                %each predictor
                                                
    low_num = AIC_matrix;
    
   while( K <= 97) %stops when K > 97

            test_matrix = [];  %test the next 10 values to 
    
            for l = K+1 : K+10
        
                AIC_pred =  AICcalculation_corrected(n, x_matrix, l, Q, d,  mu, tau, sigma2,c);

                AIC_matrix = vertcat(AIC_matrix, AIC_pred);
                
                test_matrix = vertcat(test_matrix, AIC_pred);
                 
            end
            %the point of this is to find the lowest number - minimize the
            %AIC.
            if all(low_num < test_matrix) == 1
                %if the low number is smaller than all the next 10 values,
                %then stop
                              
                break
                
            else
                %find the position of the minimum K and 
                %set the new K as the l value since the smallest K was found in
                %the set of 10 values.
                [minKvalue, minKposition] = min(test_matrix);
                low_num = minKvalue; %set the new low number
                K = minKposition + K;
                continue
                
            end                  
        
   end
  
    %find the minimum J and minimum k

        Final_AIC{d} = AIC_matrix;
        [~, finalminKvalue] = min(AIC_matrix); %all the positions of the minimum Ks.
        minK(d,:) = finalminKvalue;
       optimalJ(d,:) = minK(d,:) + 3;
end









end