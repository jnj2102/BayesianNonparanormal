function [optimalJ, minK, Final_AIC] = AICloop(n,p, x_matrix,  mu, tau, sigma2)
%
%This is a function that calculates the AIC for multiple choices of J.
%It starts with J=5 and calculates the AIC for J=5 and for each single
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
%Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%specify the order of the B-spline 
Q = 4;
 
%specify the number of knot intervals with K. The J are the number of
%Basis functions.  

 %The J increments from K=2 until it reaches a K that makes the AIC
 %increase after having decreased.

 
%initialize the loop with the first AIC calculation for J=5 (K=2)

 %I could horzcat each of these AIC_pred for each J.
%Now do K=3 and beyond in a loop.  I basically need to initialize with K=2
%because I'm not stopping with just the first K=2 and I need to compare to
%it.

    


Final_AIC = cell([p,1]);

tic
for d = 1:p
    
    AIC_matrix = [];
    K = 2; %we will start with K=2 and incrementally add to the K in the loop.

    AIC_matrix = AICcalculation(n, x_matrix, K, Q, d,  mu, tau, sigma2); %initial calculation for 
                                                %each predictor
                                                
    low_num = AIC_matrix;
    
   while( K <= 17)

    K = K+1;

    AIC_pred =  AICcalculation(n, x_matrix, K, Q, d, mu, tau, sigma2);
    
    AIC_matrix = vertcat(AIC_matrix, AIC_pred);


        if AIC_pred < low_num
            
            low_num = AIC_pred;
            continue
        else  %do a for loop to run a couple more AICs to make sure it's trending up
            
            test_matrix = [];
    
            for l = K+1 : K+5
        
                AIC_pred =  AICcalculation(n, x_matrix, l, Q, d,  mu, tau, sigma2);

                AIC_matrix = vertcat(AIC_matrix, AIC_pred);
                
                test_matrix = vertcat(test_matrix, AIC_pred);
                 
            end
            %the point of this is to find the lowest number - minimize the
            %AIC.
            if all(low_num < test_matrix) == 1
                              
                break
                
            else
                             
                K = l;
                continue
                
            end                  
        end
   end
    %Keep doing the loop for a maximum of 20 J.  If it stops at 20 = J,
    %then that's the number of basis functions that I am using.  That's the
    %maximum smoothness that we want.
    
    %My code is written this way to avoid searching all the up to 20 if the
    %AIC values are trending up after 5 values.  It stops earlier than 20 if
    %the AIC values are trending up, to save time.
        Final_AIC{d} = AIC_matrix;

end

toc

%Find the minimum number in each of the final AICs and find the position of
%that minimum number - that position is how we get K. K = position + 1.
%Because I started with K=2 not K=1.  Or Save that K in another vector for
%later.

%plot the AICs

minK = zeros([p,1]);

optimalJ = zeros([p,1]);

for d = 1:p
    
    [~, iK] = min(Final_AIC{d});
    
    minK(d,:) = iK+1;
        
    optimalJ(d,:) = minK(d,:) + 3;
    
    size_tmp = size(Final_AIC{d});
    
   x_seq_end = size_tmp(1);
   
   x_plot = (2+3) : (x_seq_end + 4); %I am using k=2 to begin and the index 
                                       %at the end needs 1 more added to it
                                        % plus the 3 for the J total, which
                                        % is why four is being added
   
 %  %H =  plot(x_plot', Final_AIC{d});
   plot(x_plot', Final_AIC{d});
   title(sprintf('Plot of Number of Basis Splines and AIC for Predictor %d',d));
   xlabel('Number of Basis Splines');
   ylabel('AIC');
    
   % %s is for strings %d is for integers signed
   
   %I don't need to plot these AICs anymore.  I can always plot them later
   %since I save the Final AICs for each iteration.
   
  % %  saveas(H, sprintf('Plot of AIC for Predictor %d for Replication %d n=%d p=%d %s.jpg', d, iters, n, p, sparsity_string))
end









end