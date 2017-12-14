function [x_matrix, true_f, gridpoint, MLE_norm, MLE_ev, MLE_log, MLE_stable,power] = xfunctions(p, y_true, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generate the x matrix and the true y functions using MLE
%Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the maximum likelihood estimate for each column of Y
%(corresponding to the predictor) using the distribution that corresponds
%to the cdf


MLE_norm =  cell([p,1]);
power = cell([p,1]);
%MLE_t =  cell([p_grp,1]);
MLE_log =  cell([p,1]);
%MLE_skew = cell([p,1]);
MLE_ev = cell([p,1]);
MLE_stable = cell([p,1]);

x_matrix = [];

d = 0;

while (1)
    d = d+1;
    
        if d > p
            break
        else
         
MLE_norm{d} = mle(y_true(:, d));
x_matrix = horzcat(x_matrix, normcdf(y_true(:, d), MLE_norm{d}(1), MLE_norm{d}(2)));

        end
        
%t distribution is no good
%MLE_t{d} = mle(y_true(:, y_grp{d}(2)),'distribution', 'tlocationscale');
  d = d+1;
  
  if d > p
       break
  else
         
MLE_ev{d} = mle(y_true(:, d),'distribution','Extreme Value');
x_matrix = horzcat(x_matrix, cdf('Extreme Value', y_true(:, d), MLE_ev{d}(1), MLE_ev{d}(2)));
  end
  
  d = d+1;
    
  if d > p
      break
  else
MLE_log{d} = mle(y_true(:, d),'distribution', 'logistic');
x_matrix = horzcat(x_matrix, cdf('logistic', y_true(:, d), MLE_log{d}(1), MLE_log{d}(2)));
  end
  
  d = d+1;
  
  if d > p
      break
  else
%the beta is on the boundary - this may or may not be ok for stable.
MLE_stable{d} = mle(y_true(:, d),'distribution', 'stable');
x_matrix = horzcat(x_matrix, cdf('stable', y_true(:, d), MLE_stable{d}(1),...
    MLE_stable{d}(2), MLE_stable{d}(3), MLE_stable{d}(4)));
  end

    d = d+1;
  
  if d > p
      break
  else

power{d} = randi(5);  %select a random integer between 1 and 5
x_matrix = horzcat(x_matrix, (normcdf(y_true(:, d))).^(1/power{d})); %this is mean=0 sigma = 1
  end
    %m is the power for the power function
  
end



gridpoint = zeros(n,p);

for d = 1:p
    column = linspace(quantile(x_matrix(:,d), .025),quantile(x_matrix(:,d), .975), n)';
gridpoint(:,d) = column;

end

true_f = [];

d = 0;

while (1)
    
    d = d+1;
    
    if d > p
            break
        else
         
true_f = horzcat(true_f, norminv(x_matrix(:,d),MLE_norm{d}(1), MLE_norm{d}(2)));

    end
d = d+1;

    if d > p
            break
    else 

true_f = horzcat(true_f, icdf('Extreme Value',x_matrix(:,d), MLE_ev{d}(1), MLE_ev{d}(2)));
    end
    
d = d+1;

    if d > p
            break
        else
         
true_f = horzcat(true_f, icdf('Logistic', x_matrix(:,d),MLE_log{d}(1), MLE_log{d}(2)));
    end
    
d = d+1;

    if d > p
            break
        else
         
true_f = horzcat(true_f, icdf('stable', x_matrix(:,d), MLE_stable{d}(1),...
    MLE_stable{d}(2), MLE_stable{d}(3), MLE_stable{d}(4)));

    end
    
    
d = d+1;

    if d > p
            break
        else
         
true_f = horzcat(true_f, norminv((x_matrix(:,d)).^power{d},0,1));

    end

end






%Old Code

%Let's create the code for 1 predictor, then I'll loop it using horzcat, so
%it's not dependent on a certain number.

%first randomly draw a mean and variance from a certain range.
%In general, you can generate N random numbers in the interval (a,b)
%with the formula r = a + (b-a).*rand(N,1).
% 
% %initialize the x matrix
% x_mat = zeros([n,p]);
% 
% %initialize the mean, var, shape, df with not a number so I know when I
% %actually use them
% final_mean_pred = NaN([p,1]);
% final_stddev_pred = NaN([p,1]);
% final_shape_pred = NaN([p,1]);
% final_df_pred = NaN([p,1]);
% 
% %set the maximum number of draws to 100
% max_draws = 1000;
% 
% %initialize the index for the pdf
% idx_pdf = 0;
% 
% for jdx = 1:p
%  
%     idx_pdf = idx_pdf + 1;
%     
%     for draws = 1:max_draws
%         %Now first check to make sure that we have a predictor that works
% 
%         %I want sigma = standard deviation, not sigma^2 = variance
%         interval_low = mean_true(jdx) -3.5*sqrt(sigma_true(jdx,jdx));
%         interval_high = mean_true(jdx) + 3.5*sqrt(sigma_true(jdx,jdx));
%         
%         y_true_tmp = y_true(:, jdx);
%         indexesInRange = y_true_tmp >= interval_low & y_true_tmp  <= interval_high;
%         
%         subVector_y_true = y_true_tmp(indexesInRange);
%         %find the derivative of cdf (pdf) at the y values using the current draw
%         %(this would cycle through each predictor and cdf) and make sure it is >.1
%         %in the interval_low and interval_high range (I can use the minimum to do
%         %this. Because we want our cdfs to be strictly increasing in this range for
%         %normal, not flat
%         
%         %do if/else statements for the 5 cdfs - I have to list them
%         
%        	%draw a random mean between 0 and 1
%         mean_draw = 0 + (3- (0))*rand(1);
%         
%        %mean_draw = rand(1);
%         
%        %draw a random variance between 2 and 5, so to find the standard
%         %deviation (sigma) then I need to take the square root of this
%         
%         stddev_draw = sqrt(10 + (20-10)*rand(1));
%         
%         %draw a random df between 3 and 100
%         df_draw = 3 + (100-3)*rand(1);
%         
%             
%         %draw a random shape between -.99 and .99
%         shape_draw = -.99 + (.99- (-.99))*rand(1);   
%         
%         if idx_pdf == 1
%                     
%         pdf_tmp =  normpdf(subVector_y_true, mean_draw, stddev_draw); %normal pdf
%              
%         elseif idx_pdf == 2
%             
%         pdf_tmp = tpdf(subVector_y_true, df_draw);
%         
%         elseif idx_pdf == 3
% 
%         pdf_tmp = pdf('Logistic',subVector_y_true, mean_draw, stddev_draw);
%             
%         elseif idx_pdf == 4
%             
%         pdf_tmp = skewtdis_pdf(subVector_y_true, df_draw, shape_draw); %skew t
%         
%         else
%             df_draw = 300; %set this for the skew normal
%             %draw a random shape between -.99 and .99
%             
%             pdf_tmp =  skewtdis_pdf(subVector_y_true,  df_draw, shape_draw); %skew normal (maybe, the largest is 300 =df
%          
%         end
%         
%         %Now find the indexes in the range of my pdf_tmp and then cut of that
%         %pdf_tmp by those indices, and then find the minimum of that.
% 
% %         indexesInRange = pdf_tmp >= interval_low & pdf_tmp <= interval_high;
% %         subVector = pdf_tmp(indexesInRange);
% 
%         %Now find the minimum of my subVector and check to make sure it's >.1
%         min_subVector = min(pdf_tmp);
% 
%         %basically I keep drawing until I get a minimum.  But perhaps I will have a
%         %loop that ends at some max iterations and if I hit it, then I tell myself
%         if min_subVector > .1
%             %save what I want and put it in the vertcat
%             final_mean_pred(jdx) = mean_draw;
%             final_stddev_pred(jdx) = stddev_draw;
%             final_shape_pred(jdx) = shape_draw;
%             final_df_pred(jdx) = df_draw;
%             
%             %I need to do another if/else statement for the cdfs now to
%             %save it for the x matrix, then horzcat it.
%          
%         
%                 if idx_pdf == 1
%             %normal cdf needs sigma = standard deviation, not variance
%             
%                 x_mat_tmp =  normcdf(y_true(:,jdx), mean_draw, stddev_draw); %normal pdf
% 
%                 elseif idx_pdf == 2
% 
%                 x_mat_tmp = tcdf(y_true(:,jdx), df_draw);
% 
%                 elseif idx_pdf == 3
% 
%                 x_mat_tmp = cdf('Logistic',y_true(:,jdx), mean_draw, stddev_draw);
% 
%                 elseif idx_pdf == 4
% 
%                 x_mat_tmp = skewtdis_cdf(y_true(:,jdx), df_draw, shape_draw); %skew t
% 
%                 else
% 
%                 x_mat_tmp =  skewtdis_cdf(y_true(:,jdx),  df_draw, shape_draw); %skew normal (maybe, the largest is 300 =df
% 
%                 end
%             
%                 %give me my x values
%                 x_mat(:, jdx) =  x_mat_tmp;
%             
%             break
%             
%         else
%             continue
%         end
%         
%     end %end the draws
% 
%    if draws == max_draws
%        %display the following  - fprint is a combo of disp and sprintf
%          fprintf('The minimum of the pdf for Predictor %d was never greater than .1', jdx);
%    end
%    
%    %reset the index for the pdf to be zero when it hits five
%    if idx_pdf == 5
%        idx_pdf = 0;
%    end
%  
% 
%    
% end  %end the jdx

%%Old Code%%
% 
% %the way I have this set up, p needs to be divisible by 5
% 
% p_grp = p/5;
% %Function to create any  number of functions for estimation.
% 
% grid_mean = linspace(1,2,p_grp);  %I want to vary the location and scale to be between 0 
%                     %and 3, using the length of p
% grid_var = linspace(3,5,p_grp); 
% grid_shape = linspace(-.99,.99, p_grp); %skew (shape) parameter for skew t.  I am using 1.1 because of 1-lambda in denominator
% %grid_param_2 = 2:1:p_grp+1; 
% %grid_mean_3 = linspace(2,3,p_grp);
% %grid_var_2 = linspace(2,3,p_grp);
% %grid_mean_4 = linspace(3,4,p_grp);
% % %grid_param_5 = 3:2:2*p_grp+1;
% % grid_df = linspace(20,50, p_grp); %degrees of freedom parameter for skew t and normal t
% % 
% % x_1 = cell([p_grp,1]);
% % x_2 = cell([p_grp,1]);
% % x_3 = cell([p_grp,1]);
% % x_4 = cell([p_grp,1]);
% % x_5 = cell([p_grp,1]);
% 


%%These Z_true's were causing me to estimate identity matrices, because
%%that was what I was giving it.

% Z_true = zeros([n,p]);
% 
% for i = 1:n
%    
%     Z_true(i, :) = sqrtm(omega_true)*(y_true(i, :)' - mean_true);
% end

% 
% for d = 1:p_grp
%     
%     
% x_1{d} = normcdf(Z_true(:, y_grp{d}(1)), grid_mean(d), grid_var(d));
% 
% x_2{d} = tcdf(Z_true(:, y_grp{d}(2)),grid_df(d));
% x_3{d} =  cdf('Logistic',Z_true(:, y_grp{d}(3)), grid_mean(d), grid_var(d));
% x_4{d} = skewtdis_cdf(Z_true(:, y_grp{d}(4)), grid_df(d), grid_shape(d)); %skew t
% x_5{d} = skewtdis_cdf(Z_true(:, y_grp{d}(5)),  300, grid_shape(d)); %skew normal (maybe, the largest is 300 =df
% %x_4{d} = pskt(y_true(:, y_grp{d}(4)), grid_mean(d), grid_var(d), grid_shape(d), grid_df(d)); %skew t
% %x_5{d} = pskt(y_true(:, y_grp{d}(5)), grid_mean(d), grid_var(d), grid_shape(d), 1e10); %skew normal
% 
% %x_4{d} = psktJJ(y_true(:, y_grp{d}(4)), grid_mean(d), grid_var(d), grid_shape(d), grid_df(d)); %skew t
% %x_5{d} = psktJJ(y_true(:, y_grp{d}(5)), grid_mean(d), grid_var(d), grid_shape(d)); %skew normal
% 
% %x_4{d} = normcdf(y_true(:, y_grp{d}(4)),grid_mean_4(d),grid_var_2(d)); %this is std dev not variance
% %x_5{d} = tcdf(y_true(:, y_grp{d}(5)), grid_param_5(d));
% end



% 
% for d = 1:p_grp
%     
%     
% x_1{d} = normcdf(y_true(:, y_grp{d}(1)), grid_mean(d), grid_var(d));
% 
% x_2{d} = tcdf(y_true(:, y_grp{d}(2)),grid_df(d));
% x_3{d} =  cdf('Logistic',y_true(:, y_grp{d}(3)), grid_mean(d), grid_var(d));
% x_4{d} = skewtdis_cdf(y_true(:, y_grp{d}(4)), grid_df(d), grid_shape(d)); %skew t
% x_5{d} = skewtdis_cdf(y_true(:, y_grp{d}(5)),  300, grid_shape(d)); %skew normal (maybe, the largest is 300 =df
% %x_4{d} = pskt(y_true(:, y_grp{d}(4)), grid_mean(d), grid_var(d), grid_shape(d), grid_df(d)); %skew t
% %x_5{d} = pskt(y_true(:, y_grp{d}(5)), grid_mean(d), grid_var(d), grid_shape(d), 1e10); %skew normal
% 
% %x_4{d} = psktJJ(y_true(:, y_grp{d}(4)), grid_mean(d), grid_var(d), grid_shape(d), grid_df(d)); %skew t
% %x_5{d} = psktJJ(y_true(:, y_grp{d}(5)), grid_mean(d), grid_var(d), grid_shape(d)); %skew normal
% 
% %x_4{d} = normcdf(y_true(:, y_grp{d}(4)),grid_mean_4(d),grid_var_2(d)); %this is std dev not variance
% %x_5{d} = tcdf(y_true(:, y_grp{d}(5)), grid_param_5(d));
% end

%Here is my X matrix.  Since I am using cdf functions, the

% %x values are already between 0 and 1.
% x_1_mat = cell2mat(x_1');
% x_2_mat = cell2mat(x_2');
% x_3_mat = cell2mat(x_3');
% x_4_mat = cell2mat(x_4');
% x_5_mat = cell2mat(x_5');
% x_matrix = horzcat(x_1_mat,x_2_mat,x_3_mat,x_4_mat,x_5_mat);

%Now Find the true Y variables (true functions) for predictors 1:5

%find the true function f_0(x) for the first predictor
%use the inverse cdf or quantile fn for the true function f_0(x) for d=1:p

%I am looping using 1, 2, 3, 4, 5 in the columns - I could do better than
%this though.
% 
% true_f1 = cell([p_grp,1]);
% true_f2 = cell([p_grp,1]);
% true_f3 = cell([p_grp,1]);
% true_f4 = cell([p_grp,1]);
% true_f5 = cell([p_grp,1]);

% 
% gridpoint_grp_mat = reshape(1:1:p,5,p/5);
% gridpoint_grp = cell([p_grp,1]);
% 
% for d = 1:p_grp
%     gridpoint_grp{d} = gridpoint_grp_mat(:,d);
% end





end