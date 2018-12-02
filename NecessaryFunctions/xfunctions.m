function [x_matrix, MLE_norm, MLE_ev, MLE_log, MLE_stable,power] = xfunctions(p, y_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generate the x matrix  using MLE
%Author: Jami Jackson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the maximum likelihood estimate for each column of Y
%(corresponding to the predictor) using the distribution that corresponds
%to the cdf


MLE_norm =  cell([p,1]);
power = cell([p,1]);
MLE_log =  cell([p,1]);
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








end