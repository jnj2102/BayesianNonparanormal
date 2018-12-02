function [star_model, AR2_model, AR4_model, circle_model, full_model, AR1_model, band_model,...
    fivepercent_model,tenpercent_model, twentypercent_model, twopercent_model,...
    fifteenpercent_model] = model_simulation(p)


%Function to create the models.  The input is the dimension.
%the output are the models.

%I could make it easier and remove the loops like Hao Wang

%star model
star_model = zeros(p,p);

for j=1:p
    for k = 1:p
        if j == k
            star_model(j,k) = 1;
        elseif j == 1
            star_model(j,k) = 0.1;
        elseif k == 1
            star_model(j,k) = 0.1;
        end
    end
end
            
%AR(2) model
AR2_model = zeros(p,p);


for j=1:p
    for k = 1:p
        if j == k
            AR2_model(j,k) = 1;
        elseif k == j-1
            AR2_model(j,k) = 0.5;
        elseif j == k-1
            AR2_model(j,k) = 0.5;
        elseif k == j-2
            AR2_model(j,k) = 0.25;
        elseif j == k-2
            AR2_model(j,k) = 0.25;
        end
    end
end


%AR(4) Model

AR4_model = zeros(p,p);

for j=1:p
    for k = 1:p
        if j == k
            AR4_model(j,k) = 1;
        elseif k == j-1
            AR4_model(j,k) = 0.2;
        elseif j == k-1
            AR4_model(j,k) = 0.2;
        elseif k == j-2
            AR4_model(j,k) = 0.2;
        elseif j == k-2
            AR4_model(j,k) = 0.2;
        elseif k == j-3
            AR4_model(j,k) = 0.2;
        elseif j == k-3
            AR4_model(j,k) = 0.2;
        elseif k == j-4
            AR4_model(j,k) = 0.1;
        elseif j == k-4
            AR4_model(j,k) = 0.1;
        end
    end
end


%circle model
circle_model = zeros(p,p);

for j=1:p
    for k = 1:p
        if j == k
            circle_model(j,k) = 2;
        elseif k == j-1
            circle_model(j,k) = 1;
        elseif j == k-1
            circle_model(j,k) = 1;
        elseif j == 1 && k == p
            circle_model(j,k) = 0.9;
        elseif j == p && k == 1 
            circle_model(j,k) = 0.9;
        end
    end
end

%full model
full_model = zeros(p,p);

for j = 1:p
    for k = 1:p
        if j == k
            full_model(j,k) = 2;
        else
            full_model(j,k) = 1;
        end
    end
end



%%%% AR(1) case
%SigTrue = toeplitz(0.7.^[0:p-1]);
% 
%Taking this inverse is problematic because of numerical erros
% SigTrue = zeros(p,p);
% 
% for j = 1:p
%     for k = 1:p
%         exponent = abs(j-k);
%         SigTrue(j,k) = 0.7^exponent;
%     end
% end
% 
% AR1_model = inv(SigTrue);

AR1_model = diag(repmat(2.9216, [p,1]));

for j=1:p
    for k = 1:p
        if j == 1 && k == 1
            AR1_model(j,k) =  1.9608;
        elseif j == p && k == p
            AR1_model(j,k) =  1.9608;
        elseif k == j-1
            AR1_model(j,k) = -1.3725;
        elseif j == k-1
            AR1_model(j,k) = -1.3725;
            
        end
    end
end



ind_noi_all = zeros(p-1,p);

for i = 1:p
       if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
       end
       ind_noi_all(:,i) = ind_noi;
       
end


%block model
% 
% block_model = eye(p);
% 
% for j=1:p/2
%     ind_noi = ind_noi_all(:,j);
%     for k = 1:(p/2-1)
%         index = ind_noi(k);
%         block_model(j,index) = 0.7;
%     end
% end
% 
% for j = (p/2+1):p
%     ind_noi = ind_noi_all(:,j);
%     for k = (p/2)+1:(p-1)
%         index = ind_noi(k);
%         block_model(j,index) = 0.7;
%     end
% end

%block_model = inv(Sigma_True);

%Sparse model
% 
% indmx = reshape(1:p^2,p,p); 
% upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
% lowerind = indmx(tril(indmx,-1)>0); 
%  
% B_matrix = zeros(p);
% %the problem is I'm going from matrix to vector
% 
% for elem = 1: length(upperind)
%     uniform = rand;
%     ind = upperind(elem);
%     if uniform <= 0.1
%         B_matrix(ind) = uniform;
%     end
% end
% 
% B_matrix = B_matrix + triu(B_matrix,1)';  %make into a symmetrix matrix
% 
% Omega = B_matrix + .2*eye(p);
% condition_number = cond(Omega); %the smaller (but greater than 1) the better

%%%Band Graph

band_model = eye(p);


for j=1:p
    for k = 1:p
        if k == j+1
            band_model(j,k) = 0.6;
        elseif j == k+1
            band_model(j,k) = 0.6;
        elseif k == j+2
            band_model(j,k) = 0.3;
        elseif j == k+2
            band_model(j,k) = 0.3;
        end
    end
end

%Percent sparsity


indmx = reshape(1:p^2,p,p); 
 upperind_diag = indmx(triu(indmx)>0);  %include the diagonal
 upperind = indmx(triu(indmx,1)>0);  %do not include the diagonal
 lowerind = indmx(tril(indmx,-1)>0);

fivepercent = zeros(p);

for i = 1:length(lowerind)
    temp = binornd(1, .05);  %5% nonzero
    index = lowerind(i);
    if temp == 1
        fivepercent(index) = normrnd(0,1);
    end
end

diag_temp = normrnd(1, .1, [p,1]);

fivepercent(logical(eye(p))) = diag_temp;

%now make it a precision matrix by LL' = (U')(U)

fivepercent_model = fivepercent*fivepercent'; %this should work



%10 percent
tenpercent = zeros(p);

for i = 1:length(lowerind)
    temp = binornd(1, .10);  % 10% nonzero
    index = lowerind(i);
    if temp == 1
        tenpercent(index) = normrnd(0,1);
    end
end


diag_temp = normrnd(1, .1, [p,1]);

tenpercent(logical(eye(p))) = diag_temp;

%now make it a precision matrix by LL' = (U')(U)

tenpercent_model = tenpercent*tenpercent'; %this should work



%20 percent
twentypercent = zeros(p);

for i = 1:length(lowerind)
    temp = binornd(1, .20);  % 20% nonzero
    index = lowerind(i);
    if temp == 1
        twentypercent(index) = normrnd(0,1);
    end
end


diag_temp = normrnd(1, .1, [p,1]);

twentypercent(logical(eye(p))) = diag_temp;

%now make it a precision matrix by LL' = (U')(U)

twentypercent_model = twentypercent*twentypercent'; %this should work



%2 percent
twopercent = zeros(p);

for i = 1:length(lowerind)
    temp = binornd(1, .02);  % 2% nonzero
    index = lowerind(i);
    if temp == 1
        twopercent(index) = normrnd(0,1);
    end
end


diag_temp = normrnd(1, .1, [p,1]);

twopercent(logical(eye(p))) = diag_temp;

%now make it a precision matrix by LL' = (U')(U)

twopercent_model = twopercent*twopercent'; %this should work



%15 percent
fifteenpercent = zeros(p);

for i = 1:length(lowerind)
    temp = binornd(1, .15);  % 15% nonzero
    index = lowerind(i);
    if temp == 1
        fifteenpercent(index) = normrnd(0,1);
    end
end


diag_temp = normrnd(1, .1, [p,1]);

fifteenpercent(logical(eye(p))) = diag_temp;

%now make it a precision matrix by LL' = (U')(U)

fifteenpercent_model = fifteenpercent*fifteenpercent'; %this should work



end
%% end of function

