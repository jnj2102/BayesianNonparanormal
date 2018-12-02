function[S_std] = SumProductsMatrix(Z,p,n)

%create the sum of products matrix
S = Z'*Z;

S_ExpCorrC = zeros(p,p);
  
for j = 1:p
    for k = 1:p
        S_ExpCorrC(j,k) = S(j,k)/(sqrt(S(j,j))*sqrt(S(k,k)));
    end
end


 S_std = S_ExpCorrC*n; %standardize the S
 
 
 
 end