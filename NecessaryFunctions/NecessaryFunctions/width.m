function [width] = width(X, alpha)
%Function to calculate the width of a (1-alpha)% interval


 width = quantile(X, (1-alpha/2)) - quantile(X, alpha/2);
    
end