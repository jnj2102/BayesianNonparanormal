function [Bayes_est] = estfunction(samples, basis_red, basis_two, W_d_trans, q_d)

%find the functions using each of the iterations of the thetas
%and use same basis functions evaluated at the grid points

%apply it across the iterations so create a function to do this.

     
 Bayes_est = horzcat(basis_red, basis_two)* vertcat(samples, W_d_trans*samples + q_d);
 
end
