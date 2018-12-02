function [Bayes_est] = estfunction(samples, basis_red, basis_two, W_d_trans, q_d)

%find the functions using the thetas
     
 Bayes_est = horzcat(basis_red, basis_two)* vertcat(samples, W_d_trans*samples + q_d);
 
end
