function [Bayes_est] = updateY_basis(basis_red,basis_two, W_mat_cell, q_vec_cell,...
                            new_sample)

 %Find the transformed variables
    
   [Bayes_est] = estfunction(new_sample,...  
                   basis_red, basis_two, W_mat_cell, q_vec_cell);
               
             

end