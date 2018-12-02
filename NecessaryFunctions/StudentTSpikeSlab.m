function[C_omega,W, pii,Sig] = StudentTSpikeSlab(ind_noi_all,sl_variance,Sig,S_std,n,p,W,...
    C_omega,lambda,c0,pii,tau, a,b,f,g,upperind)




for i = 1:p
    
      ind_noi = ind_noi_all(:,i);
 
      sl_temp = sl_variance(ind_noi,i);
      
      Sig11 = Sig(ind_noi,ind_noi);

      Sig12 = Sig(ind_noi,i);

      invC11 = Sig11 - Sig12*Sig12'/Sig(i,i);

      Ci = (S_std(i,i)+lambda)*invC11+diag(1./sl_temp);
      
       Ci = (Ci+Ci')./2;  
       
      % Ci_chol = chol(Ci);   
       
      
     % mu_i = -Ci_chol\(Ci_chol'\S_std(ind_noi,i));
      
      mu_i = -invChol_mex(Ci)*S_std(ind_noi,i);
      
       beta = mu_i+ invChol_mex(Ci)*randn(p-1,1);

     % beta = mu_i+ Ci_chol\randn(p-1,1);
        
        
          %update Omega
        C_omega(ind_noi,i) = beta;
        C_omega(i,ind_noi) = beta;
        
        a_gam = 0.5*n+1;
        b_gam = (S_std(i,i)+lambda)*0.5;
        gam = gamrnd(a_gam,1/b_gam);
      
        
        c_diag = beta'*invC11*beta;
        C_omega(i,i) = gam+c_diag;
      
        invC11beta = invC11*beta;
                
        %update Sigma
        Sig(ind_noi,ind_noi) = invC11+invC11beta*invC11beta'/gam;
        Sig12 = -invC11beta/gam;
        Sig(ind_noi,i) = Sig12;
        Sig(i,ind_noi) = Sig12';
        Sig(i,i) = 1/gam;
        
        
         tau_squared = tau(ind_noi,i);
         
        w1 = -0.5*log(c0*tau_squared) -0.5*beta.^2./(c0*tau_squared)+log(1-pii);
        w2 = -0.5*log(1*tau_squared) -0.5*beta.^2./(1*tau_squared)+log(pii);
        
        
         w_max = max([w1,w2],[],2);

            %Find the Pr(Z=1) 
        prob = exp(w2-w_max)./sum(exp([w1,w2]-repmat(w_max,1,2)),2);

        w = (rand(p-1,1)<prob);  %logical for Z = 1 based on weight w

        %convert w which is a logical array to a numeric array using double
        
        w = double(w);
                
         %change all zeros to c0.
       % w(w==0) = c0;
                
        
        %update W 
        W(ind_noi,i) = w;
        W(i,ind_noi) = w;
                
      
       
       %Lastly tau gets updated
       
       alpha_tau = (f + 1/2)*ones([p-1,1]); %constant vector
       beta_tau = 1./(g + (beta.^2.*(1/2*(w + (1-w)./c0))));
       
       tau_update = 1./(gamrnd(alpha_tau, beta_tau));
         
       %update tau
       tau(ind_noi,i) = tau_update;
       tau(i,ind_noi) = tau_update;

       %update the spike and slab variances for omega
       
       sl_variance(ind_noi,i) = w.*tau_update;
       sl_variance(i,ind_noi) = w.*tau_update;
            
       
end       
        
   
       %update the pii
       %Have to use logicals because updating the whole column vector
        %if W_vector = 1 then that's a 1 and that makes W_vector = 0
        %correspond to c0 automatically
       pii = betarnd(a + sum(W(upperind) == 1), b + sum(W(upperind) == 0));
   

end