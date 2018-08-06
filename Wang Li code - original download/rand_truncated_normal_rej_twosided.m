function  x = rand_truncated_normal_rej_twosided(mu,sig,below,above)
% Generate one sample from N(mu,sig^2)1_{above > x > below}
% Ref: 
%     Proposition 2.3  C. P. Robert 'Simulation of truncated normal variables'
mu_neg = (below - mu)/sig;
mu_pos = (above - mu)/sig;



if mu_neg>0 & mu_pos>0

    if (mu_pos-mu_neg) > 2*sqrt(exp(1))/(mu_neg+sqrt(mu_neg^2+4))*exp((mu_neg^2-mu_neg*sqrt(mu_neg^2+4))/4) % eq. (2.1)
    
        x = rand_truncated_normal_rej_below(mu,sig,below);
        
        while x > above 
        x = rand_truncated_normal_rej_below(mu,sig,below);
        end
    else
        
        x = rand(1)*(mu_pos-mu_neg)+mu_neg;
        
        while log(rand(1)) > (mu_neg^2 - x^2)/2
        x = rand(1)*(mu_pos-mu_neg)+mu_neg;            
        end
        
        x = x*sig + mu;
                 
    end
    
elseif mu_neg<0 & mu_pos<0   
      
     aa = mu_neg;
     mu_neg = -mu_pos;
     mu_pos = -aa;
     
     mu = -mu;
     
     aa = above;
     above = -below;
     below = -aa;
     
     if (mu_pos-mu_neg) > 2*sqrt(exp(1))/(mu_neg+sqrt(mu_neg^2+4))*exp((mu_neg^2-mu_neg*sqrt(mu_neg^2+4))/4) % eq. (2.1)
    
        x = rand_truncated_normal_rej_below(mu,sig,below);
        
        while x > above 
        x = rand_truncated_normal_rej_below(mu,sig,below);
        end
    else
        
        x = rand(1)*(mu_pos-mu_neg)+mu_neg;
        
        while log(rand(1)) > (mu_neg^2 - x^2)/2
        x = rand(1)*(mu_pos-mu_neg)+mu_neg;            
        end
        
        x = x*sig + mu;
                 
    end
    
    x = -x;
    
else  
    
    x = rand_tn(mu,sig,below,above);

end
