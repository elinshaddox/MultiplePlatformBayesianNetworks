function  x = rand_truncated_normal_rej_below(mu,sig,below)
% Generate one sample from N(mu,sig^2)1_{x>below}
% Ref: 
%     Proposition 2.3  C. P. Robert 'Simulation of truncated normal variables'
mu_neg = (below - mu)/sig;

if mu_neg < 0
    
    x = randn(1);
    while x < mu_neg
    x = randn(1);
    end    
    x = x * sig+mu;
else
    
    alpha = (mu_neg + sqrt(mu_neg^2+4))/2; % Optimal alpha
    
    x = exprnd(1/alpha) + mu_neg;
        
    while log(rand(1)) > (-(x-alpha)^2/2) 
         x = exprnd(1/alpha) + mu_neg;
    end
    
    x = x*sig + mu;
    
end






