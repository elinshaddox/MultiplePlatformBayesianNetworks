function [C_save, Sig_save, adj_save, Phi_save, Theta_save, ar_gamma, ar_theta, ...
    ar_zeta, ar_phi, nu_save, ar_nu, w_save, ar_w] = MCMC_multiple_graphs_SSVS_platforms(Theta, Phi, Sig, V0, V1, lambda, pii, n, SampCov, ...
    C, nu, alpha, beta, a, b, my_w, d, f, eta, kappa, my_u, burnin, nmc, disp)

% Modified from code provided by Hao Wang to allow inference on multiple graphs
% Input parameters:
%   Theta: K x K x S array of initial value for subgroup graph similarity
%   matrices
%   Phi: initial value for S x S platform graph similarity matrix
%   Sig: initial guess for Sigma p_s x p_s x K x S cell array where C inverse = Sig
%   V0:  p_s x p_s x K x S cell array of the  small variance  components;  V0(i, j, k, s) is
%   the small variance of the precision element C(i,j) for group k and
%   platform s
%   V1:  p x p x K x S cell array of the  large variance  components;  V1(i, j, k, s) is
%   the small variance of the precision element C(i,j) for group k and
%   platform s
%   *lambda: hyperparameter for the diagonal element Usually set to be 1;
%   pii: 1 x K x S cell array of prior marginal edge "inclusion probability" for each
%   group and platform
%   n: 1 x K x S cell array of sample sizes for each group and platform
%   *SampCov: p x p x K x S cell array of sample covariance matrix for each group
%   C: p x p x K x S cell array of initial precision matrices for each group
%   nu: p_s x p_s x S cell array of initial values for nu
%   w: k x k matrix of initial values for w
%   alpha: 1 x S vector of shape parameters for Theta slab for each platform
%   beta: 1 x S vector of rate parameters for Theta slab for each platform
%   a: 1 x S vector of first parameters for prior on nu for each platform
%   b: 1 x S vector of second parameter for prior on nu for each platform
%   d: first parameter for prior on w 
%   f: second parameter for prior on w
%   my_u: Bernoulli prior parameter for platform graph similarity indicators
%   eta: shape parameter for Phi slab
%   kappa: rate parameter for Phi slab
%   burnin: number of burn-in iterations
%   nmc: number of MCMC sample iterations saved
%   disp: T/F for whether to print progress to screen
% Output parameters:
%   C_save: p x p x K x S x nmc cell array sample of precision matrix
%   Sig_save: p x p x K x S x nmc cell array sample of covariance matrix
%   adj_save: p x p x K x S x nmc cell array sample of adjacency matrix
%   Phi_save: S x S x nmc sample of network similarity matrix
%   Theta_save: K x K x S x nmc cell array sample of graph similarity matrix
%   ar_gamma: acceptance rate for between-model moves at group level
%   ar_theta: acceptance rate for within-model moves at group level
%   ar_zeta: acceptance rate for between-model moves at platform level
%   ar_phi: acceptance rate for within-model moves at platform level
%   nu_save: p x p x S x nmc sample of edge-specific parameters
%   ar_nu: acceptance rate for nu
%   w_save: K x K x nmc sample of group-specific parameters
%   ar_w: acceptance rate for w


% K is number of sample groups
K = size(Theta{1,:}, 1);

% S is number of platforms
S = size(Phi, 1);

% p is number of variables
[p,foo2] = cellfun(@size,SampCov);

% Set up matrices for return values
w = zeros(K,K,nmc);
for j = 1:S
    C_save_j = zeros(p(j), p(j), K, nmc);
    C_save{j,1} = C_save_j;
    Sig_save{j,1} = C_save_j;
    adj_save{j,1} = C_save_j;
    Theta_save{j,1} = zeros(K, K, nmc);
    nu_save{j,1} = zeros(p(j), p(j), nmc);

    tau{j,1} = V1{j,1};
 
 
 
    ind_noi_all_j = zeros(p(j)-1,p(j));
 
 
    for i = 1:p(j)
       if i==1  
       ind_noi = [2:p(j)]'; 
       elseif i==p[j]
       ind_noi = [1:p(j)-1]'; 
       else
       ind_noi = [1:i-1,i+1:p(j)]';
       end
       ind_noi_all_j(:,i) = ind_noi;
       
    end
    
    ind_noi_all{j,1} = ind_noi_all_j;
 
    pii_RB_j = zeros(p(j), p(j), K);
    pii_mat_j = zeros(p(j), p(j), K);
    for i = 1:K
        pii_mat_j(:,:,i) = pii{j,:}(:,i);
    end
   
    pii_RB{j,1} = pii_RB_j;
    pii_mat{j,1} = pii_mat_j;

    % Acceptance rate for gamma based on number of between model moves at
    % group level
    ar_gamma_j = zeros(K, K);
    ar_gamma{j,1} = ar_gamma_j;
    % Acceptance rate for theta based on number of within model moves at
    % group level
    n_within_model_j = zeros(K, K);
    n_within_model{j,1} = n_within_model_j;
    ar_theta_j = zeros(K, K);
    ar_theta{j,1} = ar_theta_j;
    ar_nu_j = zeros(p(j), p(j));
    ar_nu{j,1} = ar_nu_j;

    % Initial adjacency matrices as defined by intial precision matrices
    adj = cellfun(@(x) abs(x) > 1e-5, C, 'UniformOutput', false);

    % Elementwise product zeros out elements of C that are close to 0
    C{j,:} = C{j,:} .* adj{j,:};
    
end

% Proposal parameters for MH steps at group level
alpha_prop = 1;
beta_prop = 1;
a_prop = 2;
b_prop = 4;

% Proposal parameters for MH steps at platform level
eta_prop = 1;
kappa_prop = 1;
d_prop = 2;
f_prop = 4;

% Acceptance rate for zeta based on number of between model moves at the
% platform level
ar_zeta = zeros(S, S);

%Acceptance rate for phi based on number of within model moves at the
%platform level
n_within_model_platforms = zeros(S, S);
ar_phi = zeros(S, S);
ar_w = zeros(K, K);

Phi_save = zeros(S, S, nmc);
gamma_save = zeros(K, K, S,nmc);
w_save = zeros(K, K, nmc);

% Perform MCMC sampling
for iter = 1:burnin + nmc
    if (disp && mod(iter, 500) == 0)
        fprintf('iter = %d\n', iter);
    end
    
    for cur_platform = 1:S
        
        % Update graph and precision matrix for each group using code from Hao Wang
        for cur_graph = 1:K

            % Sample Sig and C
            for i = 1:p(cur_platform)
                ind_noi = ind_noi_all{cur_platform,:}(:,i);
 
       
        tau_temp = tau{cur_platform,:}(ind_noi,i,cur_graph);
       
        Sig11 = Sig{cur_platform,:}(ind_noi,ind_noi,cur_graph); 
        Sig12 = Sig{cur_platform,:}(ind_noi,i, cur_graph);
      
        invC11 = Sig11 - Sig12*Sig12'/Sig{cur_platform,:}(i,i, cur_graph);
        
        
        Ci = ((SampCov{cur_platform,:}(i,i,cur_graph))+lambda)*invC11+diag(1./tau_temp);  

  
        Ci = (Ci+Ci')./2;       
        Ci_chol = chol(Ci);    
        mu_i = -Ci_chol\(Ci_chol'\SampCov{cur_platform,:}(ind_noi,i,cur_graph));
            epsilon = mu_i+ Ci_chol\randn(p(cur_platform)-1,1);
        
        
            C{cur_platform,:}(ind_noi, i, cur_graph) = epsilon;
            C{cur_platform,:}(i, ind_noi, cur_graph) = epsilon;
        
            a_gam = 0.5*n{cur_platform,:}(:,cur_graph)+1;   
            b_gam = (SampCov{cur_platform,:}(i,i,cur_graph)+lambda)*0.5;
            gam = gamrnd(a_gam,1/b_gam);
         
            c = epsilon'*invC11*epsilon;
        
            C{cur_platform,:}(i,i, cur_graph) = gam+c;
    
            % note epsilon is beta in original Hao Wang code
        
        
            %% Below updating Covariance matrix according to one-column change of precision matrix
            invC11epsilon = invC11*epsilon;
        
            Sig{cur_platform,:}(ind_noi,ind_noi,cur_graph) = invC11+invC11epsilon*invC11epsilon'/gam;
            Sig12 = -invC11epsilon/gam;
            Sig{cur_platform,:}(ind_noi,i,cur_graph) = Sig12;
            Sig{cur_platform,:}(i,ind_noi,cur_graph) = Sig12';
            Sig{cur_platform,:}(i,i,cur_graph) = 1/gam;
 
        
        
              
            v0 = V0{cur_platform,:}(ind_noi,i,cur_graph);
            v1 = V1{cur_platform,:}(ind_noi,i,cur_graph);
         
            mrf_sum=0;
            for foo = 1:K
                mrf_sum=mrf_sum+Theta{cur_platform,:}(cur_graph,foo)*pii_mat{cur_platform,:}(ind_noi,i,foo);
            end
         
    % Applying the MRF prior to probability of inclusion
            pii_mat{cur_platform,:}(ind_noi,i,cur_graph)=exp(pii_mat{cur_platform,:}(ind_noi,i,cur_graph).*(nu{cur_platform,:}(ind_noi,i)+2*mrf_sum))./(1+exp(nu{cur_platform,:}(ind_noi,i)+2*mrf_sum));     
                
            w1 = -0.5*log(v0) -0.5*epsilon.^2./v0+log(1-pii_mat{cur_platform,:}(ind_noi,i,cur_graph));
            w2 = -0.5*log(v1) -0.5*epsilon.^2./v1+log(pii_mat{cur_platform,:}(ind_noi,i,cur_graph));
               
            w_max = max([w1,w2],[],2);
 
            w = exp(w2-w_max)./sum(exp([w1,w2]-repmat(w_max,1,2)),2);
 
            z = (rand(p(cur_platform)-1,1)<w);
        
        
            v = v0;
            v(z) = v1(z);
        
            tau{cur_platform,:}(ind_noi, i, cur_graph) = v;        
            tau{cur_platform,:}(i, ind_noi, cur_graph) = v;
 
            adj{cur_platform,:}(ind_noi, i, cur_graph) = z;
            adj{cur_platform,:}(i, ind_noi, cur_graph) = z;
            
            end

            if iter > burnin
                pii_RB{cur_platform,:} = pii_RB{cur_platform,:} + pii_mat{cur_platform,:}/nmc; 
                C_save{cur_platform,:}(:, :, cur_graph, iter-burnin) = C{cur_platform,:}(:, :, cur_graph);
                adj_save{cur_platform,:}(:, :, cur_graph, iter-burnin) = adj{cur_platform,:}(:, :, cur_graph);
                Sig_save{cur_platform,:}(:,:,cur_graph, iter-burnin) = Sig{cur_platform,:}(:,:,cur_graph); 
            
        
            end
        end
        
   
    % Update the parameters for network relatedness
        for k = 1:K-1
            for m = k+1:K
                % Between model move
                if Theta{cur_platform,:}(k, m) == 0
                    theta_prop = gamrnd(alpha_prop, beta_prop);
                else
                    theta_prop = 0;
                end
            
                Theta_prop{cur_platform,:} = Theta{cur_platform,:};
                Theta_prop{cur_platform,:}(k, m) = theta_prop;
                Theta_prop{cur_platform,:}(m, k) = theta_prop;
            
                % Get terms that are a sum over all edges on log scale
                sum_over_edges = 0;
                for i = 1:p(cur_platform)-1
                    for j = i+1:p(cur_platform)
                        sum_over_edges = sum_over_edges + ...
                            log(calc_mrf_C(Theta{cur_platform,:}, nu{cur_platform,:}(i, j))) + 2 * ...
                            (theta_prop - Theta{cur_platform,:}(m, k)) * adj{cur_platform,:}(i, j, k) * adj{cur_platform,:}(i, j, m) - ...
                            log(calc_mrf_C(Theta_prop{cur_platform,:}, nu{cur_platform,:}(i, j)));
                    end
                end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            % To find metropolis hastings ratio term involving Phi matrix
                for platform=1:S
                    Gamma_mat{platform}=(Theta{platform}~=0);
                end
            
                platform_term=0;
                for s=1:S-1
                    for t=s+1:S
                        platform_term=platform_term + Phi(s,t)*Gamma_mat{t}(k,m);
                    end
                end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Calculate MH ratio on log scale
                if theta_prop == 0
                    log_ar = alpha_prop * log(beta_prop) - log(gamma(alpha_prop)) + ...
                        log(gamma(alpha)) - alpha * log(beta) - ...
                        (alpha - alpha_prop) * log(Theta{cur_platform,:}(m, k)) + ...
                        (beta - beta_prop) * (Theta{cur_platform,:}(m, k)) + sum_over_edges - ...
                        my_w(k,m)-2*platform_term;
                else
                    log_ar = alpha * log(beta) - log(gamma(alpha)) + ...
                        log(gamma(alpha_prop)) - alpha_prop * log(beta_prop) - ...
                        (alpha - alpha_prop) * log(theta_prop) - ...
                        (beta - beta_prop) * theta_prop + sum_over_edges + ...
                        my_w(k,m)+2*platform_term;
                end
            
                % Accept proposal with given probability
                if log_ar > log(unifrnd(0,1))
                    Theta{cur_platform,:}(m, k) = theta_prop;
                    Theta{cur_platform,:}(k, m) = theta_prop;
                
                    % Increment acceptance rate
                    ar_gamma{cur_platform,:}(k, m) = ar_gamma{cur_platform,:}(k, m) + 1 / (burnin + nmc);
                end
            
                % Within model model
                if Theta{cur_platform,:}(k, m) ~= 0
                    n_within_model{cur_platform,:}(k, m) = n_within_model{cur_platform,:}(k, m) + 1;
                    theta_prop = gamrnd(alpha_prop, beta_prop);
                    Theta_prop{cur_platform,:} = Theta{cur_platform,:};
                    Theta_prop{cur_platform,:}(k, m) = theta_prop;
                    Theta_prop{cur_platform,:}(m, k) = theta_prop;
                
                    % Get terms that are a sum over all edges on log scale
                    sum_over_edges = 0;
                    for i = 1:p(cur_platform)-1
                        for j = i+1:p(cur_platform)
                            sum_over_edges = sum_over_edges + ...
                                log(calc_mrf_C(Theta{cur_platform,:}, nu{cur_platform,:}(i, j))) + 2 * ...
                                (theta_prop - Theta{cur_platform,:}(m, k)) * adj{cur_platform,:}(i, j, k) * adj{cur_platform,:}(i, j, m) - ...
                                log(calc_mrf_C(Theta_prop{cur_platform,:}, nu{cur_platform,:}(i, j)));
                        end
                    end
                
                    % Calculate MH ratio on log scale
                    log_theta_ar = (alpha - alpha_prop) * (log(theta_prop) - log(Theta{cur_platform,:}(m, k))) + ...
                        (beta - beta_prop) * (Theta{cur_platform,:}(m, k) - theta_prop) + sum_over_edges;
                
                    % Accept proposal with given probability
                    if log_theta_ar > log(unifrnd(0,1))
                        Theta{cur_platform,:}(m, k) = theta_prop;
                        Theta{cur_platform,:}(k, m) = theta_prop;
                    
                        % Track number of proposals accepted
                        ar_theta{cur_platform,:}(k, m) = ar_theta{cur_platform,:}(k, m) + 1;
                    end
                end
            end
        end
    
        % Generate independent proposals for q from beta(a_prop, b_prop) density
        for i = 1:p(cur_platform)-1
            for j = i+1:p(cur_platform)
                q = betarnd(a_prop, b_prop);
                nu_prop = log(q) - log(1-q);
            
                % Calculate MH ratio on log scale
                % log(p(nu_prop)) - log(p(nu)) + log(q(nu)) - log(q(nu_prop))
                log_nu_ar = (nu_prop - nu{cur_platform,:}(i, j)) * (sum(adj{cur_platform,:}(i, j, :)) + a - a_prop) - ...
                    (a + b - a_prop - b_prop) * log(1 + exp(nu_prop)) - ...
                    log(calc_mrf_C(Theta{cur_platform,:}, nu_prop)) + ...
                    (a + b - a_prop - b_prop) * log(1 + exp(nu{cur_platform,:}(i, j))) + ...
                    log(calc_mrf_C(Theta{cur_platform,:}, nu{cur_platform,:}(i, j)));
            
                if log_nu_ar > log(unifrnd(0,1))
                    nu{cur_platform,:}(i, j) = nu_prop;
                    nu{cur_platform,:}(j, i) = nu_prop;
                    ar_nu{cur_platform,:}(i, j) = ar_nu{cur_platform,:}(i, j) + 1 / (burnin + nmc);
                end
            end
        end
    
        % Retain values for posterior sample
        if iter > burnin
            nu_save{cur_platform,:}(:, :, iter-burnin) = nu{cur_platform,:};
            Theta_save{cur_platform,:}(:, :, iter-burnin) = Theta{cur_platform,:};
        end
    

    % Compute acceptance rate for theta as number of proposals accepted divided
    % by number of within model moves proposed
    for k = 1:K-1
        for m = k+1:K
            ar_theta{cur_platform,:}(k, m) = ar_theta{cur_platform,:}(k, m) / n_within_model{cur_platform,:}(k, m);
        end
    end
    end
    
    
    
    % Update the parameters for platform relatedness
    
    for platform = 1:S
        Gamma_mat{platform}=(Theta{platform}~=0);
    end
    
    for s = 1:S-1
        for t = s+1:S
            % Between model move
            if Phi(s, t) == 0
                phi_prop = gamrnd(eta_prop, kappa_prop);
            else
                phi_prop = 0;
            end
            
            Phi_prop = Phi;
            Phi_prop(s, t) = phi_prop;
            Phi_prop(t, s) = phi_prop;
            
            % Get terms that are a sum over all edges on log scale
            sum_over_edges = 0;
            for i = 1:K-1
                for j = i+1:K
                    sum_over_edges = sum_over_edges + ...
                        log(calc_mrf_C(Phi, my_w(i, j))) + 2 * ...
                        (phi_prop - Phi(t, s)) * Gamma_mat{s}(i, j) * Gamma_mat{t}(i, j) - ...
                        log(calc_mrf_C(Phi_prop, my_w(i, j)));
                end
            end
            
            % Calculate MH ratio on log scale
            if phi_prop == 0
                log_ar = eta_prop * log(kappa_prop) - log(gamma(eta_prop)) + ...
                    log(gamma(eta)) - eta * log(kappa) - ...
                    (eta - eta_prop) * log(Phi(t, s)) + ...
                    (kappa - kappa_prop) * (Phi(t, s)) + sum_over_edges + ...
                    log(1 - my_u) - log(my_u);
            else
                log_ar = eta * log(kappa) - log(gamma(eta)) + ...
                    log(gamma(eta_prop)) - eta_prop * log(kappa_prop) - ...
                    (eta - eta_prop) * log(phi_prop) - ...
                    (kappa - kappa_prop) * phi_prop + sum_over_edges + ...
                    log(my_u) - log(1 - my_u);
            end
            
            % Accept proposal with given probability
            if log_ar > log(unifrnd(0,1))
                Phi(t, s) = phi_prop;
                Phi(s, t) = phi_prop;
                
                % Increment acceptance rate
                ar_zeta(s, t) = ar_zeta(s, t) + 1 / (burnin + nmc);
            end
            
            % Within model move
            if Phi(s, t) ~= 0
                n_within_model_platforms(s, t) = n_within_model_platforms(s, t) + 1;
                phi_prop = gamrnd(eta_prop, kappa_prop);
                Phi_prop = Phi;
                Phi_prop(s, t) = phi_prop;
                Phi_prop(t, s) = phi_prop;
                
                % Get terms that are a sum over all edges on log scale
                sum_over_edges = 0;
                for i = 1:K-1
                    for j = i+1:K
                        sum_over_edges = sum_over_edges + ...
                            log(calc_mrf_C(Phi, my_w(i, j))) + 2 * ...
                            (phi_prop - Phi(t, s)) * Gamma_mat{s}(i, j) * Gamma_mat{t}(i, j) - ...
                            log(calc_mrf_C(Phi_prop, my_w(i, j)));
                    end
                end
                
                % Calculate MH ratio on log scale
                log_phi_ar = (eta - eta_prop) * (log(phi_prop) - log(Phi(t, s))) + ...
                    (kappa - kappa_prop) * (Phi(t, s) - phi_prop) + sum_over_edges;
                
                % Accept proposal with given probability
                if log_phi_ar > log(unifrnd(0,1))
                    Phi(t, s) = phi_prop;
                    Phi(s, t) = phi_prop;
                    
                    % Track number of proposals accepted
                    ar_phi(s, t) = ar_phi(s, t) + 1;
                end
            end
        end
    end
    
    % Generate independent proposals for q from beta(d_prop, f_prop) density
    for i = 1:K-1
        for j = i+1:K
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To define gamma_vector of gamma_km from all platforms
            
            gamma_vector_ij=[];
            for s = 1:S
                gamma_vector_ij(s)=Gamma_mat{s}(i,j);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q_platforms = betarnd(d_prop, f_prop);
            w_prop = log(q_platforms) - log(1-q_platforms);
            
            % Calculate MH ratio on log scale
            %log(p(w_prop)) - log(p(my_w)) + log(q(my_w)) - log(q(w_prop))
            log_w_ar = (w_prop - my_w(i, j)) * (sum(gamma_vector_ij) + d - d_prop) - ...
                (d + f - d_prop - f_prop) * log(1 + exp(w_prop)) - ...
                log(calc_mrf_C(Phi, w_prop)) + ...
                (d + f - d_prop - f_prop) * log(1 + exp(my_w(i, j))) + ...
                log(calc_mrf_C(Phi, my_w(i, j)));
            
            if log_w_ar > log(unifrnd(0,1))
                my_w(i, j) = w_prop;
                my_w(j, i) = w_prop;
                ar_w(i, j) = ar_w(i, j) + 1 / (burnin + nmc);
            end
        end
    end
    
    % Retain values for posterior sample
    if iter > burnin
        w_save(:,:, iter-burnin) = my_w;
        Phi_save(:, :, iter-burnin) = Phi;
    end
end
 
%% burn in loop should end above below

% Compute acceptance rate for phi as number of proposals accepted divided
% by number of within model moves proposed
for s = 1:S-1
    for t = s+1:S
        ar_phi(s, t) = ar_phi(s, t) / n_within_model_platforms(s, t);   
    end
end
    
