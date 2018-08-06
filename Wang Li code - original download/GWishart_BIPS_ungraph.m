function [C_save,Sig_save,adj_save] = GWishart_BIPS_ungraph(b_prior,D_prior,n,S,beta,C,burnin,nmc)
%     Gibbs sampler for Bayesian GGM determination without approximating
%     normalizing constant of the G Wishart distribution
%     Input:     
%        b_prior : prior d.f.
%        D_prior :  prior location
%        n : sample size
%        S : Y'*Y matrix or corr(Y)*n or cov(Y)*n ;
%        beta: prior edge inclusion prob.
%        C: initial partial covariance matrix;
%        burnin, nmc : number of MCMC burnins and saved samples
%     Output:
%        C_save: saved precision matrices;
%        Sig_save: saved covariance matrices;
%        adj_save: saved adjacency matrices.
%     Reference: Wang and Li (2011) "Efficient Gaussian Graphical Model Determination without 
%      Approximating Normalizing Constant of the G-Wishart distribution "

[p] = size(D_prior,1); 
b_post = b_prior+n; bG = b_post;
D_post = D_prior+S; DG = D_post;

C_save = zeros(p,p,nmc);
Sig_save = C_save;
adj_save = C_save;
 

adj = abs(C)>1e-5;
C = C.*adj; 
Sig = inv(C);

for iter = 1: burnin+nmc    
            
       if(mod(iter,1000)==0)
        fprintf('iter=%d nedge = %d\n',iter,(sum(adj(:))-p)/2);
       end
     
       
       %% Sample off-diagonal elements     
for i = 1:p-1    
    for j = i+1:p
       
        
        reorder = [1:p];
        reorder([i,j])=[]; 
        reorder = [reorder,i,j];
        
        C_reorder = C(reorder,reorder);
        
        R = chol(C_reorder);
        b = R(p-1,p-1);
        
       %% If no edge
       R(p-1,p) = -sum(R(1:p-2,p-1).*R(1:p-2,p))/R(p-1,p-1); %  Ref: Eq. (7) of Atay-Kayis and Massam
       log_p_gamma0 = -S(j,j)*R(p-1,p)^2/2 - S(i,j)*R(p-1,p-1)*R(p-1,p)+log(1-beta);
       
       %% If edge  
        m_prior = -b*D_prior(i,j)/D_prior(j,j);        
        m_post = -b*D_post(i,j)/D_post(j,j);        
        
        sig_post = 1/sqrt(D_post(j,j));
       
       
        log_p_gamma1 = 1/2*(log(D_prior(j,j))-log(D_post(j,j)))...
            +D_post(j,j)/2*m_post^2 - D_prior(j,j)/2*m_prior^2+log(beta);
       
        w = [log_p_gamma0,log_p_gamma1];
        w = exp(w-max(w))./sum(exp(w-max(w)));
        
        adj(i,j) = rand(1)<w(2);
        adj(j,i) = adj(i,j);
  
        
       R(p,p) = sqrt(gamrnd(bG/2,2/DG(j,j)));  
       
        if(adj(i,j)==1)       
       R(p-1,p) = randn(1)*sig_post + m_post;        
       C_updated = R(:,end-1:end)'*R(:,end);
       C(i,j) = C_updated(1);
       C(j,i) = C_updated(1);
       C(j,j) = C_updated(2);
        else
       C_updated = R(:,end)'*R(:,end);
       C(j,j) = C_updated(1);
       C(i,j) = 0;
       C(j,i) = 0;     
        end
        

    end
end

   Sig = inv(C); %% This inversion can be replaced by low rank updates in the above loop;

IsolatedNodeId = find(sum(adj)==1); % isolated node, i.e. nodes that are not connected to any other nodes
INsize = length(IsolatedNodeId);

%%%  Sample isolated nodes
for i = 1:INsize
    
    cliqueid = IsolatedNodeId(i);    
    K_c = wishrnd(inv(DG(cliqueid,cliqueid)),bG);
    C(cliqueid,cliqueid) = K_c;
    Sig(cliqueid,cliqueid) = 1/K_c;
      
end



for i = 1:p-1
    for j = i+1:p
      
    if(adj(i,j)==1)    
     cliqueid = [i,j];
     cliquesize = 2;
     A = wishrnd(inv(DG(cliqueid,cliqueid)),bG+cliquesize-1);
            
        C_12 = C(cliqueid,:);  C_12(:,cliqueid)=[];

        Sig_12 = Sig(cliqueid,:);  Sig_12(:,cliqueid)=[];        
        Sig_22 = Sig; Sig_22(cliqueid,:)=[]; Sig_22(:,cliqueid)=[];
        invSig_11 = inv(Sig(cliqueid,cliqueid));
        invSig_11 = (invSig_11+ invSig_11')/2;
        invC_22 = Sig_22 - Sig_12'*invSig_11*Sig_12;
        
          K_c = A + C_12*invC_22*C_12';
        
        
        K_c = (K_c + K_c')/2; % Numerical stable
        
        Delta = inv( C(cliqueid,cliqueid)-K_c);
           
         
        C(cliqueid,cliqueid) = K_c;
        
       %% Rank-2 update Sigma 
        Sig_bb = Sig(cliqueid,cliqueid);        
        aa = inv(Delta- Sig_bb);
        aa = (aa+aa')/2;
        Sig = Sig + Sig(:,cliqueid)*aa*Sig(:,cliqueid)';
        
    end

    end
end


       if iter >burnin
            Sig_save(:,:,iter-burnin) = Sig; 
            C_save(:,:,iter-burnin) = C;
            adj_save(:,:,iter-burnin) = adj;
       end



end



