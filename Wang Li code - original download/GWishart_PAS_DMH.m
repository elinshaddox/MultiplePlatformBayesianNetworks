function [C_save,Sig_save,adj_save] = GWishart_PAS_DMH(b_prior,D_prior,n,S,C,beta,burnin,nmc)
%     Sample C from Gwishart distribution with  density:
%          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
%     where     
%      (1)  bG : d.f.
%      (2)  DG: location
%     C: initial partial covariance matrix;
%     burnin, nmc : number of MCMC burnins and saved samples

[p] = size(D_prior,1); 
b_post = b_prior+n; 
D_post = D_prior+S; 
C_save = zeros(p,p,nmc);
Sig_save = C_save;
adj_save = C_save;
 

adj = abs(C)>1e-5;
nedge = (sum(adj(:))-p)/2;
C = C.*adj; 

    

for iter = 1: burnin+nmc    
            
       if(mod(iter,1000)==0 & iter>burnin)
        fprintf('iter=%d nedge = %d\n',iter,nedge); 
        mean(adj_save(:,:,1:iter-burnin-1),3)        
       end
  
       
       %% Sample off-diagonal elements     
for i = 1:p-1    
    for j = i+1:p
       
        nedge = (sum(adj(:))-p)/2;
%         %%%%%%%%%%% If no edge
%        adj0 = adj; adj0(i,j)=0; adj0(j,i)=0;
%        
%        if mcard(adj0)==1
%            log_c0 = log_hiwishart_InvA_const(makedecompgraph(adj0),b_prior,D_prior); 
%        else
%            log_c0 = -log_GWishart_ud_const_mc(b_prior,D_prior,adj0,50);
%        end
%        
%        
%        %%%%%%%%%%%%% If edge  
% 
%        adj1 = adj; adj1(i,j)=1; adj1(j,i)=1;
%        if mcard(adj1)==1
%            log_c1 = log_hiwishart_InvA_const(makedecompgraph(adj1),b_prior,D_prior); 
%        else
%            log_c1 = -log_GWishart_ud_const_mc(b_prior,D_prior,adj1,50);
%        end
         
        
        
        
        w = log_H(b_prior,D_prior,n,S,C,i,j)+log(1-beta)-log(beta);%+log_c0-log_c1;
        
        
        w = 1/(exp(w)+1);

%         current_ij = rand(1)<w;
%         adj(i,j) = current_ij; adj(j,i) = current_ij;
%         
%         [C,Sig] = GWishart_NOij_Gibbs(b_post,D_post,adj,C,i,j,current_ij,0,0);
        
        current_ij = adj(i,j);        
        propose_ij = rand(1)<w;
 
        if(propose_ij ~= current_ij)
        
        [C_prop,Sig_prop] = GWishart_NOij_Gibbs(b_prior,D_prior,adj,C,i,j,propose_ij,0,1);
        
        r2 = log_GWishart_NOij_pdf(b_prior,D_prior,C_prop,i,j,current_ij) ...
             -log_GWishart_NOij_pdf(b_prior,D_prior,C_prop,i,j,propose_ij);

         
         if( log(rand(1))<r2 )
         
         adj(i,j) = propose_ij; adj(j,i) = propose_ij;
         current_ij = propose_ij;
        end
         
       end % Second MH
        
        [C,Sig] = GWishart_NOij_Gibbs(b_post,D_post,adj,C,i,j,current_ij,0,0);
        
        
    end
end

 

 %%%%%%%%%  Update C and Sigma given graph
      [C,Sig] = GWishart_BIPS_maximumClique(b_post,D_post,adj,C,0,1);



       if iter >burnin
            Sig_save(:,:,iter-burnin) = Sig; 
            C_save(:,:,iter-burnin) = C;
            adj_save(:,:,iter-burnin) = adj;
       end



end
