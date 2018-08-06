function [C_save,Sig_save] = GWishart_BIPS_maximumClique(bG,DG,adj,C,burnin,nmc)
%     Sample C from Gwishart distribution with  density:
%          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
%     where     
%      (1)  bG : d.f.
%      (2)  DG: location
%      (3)  adj: adjacency matrix 
%     C: initial partial covariance matrix;
%     burnin, nmc : number of MCMC burnins and saved samples

adj0 = adj - diag(diag(adj));% Create adjacency matrix with diagonal elements zero
cliqueMatrix = maximalCliquesBK(adj0);

numberofcliques = size(cliqueMatrix,2);


[p] = size(DG,1); 
C_save = zeros(p,p,nmc);
Sig_save = C_save;

 
C = C.*adj; 
for iter = 1: burnin+nmc    
            
%        if(mod(iter,5000)==0)
%         fprintf('iter=%d\n',iter);
%        end
       
       
for i = 1:numberofcliques
      
        
   cliqueid = find(cliqueMatrix(:,i)==1);
   cliquesize = length(cliqueid);
   A = wishrnd(inv(DG(cliqueid,cliqueid)),bG+cliquesize-1);

        
    
            
        C_12 = C(cliqueid,:);  C_12(:,cliqueid)=[];
        C_22 = C; C_22(cliqueid,:)=[]; C_22(:,cliqueid)=[];
        
        C121 = C_12*inv(C_22)*C_12'; C121 = (C121+C121')/2;
  
            
        K_c = A + C121;
        
        
        C(cliqueid,cliqueid) = (K_c+K_c')/2;
        

end


       if iter >burnin
            Sig_save(:,:,iter-burnin) = inv(C); 
            C_save(:,:,iter-burnin) = C;
       end



end



