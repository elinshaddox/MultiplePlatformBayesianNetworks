function [C_save,Sig_save] = GWishart_BIPS_pairwise(bG,DG,adj,C,burnin,nmc)
%     BIPS algorithm for sampling from G-Wishart distribution with  density:
%          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
%     where     
%      (1)  bG : d.f.
%      (2)  DG: location
%      (3)  adj: adjacency matrix 
%     C: initial precision matrix;
%     burnin, nmc : number of MCMC burnins and saved samples
%     Reference: Wang and Li (2011) "Efficient Gaussian Graphical Model Determination without 
%      Approximating Normalizing Constant of the G-Wishart distribution "


[p] = size(DG,1); 
C_save = zeros(p,p,nmc);
Sig_save = C_save;

C = C.*adj; 
Sig = inv(C);

IsolatedNodeId = find(sum(adj)==1); % isolated node, i.e. nodes that are not connected to any other nodes
INsize = length(IsolatedNodeId);

for iter = 1: burnin+nmc    
            
       if(mod(iter,1000)==0)
        fprintf('iter = %d \n',iter);
       end


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
        
       %% Rank 2 update Sig 
        Sig_bb = Sig(cliqueid,cliqueid);        
        aa = inv( Delta- Sig_bb);
        aa = (aa+aa')/2;
        Sig = Sig + Sig(:,cliqueid)*aa*Sig(:,cliqueid)';
        
    end

    end
end


       if iter >burnin
           
            Sig_save(:,:,iter-burnin) = Sig; 
            C_save(:,:,iter-burnin) = C;
       end



end



