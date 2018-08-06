function [C,Sig] = GWishart_NOij_Gibbs(bG,DG,adj,C,i,j,edgeij,burnin,nmc)

%   Start from C\(i,j)  
%
%     Sample C from Gwishart distribution with  density:
%          p(C) \propto |C|^{(bG-2)/2} exp(-1/2 tr(C DG))
%     where     
%      (1)  bG : d.f.
%      (2)  DG: location
%      (3)  adj: adjacency matrix 
%     C: initial partial covariance matrix;
%     burnin, nmc : number of MCMC burnins and saved samples



[p] = size(DG,1); 
%C_save = zeros(p,p,nmc);
%Sig_save = C_save;

% C = C.*adj; 
% Sig = inv(C);








        
       %% If no edge
       
       if(edgeij==0)
       
% %        adj(i,j)=0;
% %        adj(j,i)=0;
% %        
% %        R(p-1,p) = -sum(R(1:p-2,p-1).*R(1:p-2,p))/R(p-1,p-1); %  Ref: Eq. (7) of Atay-Kayis and Massam
% %        R(p,p) = sqrt(gamrnd(bG/2,2/DG(j,j)));  
% %        C_updated = R(:,end)'*R(:,end);
% %        C(j,j) = C_updated(1);
% %        C(i,j) = 0;
% %        C(j,i) = 0; 
 


       C(i,j) = 0;  C(j,i) = 0;
        l = inv(DG(j,j));
        A = wishrnd(l,bG);
            
       C_12 = C(j,:);  C_12(:,j)=[];
       C_22 = C; C_22(j,:)=[]; C_22(:,j)=[];
       invC_22 = inv(C_22);   

        c = C_12*invC_22*C_12';
        
               
        C(j,j) = A + c;

        
        
        
       
       else
           
        reorder = [1:p];
        reorder([i,j])=[]; 
        reorder = [reorder,i,j];
        
        C_reorder = C(reorder,reorder);
        
        R = chol(C_reorder);
        b = R(p-1,p-1);


        m_post = -b*DG(i,j)/DG(j,j);               
        sig_post = 1/sqrt(DG(j,j));
       
           
           
       adj(i,j)=1;
       adj(j,i)=1;
       
       R(p-1,p) = randn(1)*sig_post + m_post; 
       R(p,p) = sqrt(gamrnd(bG/2,2/DG(j,j)));  
       C_updated = R(:,end-1:end)'*R(:,end);
       C(i,j) = C_updated(1);
       C(j,i) = C_updated(1);
       C(j,j) = C_updated(2);
    
%      cliqueid = [i,j];
%      cliquesize = 2;
%      l = inv(DG(cliqueid,cliqueid));
%      l = (l+l')/2;
%      A = wishrnd(l,bG+cliquesize-1);
%             
%       
%         C_12 = C(cliqueid,:);  C_12(:,cliqueid)=[];
%         C_22 = C; C_22(cliqueid,:)=[]; C_22(:,cliqueid)=[];
%         invC_22 = inv(C_22);
%                 
%         K_c = A + C_12*invC_22*C_12';        
%         K_c = (K_c + K_c')/2; % Numerical stable
%         
%          
%         C(cliqueid,cliqueid) = K_c;
        
        
           
       end
       
     
    Sig = inv(C);
     

    
    
IsolatedNodeId = find(sum(adj)==1); % isolated node, i.e. nodes that are not connected to any other nodes
INsize = length(IsolatedNodeId);
   
    
for iter = 1: burnin+nmc    
            


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
     l = inv(DG(cliqueid,cliqueid));
     l = (l+l')/2;
     A = wishrnd(l,bG+cliquesize-1);
            
        C_12 = C(cliqueid,:);  C_12(:,cliqueid)=[];

        Sig_12 = Sig(cliqueid,:);  Sig_12(:,cliqueid)=[];        
        Sig_22 = Sig; Sig_22(cliqueid,:)=[]; Sig_22(:,cliqueid)=[];
        invSig_11 = inv(Sig(cliqueid,cliqueid));
        invSig_11 = (invSig_11+ invSig_11')/2;
        invC_22 = Sig_22 - Sig_12'*invSig_11*Sig_12;
        
          K_c = A + C_12*invC_22*C_12';
        
   %      K_c = A + C_12*inv(C_22)*C_12';
        
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

    %  Sig = inv(C); 




end
