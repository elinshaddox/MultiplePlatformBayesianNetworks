function [K_save,Sigma_save,acc] = wishart_InvA_rnd_ud_rw(df,S,adj,sig0,nmc,burnin);
%% Section 2.2 of Dobra,Lenkoski,Abel


T = chol(inv(S));

T1 = T*diag(1./diag(T));

p=size(S,1); 

indmx= reshape([1:p^2],p,p); % Uppermatrix index
upperind=indmx(triu(indmx,1)>0);       
acc = zeros(p);

A = zeros(p);
A(upperind) = adj(upperind);
nu = sum(A,2);

K_save= zeros(p,p,nmc);
Sigma_save = zeros(p,p,nmc);
Psi_save = K_save;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initial Values %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi = diag(sqrt(chi2rnd(df+nu)));
Psi(find(A==1))  = randn(sum(nu),1);
   
 %%%%%% if i = 1;   
   for j=2:p
    if(A(1,j)==0)
    Psi(1,j) = -sum(Psi(1,1:j-1)'.*T1(1:j-1,j));
%    logJeach = logJeach - Psi(1,j)^2/2;
    end
  end


for i=2:p
    for j=i+1:p
   
    if(A(i,j)==0)
    Psi(i,j) = -sum(Psi(i,i:j-1)'.*T1(i:j-1,j));
    
      
    for r = 1: i-1
        
    Psi(i,j) = Psi(i,j) - 1/Psi(i,i)*(Psi(r,i)+ sum(Psi(r,r:i-1).*T1(r:i-1,i)'))* ...
        (Psi(r,j)+ sum(Psi(r,r:j-1).*T1(r:j-1,j)'));
        
    end
 %   logJeach = logJeach - Psi(i,j)^2/2;   
    end
    end
end
    
    logR = -sum(Psi(:).^2)/2; 

    
%%%%%%%%%%%%%%%%%%%%%%%%    
%%%  Below MCMC %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:(burnin+nmc)
    if(mod(iter,1000)==0)
        fprintf('iter=%d\n',iter);
    end 

  for k = 1:p
      Psi_star_k = rand_truncated_normal_rej_below(Psi(k,k),sig0,0);
      Psi_star = Psi;
      Psi_star(k,k) = Psi_star_k; 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN Completion Lemma 2 Atay-Kayis and Massam 2005 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for j=2:p
    if(A(1,j)==0)
    Psi_star(1,j) = -sum(Psi_star(1,1:j-1)'.*T1(1:j-1,j));
    end
   end
   
   for i=2:p-1
    for j=i+1:p
   
    if(A(i,j)==0)
    Psi_star(i,j) = -sum(Psi_star(i,i:j-1)'.*T1(i:j-1,j));
    
      
    for r = 1: i-1
        
    Psi_star(i,j) = Psi_star(i,j) - 1/Psi_star(i,i)*(Psi_star(r,i)+ sum(Psi_star(r,r:i-1).*T1(r:i-1,i)'))* ...
        (Psi_star(r,j)+ sum(Psi_star(r,r:j-1).*T1(r:j-1,j)'));
        
    end
    end
    end
   end
   
   logR_star = -sum(Psi_star(:).^2)/2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END Completion Lemma 2 Atay-Kayis and Massam 2005 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   prior_ratio = log(normcdf(Psi(k,k),0,sig0))-log(normcdf(Psi_star(k,k),0,sig0));   
   post_ratio = (nu(k)+df-1)*(log(Psi_star(k,k))-log(Psi(k,k))) + logR_star - logR;
   
   if(log(rand(1)) < (post_ratio + prior_ratio))
     
     if(iter>burnin)  
      acc(k,k) = acc(k,k)+1;    
     end
      Psi = Psi_star;
   logR = logR_star;   
   end


  end

%%%  Sample Free Off-Diag elements in Psi
 for k = 1:p-1
   for s = k+1:p   
       
    if(A(k,s)==1)    
      Psi_star_ks = Psi(k,s)+randn(1)*sig0;
      Psi_star = Psi;
      Psi_star(k,s) = Psi_star_ks; 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN Completion Lemma 2 Atay-Kayis and Massam 2005 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for j=s+1:p
    if(A(1,j)==0)
    Psi_star(1,j) = -sum(Psi_star(1,1:j-1)'.*T1(1:j-1,j));
    end
   end
   
   for i=k:p-1
    for j=i+1:p
   
    if(A(i,j)==0)
    Psi_star(i,j) = -sum(Psi_star(i,i:j-1)'.*T1(i:j-1,j));
    
      
    for r = 1: i-1        
    Psi_star(i,j) = Psi_star(i,j) - 1/Psi_star(i,i)*(Psi_star(r,i)+ sum(Psi_star(r,r:i-1).*T1(r:i-1,i)'))* ...
    (Psi_star(r,j)+ sum(Psi_star(r,r:j-1).*T1(r:j-1,j)'));       
    end
    end
    
    end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END Completion Lemma 2 Atay-Kayis and Massam 2005 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   logR_star = -sum(Psi_star(:).^2)/2; 
   post_ratio = logR_star - logR;
   
   if(log(rand(1)) < post_ratio)
       if(iter>burnin)  
           acc(k,s) = acc(k,s)+1;    
       end
   Psi = Psi_star;
   logR = logR_star;   
   end

    end
  
   end
 end

Phi = Psi*T;
K = Phi'*Phi;

K = K.*adj;

Sigma = inv(K);


if(iter>burnin)
Psi_save(:,:,iter-burnin) = Psi;
K_save(:,:,iter-burnin) = K;
Sigma_save(:,:,iter-burnin) = Sigma;
end

end