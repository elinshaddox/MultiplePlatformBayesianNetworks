%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% (I) Sample from G-Wishart Distribution for given graphs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% e.g. (1) Sparse Circle graph:  Dobra Lenkoski and Abel (2011, JASA) Example      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = 10;  bG = 103;
A = toeplitz([1,0.5,zeros(1,p-2)]);
A(1,p)=0.4; A(p,1) = 0.4;

D = eye(p) + 100*inv(A);
adj = abs(A)>0.001;
DG = MLE_GGM(D,adj,200,0.0001,0); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% e.g. (2) dense graph      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p=10; alpha = 0.3; J0 = 0.5*eye(p); B0 = rand(p,p)<alpha;
    for i = 1:p, for j = 1:i, if B0(i,j) && i~=j, J0(i,j) = 0.5; end; end; end
    J0 = J0 + J0'; w = eig(J0); delta = (p*min(w)-max(w))/(1-p); J0 = J0 + delta*eye(p);
    D = eye(p) + 100*inv(J0);
    DG = D;
    adj = abs(J0)>1e-4; bG= 103;
    
    


burnin = 1000; nmc = 1000;
%%% algorithm (i) Edgewise Gibbs algorithm  %%%%%%%%%%%%%%%%
    C = eye(p); % Initial value
 [C_save_edgewise,Sig_save_edgewise] = GWishart_BIPS_pairwise(bG,DG,adj,C,burnin,nmc);


    
%%%  Maximum clique algorithm  %%%%%%%%%%%%%%%%
    C = eye(p); % Initial value
 [C_save_maxC,Sig_save_maxC] = GWishart_BIPS_maximumClique(bG,DG,adj,C,burnin,nmc);

    
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% (II) Graphical Model Deterimation Without Approximating Normalizing 
%%%%% of G-Wishart Distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% e.g. 1:  a 6-node circle graph example 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    p = 6; tedge = p*(p-1)/2; beta = 0.5; 
    indmx = reshape([1:p^2],p,p); 
    upperind = indmx(triu(indmx,1)>0); 

    b_prior = 3; D_prior = eye(p); n = 3*p; b_post = b_prior+n;
    A = toeplitz([1,0.5,zeros(1,p-2)]); A(1,p)=0.4; A(p,1) = 0.4;
    S = inv(A)*n; 
    D_post = D_prior + S;
    adjTrue = abs(A)>0.001;  
    
 


burnin  = 3000; nmc = 3000;     C = eye(p); 
tic
[C_save,Sig_save,adj_save] = GWishart_PAS_DMH(b_prior,D_prior,n,S,C,beta,burnin,nmc);
toc



