function l = log_mY_NOij_pdf(S,n,b_prior,D_prior,C,i,j,edgeij)

% Compute p(Y | C\c(i,j) ) the normalizing constant of G-Wishart



if edgeij ==0
   
        C(i,j) = 0; C(j,i)=0; %adj0=adj; adj0(i,j) = 0; adj0(j,i)=0; 
 
       
        C_12 = C(j,:);  C_12(:,j)=[];
        C_22 = C; C_22(j,:)=[]; C_22(:,j)=[];
        invC_22 = inv(C_22);   
        c = C_12*invC_22*C_12';

        %        A = 1; 
%        C(j,j) = A + c; 
        
%        log_f = n/2*log(det(C))-trace(S*C)/2; 
%        logKjj_prior = log_GWishart_pdf(A,b_prior,D_prior(j,j),1,1,1);
% 
%        logKjj_post = log_GWishart_pdf(A,b_prior+n,D_prior(j,j)+S(j,j),1,1,1);
% 
%        l = log_f + logKjj_prior - logKjj_post;
       
       C_new = C; C_new(j,j) = c;
       l = log_iwishart_InvA_const(b_prior,D_prior(j,j))-log_iwishart_InvA_const(b_prior+n,D_prior(j,j)+S(j,j)) ...
       +n/2*log(det(C_22)) - 1/2*trace(S*C_new);

else
    
    
        cliqueid = [i,j]; %adj1=adj; adj1(i,j)=1; adj1(j,i)=1;
       
    
        C_12 = C(cliqueid,:);  C_12(:,cliqueid)=[];
        C_22 = C; C_22(cliqueid,:)=[]; C_22(:,cliqueid)=[];
        invC_22 = inv(C_22);
        A = C(cliqueid,cliqueid)-C_12*invC_22*C_12';        A = (A'+A)/2;                    

        
        log_f = n/2*log(det(C))-trace(S*C)/2;

        logK2by2_prior = log_GWishart_pdf(A,b_prior,D_prior(cliqueid,cliqueid),ones(2),10,1);
        
        
        V = inv(D_prior(cliqueid,cliqueid)); D_priorii = 1/V(1,1);
        logKii_prior = log_GWishart_pdf(A(1,1),b_prior+1,D_priorii,1,10,1);
%        logKii_prior = log_GWishart_pdf(A(1,1),b_prior,D_prior(i,i),1,10,1);

        logK2by2_post = log_GWishart_pdf(A,b_prior+n,D_prior(cliqueid,cliqueid)+S(cliqueid,cliqueid),ones(2),10,1);

        V = inv(D_prior(cliqueid,cliqueid)+S(cliqueid,cliqueid)); D_postii = 1/V(1,1);
        logKii_post = log_GWishart_pdf(A(1,1),b_prior+n+1,D_postii,1,10,1);        
%        logKii_post = log_GWishart_pdf(A(1,1),b_prior+n,D_prior(i,i)+S(i,i),1,10,1);

        l = log_f+logK2by2_prior-logKii_prior-logK2by2_post+logKii_post;

    
    
end
