function l = log_GWishart_NOij_pdf(b_prior,D_prior,C,i,j,edgeij)

% Compute log p(C\c(i,j) ) upto the normalizing constant of G-Wishart



if edgeij ==0

    
    
        C(i,j) = 0; C(j,i)=0; %adj0=adj; adj0(i,j) = 0; adj0(j,i)=0; 
 
       
        C_12 = C(j,:);  C_12(:,j)=[];
        C_22 = C; C_22(j,:)=[]; C_22(:,j)=[];
        invC_22 = inv(C_22);   

        c = C_12*invC_22*C_12';

%        A = 1; 
%        C(j,j) = A + C_12*invC_22*C_12';   
%         
%        log_Joint = (b_prior-2)/2*log(det(C))-trace(D_prior*C)/2; 
%        logKii = log_GWishart_pdf(A,b_prior,D_prior(j,j),1,1,1);
%        l = log_Joint-logKii;
       
       C_new = C; C_new(j,j)=c;
       l = -log_iwishart_InvA_const(b_prior,D_prior(j,j))+(b_prior-2)/2*log(det(C_22))-trace(D_prior*C_new)/2;  



else
    
    
        cliqueid = [i,j]; %adj1=adj; adj1(i,j)=1; adj1(j,i)=1;
       
    
        C_12 = C(cliqueid,:);  C_12(:,cliqueid)=[];
        C_22 = C; C_22(cliqueid,:)=[]; C_22(:,cliqueid)=[];
        invC_22 = inv(C_22);
        A = C(cliqueid,cliqueid)-C_12*invC_22*C_12';        A = (A'+A)/2;                    

        
        log_Joint = (b_prior-2)/2*log(det(C))-trace(D_prior*C)/2;

        logK2by2 = log_GWishart_pdf(A,b_prior,D_prior(cliqueid,cliqueid),ones(2),10,1);
        
        
        V = inv(D_prior(cliqueid,cliqueid)); D_priorii = 1/V(1,1);
        logKii = log_GWishart_pdf(A(1,1),b_prior+1,D_priorii,1,10,1);


        l = log_Joint+logKii-logK2by2;

    
    
end
