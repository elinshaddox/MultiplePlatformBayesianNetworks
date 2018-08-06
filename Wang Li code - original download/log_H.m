function h = log_H(b_prior,D_prior,n,S,C,i,j)


%%%%%%%%%  (i,j) = 0
        e = j; C0 = C; C0(i,j)=0; C0(j,i)=0;      
        C_12 = C0(e,:);  C_12(:,e)=[];
        C_22 = C0; C_22(e,:)=[]; C_22(:,e)=[];
        invC_22 = inv(C_22);
        c = C_12*invC_22*C_12';
        
       
        C0_ij = [C(i,i), 0; 0, c];        
  %      detC0j = det(C_22);

        
%%%%%%%%%  (i,j) = 1
        e = [i,j];     
        C_12 = C(e,:);  C_12(:,e)=[];
        C_22 = C; C_22(e,:)=[]; C_22(:,e)=[];
        invC_22 = inv(C_22);
        Ce = C_12*invC_22*C_12';
        A = C(e,e)-Ce;       % A = (A'+A)/2;                    

                
        a11 = A(1,1);
        
        C1_ij = Ce; 
        
        
        
        
        b_post = b_prior+n; D_post = D_prior+S;
        h = -log_iwishart_InvA_const(b_post,D_post(j,j))-log_J(b_post,D_post(e,e),a11) ... 
            +(n+b_prior-2)/2*(log(a11))  - trace(D_post(e,e)*(C0_ij-C1_ij))/2;
        
        

   