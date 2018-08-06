function y = log_J(h,B,a11)


% h = 3;
% 
% A = wishrnd(eye(2),5); a11 = A(1,1);
% 
% B = wishrnd(eye(2),5)

y = 1/2*log(2*pi/B(2,2))-log_iwishart_InvA_const(h,B(2,2))+(h-1)/2*log(a11)-1/2*(B(1,1)-B(1,2)^2/B(2,2))*a11;

% p = -y+(h-2)/2*log(det(A))-trace(B*A)/2
% 
% 
% 
% 
% V = inv(B); D_priorii = 1/V(1,1);
% log_GWishart_pdf(A,h,B,ones(2),10,1) - log_GWishart_pdf(A(1,1),h+1,D_priorii(1,1),1,10,1)