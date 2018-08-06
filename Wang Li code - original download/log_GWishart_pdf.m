function log_pdf = log_GWishart_pdf(K,b,D,adj,N,IsDec)
% log density function of G-Wishart(b,D,adj):
%
%  p(K) = I(b,D)^{-1} |K|^{(d-2)/2} exp(-trace(K D)/2) using Monte Carlo sample of size N;
%
% Reference: Atay-Kayis and Massam 2005 Biometrika 
%
% INPUT:   b,D,adj,IsDec
% OUTPUT:  log( p(K) ) = log[ I(b,D)^{-1} |K|^{(d-2)/2} exp(-trace(K D)/2) ]; using Monte Carlo sample of size N;

log_pdf_unnormalized = (b-2)/2*log(det(K))-trace(D*K)/2;
if IsDec  % If graph is decomposable 
    
log_pdf = log_pdf_unnormalized+log_hiwishart_InvA_const(makedecompgraph(adj),b,D);    

else
    
log_pdf = log_pdf_unnormalized - log_GWishart_ud_const_mc(b,D,adj,N);
    
    
end