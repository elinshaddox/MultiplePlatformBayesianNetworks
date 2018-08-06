function f=log_multi_gamma(p,n)
%  MULTI_GAMMA returns multivariate gamma function at value n,p 
%   \Gamma_p(n/2)= \pi^{p(p-1)/4}\Pi_{j=1}^p \Gamma\left[ (n+1-j)/2\right]. 

f=(p*(p-1)/4)*log(pi);
for j=1:p
f=f+gammaln(n+(1-j)/2);
end
