function f = log_hiwishart_InvA_const(G,df,S)
%  Generates the normalizing constant "cons."  for a HyperInv-Wishart(G,df,S).
% Nonsingular  pdf is p(clique) = cons. |K|^{-(df+2|p|)/2} exp(-trace(inv(K) S)/2) 

[p,p]=size(S);

cliques=G{1};
separators=G{2};
numberofcliques=length(cliques);
f=log_iwishart_InvA_const(df,S(cliques(1).ID,cliques(1).ID));

if numberofcliques > 1
for i=2:numberofcliques
    f=f+log_iwishart_InvA_const(df,S(cliques(i).ID,cliques(i).ID))-log_iwishart_InvA_const(df,S(separators(i).ID,separators(i).ID));
end

elseif numberofcliques==1
    f=f;
end