function [G] = makedecompgraph(A,nodesIDs,nodenames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  an adjacency  matrix A of a decomposable graph G
% Output: cell array G containing the cliques and separators of G   
%        nodeIDs and nodenames are optional inputs                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Adj = full(A); p=size(Adj,1);
if nargin<3, nodenames=int2str((1:p)'); end
if nargin>1, order = nodesIDs; else, order = 1:p; end
i=1;
while i<p
    [a,b] = max(sum(Adj(1:i,i+1:end),1)); 
    order([i+1 b+i]) = order([b+i i+1]);
    i=i+1;
    Adj = full(A); Adj=Adj(order,order);
end
numberofcliques = 1; C(1).ID = [1]; i=2;
while i<=p
    if(sum(Adj(i,C(numberofcliques).ID))==length(C(numberofcliques).ID))
      C(numberofcliques).ID = [C(numberofcliques).ID i];
    else
        C(numberofcliques).dim = length(C(numberofcliques).ID);
        numberofcliques = numberofcliques + 1;
        C(numberofcliques).ID = union(i,find(Adj(i,1:i)==1));
    end
    i=i+1;
end
C(numberofcliques).dim = length(C(numberofcliques).ID);
for i=1:numberofcliques
    C(i).ID = sort(order(C(i).ID));
    C(i).names = nodenames(C(i).ID,:);
end
UN = C(1).ID; S(1).ID=[]; S(1).names=[]; S(1).dim=[];
for i=2:numberofcliques
    S(i).ID    = intersect(UN,C(i).ID);
    S(i).dim   = length(S(i).ID);
    S(i).names = nodenames(S(i).ID,:);
    S(i).names = [];
    UN = union(UN,C(i).ID);
end
C(1).names = nodenames(C(1).ID,:);  S(1).names = []; 
G{1}=C; G{2}=S;


















