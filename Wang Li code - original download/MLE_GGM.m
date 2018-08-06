function [W,M]=MLE_GGM(S,adj,iterMax,epsilon,IsDec)
%% A modified regression algorithm for estimation of an Undirected GGM with
%% known structure adjacency matrix adj; (reference: the elements of statistical
%% learning; 2rd Algorithm 19.1)


if IsDec==1 % decomposable graph %
    G = makedecompgraph(adj);
    [W] = Gcompletion(S,adj,G);
    M = inv(W);
    
else
    
adj = adj - diag(diag(adj)); % set diagnoal elements of adjacency matrix to be 0;
[p]= size(S,1);
    
W = S;
M = inv(W);
for iter = 1:iterMax
     
    W_previous = W;
    for j=1:p
        ind = find(adj(:,j)==1);
        W11 = W ;
        W11(j,:) = [];
        W11(:,j)=[];
        
        W11star = W(ind,ind);
  
        s12star = W(ind,j);        
        betahat = zeros(p,1);
        betastar = inv(W11star)*s12star;
        betahat(ind) = betastar;
        betahat(j)=[];
        w12 = W11*betahat;
        
        if(j>1)
        W(1:j-1,j)=w12(1:j-1);
        W(j,1:j-1)=w12(1:j-1)';      
        end
        
        if(j<p)      
        W(j+1:end,j)=w12(j:end);
        W(j,j+1:end)=w12(j:end)';
        end
    end
    W_current = W;
    
    Wdiff = abs(W_current - W_previous);
    if(max(Wdiff(:))<epsilon);
      % Solve for precision matrix;   
        break;
    end
        
end % non-decomposable graph 
end
M = inv(W);
