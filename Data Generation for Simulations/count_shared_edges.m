function [ count_shared ] = count_shared_edges(adj1, adj2)
p = size(adj1);
count_shared=0;
for i =1:p(1)
    for j =1:p(2)
        if (adj1(i,j)==adj2(i,j)) && (adj1(i,j)==1)
            count_shared = count_shared+1;
        end
    end
end
count_shared =(count_shared-p(1))/2
end



