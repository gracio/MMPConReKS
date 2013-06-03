function [val mi mj]=cluster_matching(clusterList1,clusterList2)
% BIPARTITE_MATCHING Solve a maximum weight bipartite matching problem
%
% [val m1 m2]=bipartite_matching(A) for a rectangular matrix A 
% [val m1 m2 mi]=bipartite_matching(x,ei,ej,n,m) for a matrix stored
% in triplet format.  This call also returns a match


% construct A pairwise 
for i=1:length(clusterList1)
    for j=1:length(clusterList2)
        %% using Tanimono set-similarity
        %Tanimono(i,j) = 1 - (length(clusterList1{i}) + length(clusterList2{j}) - 2*length(intersect(clusterList1{i},clusterList2{j})) )/(length(clusterList1{i}) + length(clusterList2{j}) - length(intersect(clusterList1{i},clusterList2{j})) );
        % alternatively using Jaccard distance (Same as Tanimono
        % set-similarity in Yiannis' paper)
        Jaccard(i,j) = length(intersect(clusterList1{i},clusterList2{j}))/length(union(clusterList1{i},clusterList2{j}));
    end
end



% run bipartite_matching to figure out best match, and cost
[val mi mj]=bipartite_matching(Jaccard);