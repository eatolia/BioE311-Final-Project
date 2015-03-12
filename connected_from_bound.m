function [connected cluster_list numClusters] = connected_from_bound(bound)
%this function takes in a pairwise matrix that indicates which pairs of
%cells are directly bound to each other(analogous to an adjacency matrix in graph theory), and
%returns a pairwise matrix that indicates which cells are connected, even
%indirectly (connected matrix)
%it uses a tree traversal strategy with an iterative subfunction to do this
%cluster_list returns a list of the clusters that exist (max number of
%clusters = number of cells, only when no cells are bound to each other)
numCells=length(bound);
cluster_list=zeros(numCells);
unsearched=ones(numCells, 1);
numClusters=0;
%iteratively searches bound matrix for cells that are connected to
%cells that are connected to.... until all cells have been searched and
%assigned to a cluster
for index=1:numCells
if unsearched(index)
unsearched(index)=0;
numClusters=numClusters+1;
[cluster unsearched]=find_bound(bound, index, unsearched);
cluster_list(cluster, numClusters)=1;
end
end
%from these clusters, generate the connected matrix, by convention,
%cells are not connected to themselves.
for index=1:numClusters
a=find(cluster_list(:,index));
connected(a, a)=1;
end
connected=connected-eye(numCells);
end