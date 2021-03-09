function [centroids, n_cluster_values, dimension] = ...
    getClusterCentroids(data2cluster, clusterSolution)

n_cluster_values = unique(clusterSolution);

% as unique function considers NaN values as unique the following code was
% added to remove the NaN:
n_cluster_values = (n_cluster_values(~isnan(n_cluster_values)));


ncluster = numel(n_cluster_values);
dimension = size(data2cluster,2);


centroids = zeros(ncluster,dimension);
for cc = 1:ncluster
    data_each_cluster = data2cluster(clusterSolution==n_cluster_values(cc),:);
    if size(data_each_cluster,1)>1
        centroids(cc,:) = nanmean(data_each_cluster);
    else
        centroids(cc,:) = data_each_cluster;
    end
end


end