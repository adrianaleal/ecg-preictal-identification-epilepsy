function [overallDeviation, intraClusterVariance] = compactness(data2cluster, ...
    clusteringSolutions, n_clustersInClusterSolution)


% Inputs:
% - data2cluster: matrix of data to which clustering was performed, the
%   rows corresponding to samples/observations and the columns to the
%   features considered
% - clusteringSolutions: a matrix which rows correspond to the number of
%   samples/observations and columns to the number of clustering solutions
% - n_clustersInClusterSolution: number of clusters in each cluster 
%   solution, if not provided it is computed in this function

% Outputs:
% - overallDeviation: value of compactness for each clustering solution
% - intraClusterVariance: value of compactness for each clustering solution


% Interpretation:
% - Validity index based on cluster prototypes.
% - Used to capture intracluster variation.
% - Useful when clustering methods, such as k-means, return round shaped 
%   cluster solutions.
% - The lowest the overallDeviation of the clustering solution, the more 
%   compact are the obtained clusters and the highest is the number of 
%   clusters.


[n_row, n_col] = size(clusteringSolutions);
if n_col>n_row
    clusteringSolutions = clusteringSolutions';
end

[n_row, n_col] = size(data2cluster);
if n_col>n_row
    data2cluster = data2cluster';
end

n_clusterSolutions = size(clusteringSolutions,2);

% check if the number of clusters is the same as the number of clustering
% solutions:
if ~isempty(n_clustersInClusterSolution)
    if numel(n_clustersInClusterSolution)~=n_clusterSolutions
        disp(['n_clustersInClusterSolution must have the same size as ' ...
            'the number of the clustering solution'])
        return
    end
end

% initialize overallDeviation and intraClusterVariance for each clustering 
% solution:
overallDeviation = zeros(n_clusterSolutions,1);
intraClusterVariance = overallDeviation;


for ss = 1:n_clusterSolutions % go through each clustering solution
    
    clusteringSolution = clusteringSolutions(:,ss);
    if isempty(n_clustersInClusterSolution)
        n_clusters = numel(unique(clusteringSolution));
    else
        n_clusters = n_clustersInClusterSolution(ss);
    end
    
    deviations = zeros(1,n_clusters);
    variances = deviations;

    for kk = 1:n_clusters
        cluster = data2cluster(clusteringSolution==kk,:);
        if size(cluster,1)>1
            centroid = mean(cluster);
        else
            centroid = cluster;
        end
        deviations(kk) = sum(pdist2(cluster,centroid));
        variances(kk) = sum((pdist2(cluster,centroid)).^2);
    end

    overallDeviation(ss) = sum(deviations);

    ninstances = size(data2cluster,2);
    intraClusterVariance(ss) = (1/ninstances)*sum(variances);

end

end