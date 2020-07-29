function DI = dunns(clusteringSolutions, dissimilarityMat, ...
    n_clustersInClusterSolution)


% Inputs:
% - dissimilarityMat: dissimilarity matrix
% - clusteringSolutions: a matrix which rows correspond to the number of
%   samples/observations and columns to the number of clustering solutions
% - n_clustersInClusterSolution: number of clusters in each cluster 
%   solution, if not provided it is computed in this function


% Outputs:
% - DI: value of Dunn's Index for each clustering solution


% Interpretation:
% - Validity index based on cluster labels.
% - If the cluster solution represents compact and well-separated clusters,
%   a high value of Dunn Index will be obtained.
% - DI = 0 is obtained when the clustering solution constains one cluster


[n_row, n_col] = size(clusteringSolutions);
if n_col>n_row
    clusteringSolutions = clusteringSolutions';
end


n_clusterSolutions = size(clusteringSolutions,2);

% check if the number of clusters is the same as the number of clustering
% solutions:
if ~isempty(n_clustersInClusterSolution)
    if numel(n_clustersInClusterSolution)~= n_clusterSolutions
        disp('n_clustersInClusterSolution must have the same size as the number of the clustering solution')
        return
    end
end

% initialize DI for each clustering solution:
DI = zeros(n_clusterSolutions,1);


for ss = 1:n_clusterSolutions % go through each clustering solution
    
    clusteringSolution = clusteringSolutions(:,ss);
    if isempty(n_clustersInClusterSolution)
        n_clusters = numel(unique(clusteringSolution));
    else
        n_clusters = n_clustersInClusterSolution(ss);
    end
    
    if n_clusters>1
        
        denominator = [];
        
        for kk = 1:n_clusters
            indi = find(clusteringSolution==kk);
            indj = find(clusteringSolution~=kk);
            temp = dissimilarityMat(indi,indj);
            denominator = [denominator; temp(:)];
        end
        
        num = min(min(denominator));
        neg_obs = zeros(size(dissimilarityMat,1),size(dissimilarityMat,2));
        
        for kk = 1:n_clusters
            indxs = find(clusteringSolution==kk);
            neg_obs(indxs,indxs) = 1;
        end
        
        dem = neg_obs.*dissimilarityMat;
        dem = max(max(dem));
        
        DI(ss) = num/dem;
        
    else
        DI(ss) = NaN;
    end
end

end
