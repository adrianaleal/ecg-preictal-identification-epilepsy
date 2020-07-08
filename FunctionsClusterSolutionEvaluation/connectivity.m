function [connectedness, connection_matrix] = connectivity(clusteringSolutions, ...
    nn_matrix, L, n_clustersInClusterSolution)


% Inputs:
% - nn_matrix: nearest neighbours matrix
% - clusteringSolutions: a matrix which rows correspond to the number of
%   samples/observations and columns to the number of clustering solutions
% - L: number of nearest neighbours (default is 10)
% - n_clustersInClusterSolution: number of clusters in each cluster
%   solution, if not provided it is computed in this function


% Outputs:
% - connectedness: value of connectivity for each clustering solution
% - connection_matrix:


% Interpretation:
% - Validity index based on cluster labels.
% - Useful to identify long-shaped clusters.
% - The lowest the cluster connectedness, the more connected are the
%   obtained clusters.


[n_row, n_col] = size(clusteringSolutions);
if n_col>n_row
    clusteringSolutions = clusteringSolutions';
end

if isempty(L)
    L = 10;
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

% initialize connectedness for each clustering solution:
connectedness = zeros(n_clusterSolutions,1);


for ss = 1:n_clusterSolutions % go through each clustering solution
    
    clusteringSolution = clusteringSolutions(:,ss);
    if isempty(n_clustersInClusterSolution)
        n_clusters = numel(unique(clusteringSolution));
    else
        n_clusters = n_clustersInClusterSolution(ss);
    end
    
    connection_matrix = zeros(length(nn_matrix),L);
    
    for cc = 1:length(nn_matrix) % go through each samples/observation
        for ind_j_in_nn_matrix = 2:L+1
            count = 0;
            i_val = clusteringSolution(cc);
            j_val = clusteringSolution(nn_matrix(ind_j_in_nn_matrix,cc));
            for kk = 1:n_clusters
                if i_val==kk && j_val==kk
                    count = count+1;
                end
            end
            if count==0
                connection_matrix(cc,ind_j_in_nn_matrix-1) = 1/(ind_j_in_nn_matrix-1);
            end
        end
    end
    
    connectedness(ss) = sum(sum(connection_matrix));
    
end


end