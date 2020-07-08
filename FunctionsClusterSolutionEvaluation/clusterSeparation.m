function CS = clusterSeparation(data2cluster, clusteringSolutions, ...
    clusterSolutionPrototype)


% Inputs:
% - data2cluster: matrix of data to which clustering was performed, the
%   rows corresponding to samples/observations and the columns to the
%   features considered
% - clusteringSolutions: a matrix which rows correspond to the number of
%   samples/observations and columns to the number of clustering solutions
% - clusterSolutionPrototype: a matrix of 


% Outputs:
% - CS: value of cluster separation for each clustering solution 


% Interpretation:
% - Validity index based on cluster prototypes.
% - Quantifies the intercluster distance by averaging the distance among  
%   the clusters centres
% - Well-separated clusters are characterized by high values of CS.


[n_row, n_col] = size(clusteringSolutions);
if n_col>n_row
    clusteringSolutions = clusteringSolutions';
end


n_clusterSolutions = size(clusteringSolutions,2);


% initialize CS for each clustering solution:
CS = zeros(n_clusterSolutions,1);


for ss = 1:n_clusterSolutions % go through each clustering solution
    
    if isempty(clusterSolutionPrototype{ss})
        clusteringSolution = clusteringSolutions(:,ss);
        cluster_val = unique(clusteringSolution);
        n_clusters = numel(cluster_val);
        clusterSolutionPrototype = zeros(n_clusters,size(data2cluster,2));
        for cc = 1:n_clusters
            cluster = data2cluster(clusteringSolution==cluster_val(cc),:);
            if size(cluster,1)>1
                clusterSolutionPrototype(cc,:) = mean(cluster);
            else
                clusterSolutionPrototype(cc,:) = cluster;
            end
        end
    else
        clusterSolutionPrototype = clusterSolutionPrototype{ss};
        n_clusters = size(clusterSolutionPrototype,1);
    end
    
    if n_clusters>1
        D_mat = dist(clusterSolutionPrototype');
        upper_values = triu(ones(n_clusters),1);
%         [ind_row_upper_values, ind_col_upper_values] = find(upper_values);
        getDist = D_mat(logical(upper_values));
        CS(ss) = 2/(n_clusters*(n_clusters-1))*sum(getDist);
    else
        CS(ss) = NaN;
    end
end

end