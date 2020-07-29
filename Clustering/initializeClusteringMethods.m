function [clustering_methods, ind_method_dbscan, n_clustering_methods, ...
    k_cluster2test_vec, n_epsilon4dbscan] = initializeClusteringMethods(k_cluster2test_vec)

% INDICATE THE CLUSTERING METHODS; 

% A DEFINIR: *************************************************************
ind_pat = 7;
k_cluster2test_vec = 2:k_cluster2test_vec;
n_k_cluster2test = numel(k_cluster2test_vec);
n_epsilon4dbscan = 4;


% K-means clustering:
k_means_clustering = cellstr(strcat(repmat('kmeans_k',n_k_cluster2test,1), ...
    num2str(k_cluster2test_vec')));

% Agglomerative hierarchical clustering:
agglo_hier_clustering = cellstr(strcat(repmat('agglo_hier_k',n_k_cluster2test,1), ...
    num2str(k_cluster2test_vec')));

% DBSCAN
% 4 different epsilons/distances are tested in dbscan
dbscan_clustering = cellstr(strcat(repmat('dbscan_d',n_epsilon4dbscan,1), ...
    num2str((1:n_epsilon4dbscan)')));

% Gaussian mixture model clustering:
gmm_clustering = cellstr(strcat(repmat('gmm_k',n_k_cluster2test,1), ...
    num2str(k_cluster2test_vec')));

% clustering_methods = [k_means_clustering; agglo_hier_clustering; ...
%     dbscan_clustering]; %'moea'


clustering_methods = [k_means_clustering; agglo_hier_clustering; ...
    dbscan_clustering; gmm_clustering]; % ; 'moea'

IndexC = strfind(clustering_methods,'dbscan');
ind_method_dbscan = find(not(cellfun('isempty',IndexC)));


% clustering_methods = cellstr(strcat(repmat('dbscan_d',n_epsilon4dbscan,1), ...
%     num2str((1:n_epsilon4dbscan)')));

n_clustering_methods = numel(clustering_methods);


end