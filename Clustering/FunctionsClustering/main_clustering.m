% Example of how to run clustering methods on the supplied feature data

fclose all; clear; close all; clc;




%% Load extracted features for all 4 seizures of patient 5 (21902)
% (that were used in Figure 2 for specific 3-by-3 feature combinations)

% defined number of features to combine(e.g. 3-by-3 or 2-by-2)
feat_comb = 3;


feature_dataset = load('feature_dataset_240min_before_seizure_pat_21902.mat');
feature_dataset = feature_dataset.feature_dataset_240min_before_seizure_pat_21902;

% get feature names
feat_names = fieldnames(feature_dataset{1});
feat_names = feat_names(4:end);


% indexes of features for each feature combination:
feat_combs = nchoosek(1:numel(feat_names), feat_comb);

% choose the feature combinations for each seizure
feat_comb_seizures = [{'RRMin', 'LF_POWER', 'SD2'}; ...
    { 'LFtoHF', 'SD1', 'RQA_Lmax'}; ...
    {'NN50', 'ApEn', 'DFA_alpha1'}; ...
    {'RRMax', 'HF_NORM', 'RQA_L'}];

% choose the clustering method for each seizure
clustering_methods = {'dbscan_d3', 'dbscan_d3', 'dbscan_d3', 'gmm_k2'};
% options: 'kmeans_k2', 'agglo_hier_k2', 'dbscan_d1', 'dbscan_d2', 
% 'dbscan_d3', 'dbscan_d4', 'gmm_k2'

epsilon_vec = [0.1 0.15 0.2 0.3]; % epsilon values used for DBSCAN
ncluster = 2; % number of clusters for kmeans, agglomerative hierarchical 
% and gaussian mixture model
plotClusteringSolution = 0; % flag to plot exclusively the clustering 
% solutions
opts_KMEANS = statset('Display','final'); % options to compute kmeans


n_seiz = size(feature_dataset,1);
for ss = 1:n_seiz
    disp(['Seizure ' num2str(ss) ' **************************************'])
    
    feat_names_seiz = feat_comb_seizures(ss,:)';
    ind_feat_comb = find(sum(ismember(feat_names(feat_combs), feat_comb_seizures(ss,:)), 2)==3)
    
    for ff = 1:numel(feat_names_seiz)
        feat_data_seiz.(feat_names_seiz{ff}) = vertcat(feature_dataset{ss}.(feat_names_seiz{ff}));
    end
    feat_data_seiz.segment_size_seconds = vertcat(feature_dataset{ss}.segment_size_seconds);
    feat_data_seiz.percent_detected_noise = vertcat(feature_dataset{ss}.percent_detected_noise);
    
    feat_data_seiz = editInvalidWindows(feat_names_seiz, feat_data_seiz);
    
    feat_data = [vertcat(feat_data_seiz.(feat_names_seiz{1})), ...
        vertcat(feat_data_seiz.(feat_names_seiz{2})), ...
        vertcat(feat_data_seiz.(feat_names_seiz{3}))];
    
    dimension = size(feat_data,2);
    n_wins = length(feat_data);
    time_min = linspace(5,240,n_wins);
    
    %% normalize data
    % The reason is that normalization gives the same importance to all
    % the variables. The standard example is considering age (in year)
    % and height (in cm). The age may range in [18 50], while the
    % height may range in [130 180]. If you use the classical Euclidean
    % distance, the height will have dispoportionately more importance
    % in its computation with respect to the age.
    feat_data_norm = (feat_data-min(feat_data))./(max(feat_data)-min(feat_data));
    
    ind_NaN = any(isnan(feat_data_norm),2);
    feat_data_norm_no_NaN = feat_data_norm(~ind_NaN,:);
    
    
    
    
    %% Run K-means clustering
    
    kmeans_name = 'kmeans_k';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss),kmeans_name)))
        
        [clusteringSolutionNoNaN, centroids] = kmeans(feat_data_norm_no_NaN, ...
            ncluster, 'Distance', 'sqeuclidean', 'Options', opts_KMEANS, ...
            'MaxIter', 10);
        
        % add the cluster values to the proper non NaN positions:
        clusteringSolution = NaN(size(ind_NaN));
        clusteringSolution(~ind_NaN) = clusteringSolutionNoNaN;
        
    end
    
    %% Run agglomerative hierarchical clustering
    
    agglo_name = 'agglo_hier_k';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss),agglo_name)))
        
        Z = linkage(feat_data_norm_no_NaN,'ward','euclidean');
        
        clusteringSolutionNoNaN = cluster(Z, 'maxclust', ncluster);
        
        % add the cluster values to the proper non NaN positions:
        clusteringSolution = NaN(size(ind_NaN));
        clusteringSolution(~ind_NaN) = clusteringSolutionNoNaN;
        
    end
    
    %% Run DBSCAN Clustering
    dbscan_name = 'dbscan_d';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss),dbscan_name)))
        
        epsilon = epsilon_vec(str2double(clustering_methods{ss}(end)));
        clusteringSolutionNoNaN = getDBSCANClustering(dimension, ...
            feat_data_norm_no_NaN, epsilon);
        
        % add the cluster values to the proper non NaN positions:
        clusteringSolution = NaN(size(ind_NaN));
        clusteringSolution (~ind_NaN) = clusteringSolutionNoNaN;
        
    end
    
    %% Run gaussian mixture model clustering
    
    gmm_name = 'gmm_k';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss), gmm_name)))
        
        plotClusteringSolutionGMM = 0; % to plot the four different 
        % clustering solutions that are chosen within GMM
        clusteringSolutionNoNaN = getGMMClustering(feat_data_norm_no_NaN, ...
            ncluster, plotClusteringSolutionGMM);
        
        % add the cluster values to the proper non NaN positions:
        clusteringSolution = NaN(size(ind_NaN));
        clusteringSolution(~ind_NaN) = clusteringSolutionNoNaN;
    end
    
    
    n_cluster_values = unique(clusteringSolution);
    % as unique function considers NaN values as unique the following
    % code was added to remove the NaN:
    n_cluster_values = (n_cluster_values(~isnan(n_cluster_values)));
    
    
    if any(~cellfun(@isempty,strfind(clustering_methods(ss),kmeans_name)))
        clust_method = ['KM, k = ' num2str(ncluster)];
    elseif any(~cellfun(@isempty,strfind(clustering_methods(ss),agglo_name)))
        clust_method = ['AH, k = ' num2str(ncluster)];
    elseif any(~cellfun(@isempty,strfind(clustering_methods(ss),dbscan_name)))
        clust_method = ['DBSCAN_' clustering_methods{ss}(end)];
    elseif any(~cellfun(@isempty,strfind(clustering_methods(ss), gmm_name)))
        clust_method = ['GMM, k = ' num2str(ncluster)];
    end
    
    title_string = ['seizure ' num2str(ss) ', ' clust_method];
    
    n_samples_smaller_cluster = plotFigureClusteringResults( ...
        ss, clusteringSolution, n_cluster_values, n_seiz, ...
        feat_names_seiz, feat_data_norm, title_string)
    
    if plotClusteringSolution==1
        [centroids, n_cluster_values, dimension] = getClusterCentroids( ...
            feat_data_norm_no_NaN, clusteringSolutionNoNaN);
        distances = [];
        figure()
        plotClusterSolution(feat_data_norm_no_NaN, clusteringSolutionNoNaN, ...
            centroids, distances, dimension, 1, feat_names, [], [], []);
    end
    
end