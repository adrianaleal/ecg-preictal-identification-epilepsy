% Example of how to run clustering methods on the supplied feature data

fclose all; clear; close all; clc;




%% Load extracted features for all 4 seizures of patient 5 (21902)
% (that were used in Figure 2 for specific 3-by-3 feature combinations)

% defined number of features to combine(e.g. 3-by-3 or 2-by-2)
feat_comb = 3;


feature_dataset = load('feature_dataset_240min_before_seizure_pat_21902.mat');
feature_dataset = feature_dataset.feature_dataset_240min_before_seizure_pat_21902;

% get feature names
feat_names = fieldnames(feature_dataset);
feat_names = feat_names(3:end);


% indexes of features for each feature combination:
feat_combs = nchoosek(1:numel(feat_names), feat_comb);

% choose the feature combinations for each seizure
feat_comb_seizures = [{'SD2', 'LF_POWER', 'RRMin'}; ...
                      {'RQA_Lmax', 'SD1', 'LFtoHF'}; ...
                      {'DFA_alpha1', 'ApEn', 'NN50'}; ...
                      {'RQA_L', 'HF_NORM', 'RRMax'}]

% choose the clustering method for each seizure
clustering_methods = {'dbscan_d3', 'dbscan_d3', 'dbscan_d3', 'gmm_k'}
% options: 'kmeans_k2', 'agglo_hier_k2', 

epsilon_vec = [0.1 0.15 0.2 0.3];
ncluster = 2; % number of clusters for kmeans and agglomerative
plotClustering = 1;
opts_KMEANS = statset('Display','final');
 
figure()               
n_seiz = size(feature_dataset.NN50,1);                 
for ss = 1:n_seiz
    
    feat_names_seiz = feat_comb_seizures(ss,:)';
    ind_feat_comb = find(sum(ismember(feat_names(feat_combs), feat_comb_seizures(ss,:)), 2)==3);
    
    for ff = 1:numel(feat_names_seiz)
        feat_data_seiz.(feat_names_seiz{ff}) = feature_dataset.(feat_names_seiz{ff})(ss,:);
    end
    feat_data_seiz.segment_size_seconds = feature_dataset.segment_size_seconds(ss,:);
    feat_data_seiz.percent_detected_noise = feature_dataset.percent_detected_noise(ss,:);
    
    feat_data_seiz = editInvalidWindows(feat_names_seiz, feat_data_seiz);
    
    feat_data = [vertcat(feat_data_seiz.(feat_names_seiz{1})); ...
        vertcat(feat_data_seiz.(feat_names_seiz{2})); ...
        vertcat(feat_data_seiz.(feat_names_seiz{3}))]';
    
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
    
    
    
    
    %% Run K-means clustering algorithm for 2 clusters
    % ncluster2kmeans = 4;
    kmeans_name = 'kmeans_k';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss),kmeans_name)))

        [clusteringSolutionKmeans, centroids] = kmeans(feat_data_norm_no_NaN, ...
            ncluster, 'Distance', 'sqeuclidean', 'Options', opts_KMEANS, ...
            'MaxIter', 10);
        
        % add the cluster values to the proper non NaN positions:
        clusteringSolution = NaN(size(ind_NaN));
        clusteringSolution(~ind_NaN) = clusteringSolutionKmeans;
    end
    
    %% Run agglomerative hierarchical clustering
    
    agglo_name = 'agglo_hier_k';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss),agglo_name)))
        
        Z = linkage(feat_data_norm_no_NaN,'ward','euclidean');
        
        clusteringSolutionAgglo = cluster(Z, 'maxclust', ncluster);
        
        % add the cluster values to the proper non NaN positions:
        clusteringSolution = NaN(size(ind_NaN));
        clusteringSolution(~ind_NaN) = clusteringSolutionAgglo; 
        
    end
    
    %% Run DBSCAN Clustering Algorithm on the features intervals
    dbscan_name = 'dbscan_d';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss),dbscan_name)))
        
        epsilon = epsilon_vec(str2double(clustering_methods{ss}(end)));
        [centroids, clusteringSolutionDBSCAN] = getDBSCANClustering( ...
            dimension, feat_data_norm_no_NaN, epsilon, feat_names, ...
            time_min, plotClustering);
        
        % add the cluster values to the proper non NaN positions:
        clusteringSolution = NaN(size(ind_NaN));
        clusteringSolution (~ind_NaN) = clusteringSolutionDBSCAN;
        
    end
    
    %% Run gaussian mixture model clustering
    
    gmm_name = 'gmm_k';
    if any(~cellfun(@isempty,strfind(clustering_methods(ss), gmm_name)))
        
        [centroids, clusteringSolution] = getGMMClustering( ...
            feat_data_norm_no_NaN, ind_NaN, plotClustering, feat_names);
        
    end
    
    if plotClustering==1
        subplot(2,2,ss)
        
        n_cluster_values = unique(clusteringSolution);
        % as unique function considers NaN values as unique the following 
        % code was added to remove the NaN:
        n_cluster_values = (n_cluster_values(~isnan(n_cluster_values)));
        
        [centroids, n_cluster_values, dimension] = ...
                getClusterCentroids(feat_data_norm, clusteringSolution);
            
        [percent_missing_values, frequency_smaller_cluster] = plotFigureClusteringResults( ...
            ss, clusteringSolution, n_cluster_values, numel(n_cluster_values), ...
            n_seiz, selected_feat, feat_names2analyse, patFolderPath, ...
            name2save, figure_title, feat_data_norm, ...
            time_vec_original{ind_seizures_pat(ss)}, ...
            percent_noise_original{ind_seizures_pat(ss)}, title_name);
        
        
        if any(~cellfun(@isempty,strfind(clustering_methods(ss),kmeans_name)))
            clust_method = ['KM, k = ' num2str(ncluster)];
        elseif any(~cellfun(@isempty,strfind(clustering_methods(ss),agglo_name)))
            clust_method = ['AH, k = ' num2str(ncluster)];
        elseif any(~cellfun(@isempty,strfind(clustering_methods(ss),dbscan_name)))
            clust_method = ['DBSCAN_' dbscan_name(end)];
        elseif any(~cellfun(@isempty,strfind(clustering_methods(ss), gmm_name)))
            clust_method = 'GMM';
        end
        
        title(['seizure ' num2str(ss) clust_method])
        
    end
        
end