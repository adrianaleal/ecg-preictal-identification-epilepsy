function [centroids_matrix, clustering_solutions_matrix] = getDBSCANClustering(dimension, ...
    feat_data_norm, feat_names, time_min, plotClustering, n_epsilon4dbscan)


centroids_matrix = cell(n_epsilon4dbscan,1);
clustering_solutions_matrix = cell(n_epsilon4dbscan,1);


ind_NaN = any(isnan(feat_data_norm),2);
feat_data_norm_no_NaN = feat_data_norm(~ind_NaN,:);

%% Run DBSCAN Clustering Algorithm on the features intervals

% since data is normalized it is possible do define different
% epsilon measures in the range [0 1], without needing to find an
% adequate epsilon
if dimension==2
    % For two-dimensional data: use default value of minPts=4
    % (Ester et al., 1996)
    MinPts_DBSCAN = 4;
elseif dimension>2
    % For more than 2 dimensions: minPts=2*dim
    % (Sander et al., 1998)
    MinPts_DBSCAN = 2*dimension;
end

count_cluster_methods = 1;

count_dbscan = 1;
stop = 1;
% eps: Two points are considered neighbors if the distance
% between the two points is below the threshold epsilon.
% epsilon_DBSCAN_start = 0.05;
%
% firstPassInUnique = 0;
% save_eps = [];
% distances_matrix = pdist2(feat_data_norm_no_NaN,feat_data_norm_no_NaN);
% [indexes_sorted_distances, sorted_distances] = nearest_neighbours(distances_matrix);
% mean_sorted_distances = mean(sorted_distances,2);
% diff_mean_sorted_distances = diff(mean_sorted_distances);
% [sorted_diff_mean_sorted_distances, ind_sorted] = sort(diff_mean_sorted_distances);
% mean_sorted_distances = mean_sorted_distances(ind_sorted<(ind_sorted(end)-3))

% uma forma
thresholds = [0.1 0.15 0.2 0.3];
eps_vec = thresholds;
% inds = zeros(numel(thresholds),1);
% for tt = 1:numel(thresholds)
%     inds_vec = find(mean_sorted_distances>thresholds(tt));
%     inds(tt) = inds_vec(1);
% end

% figure(31)
% plot(mean_sorted_distances())
% hold on
% plot(diff_mean_sorted_distances)
% plot(inds,mean_sorted_distances(inds),'r*')
% inds2plot = ind_sorted(end-20:end)+3;
% inds2plot_final = inds2plot(inds2plot<numel(mean_sorted_distances));
% plot(inds2plot_final, mean_sorted_distances(inds2plot_final), 'g*')
% hold off
% axis tight

while stop
    % epsilon = epsilon_vec(ndist);
    epsilon_DBSCAN_start = eps_vec(count_dbscan);
    
    plotFigure = 0;
    % the cluster solution of DBSCAN has zeros corresponding to
    % noisy points
    clusteringSolutionDBSCAN = DBSCAN(feat_data_norm_no_NaN, ...
        epsilon_DBSCAN_start, MinPts_DBSCAN, plotFigure);
    
    
    %     ind_zeros = clusteringSolutionDBSCAN==0;
    %     check_NaNs = isequal(ind_zeros, ind_NaN);
    %
    %     if ~isempty(ind_NaN) % when the vector of data with NaN
    %         is given to the cluster method
    %         clusteringSolutionDBSCAN(ind_NaN) = NaN;
    %     end
    if plotClustering==1
        [centroids, n_cluster_values, dimension] = getClusterCentroids(feat_data_norm_no_NaN, ...
            clusteringSolutionDBSCAN);
        
        distances = [];
        ncluster = numel(n_cluster_values);
        figure(67)
        plotClusterSolution(feat_data_norm_no_NaN, ...
            clusteringSolutionDBSCAN, centroids, ...
            distances, n_cluster_values, ncluster, ...
            dimension, 1, feat_names, time_min, [], [])
        
    else
        [centroids, n_cluster_values, dimension] = getClusterCentroids(feat_data_norm_no_NaN, ...
            clusteringSolutionDBSCAN);
    end
    
    frequency_clusters = histc(clusteringSolutionDBSCAN, n_cluster_values);
    
    
%     if numel(unique(clusteringSolutionDBSCAN))==1 || min(frequency_clusters)<10 || firstPassInUnique>0
%         firstPassInUnique = firstPassInUnique + 1;
%         if firstPassInUnique==1
%             epsilon_DBSCAN_start = save_eps(1);
%         end
%         epsilon_DBSCAN_increment = epsilon_DBSCAN_start/2;
%         epsilon_DBSCAN_start = epsilon_DBSCAN_start-epsilon_DBSCAN_increment;
%     else
%         epsilon_DBSCAN_start = epsilon_DBSCAN_start+epsilon_DBSCAN_increment;
%     end
    
    
    count_dbscan = count_dbscan+1;
    % save clustering solutions and centroids *****************
    % save_clustering_solutions(:,count_cluster_methods) = clusteringSolutionDBSCAN;
    
    centroids_matrix(count_cluster_methods) = {centroids};
    clusteringSolutionDBSCANFinal = NaN(size(ind_NaN));
    clusteringSolutionDBSCANFinal(~ind_NaN) = clusteringSolutionDBSCAN;
    clustering_solutions_matrix(count_cluster_methods) = {clusteringSolutionDBSCANFinal};
    count_cluster_methods = count_cluster_methods+1;
    
    % clustering_solutions.([dbscan_name num2str(count_dbscan)])(ss,:) = clusteringSolutionDBSCAN;
    
    
    if count_dbscan==(n_epsilon4dbscan+1)
        stop = 0;
    end
    
%     figure()
%     PlotClusteringResultDBSCAN(feat_data_norm, clusteringSolutionDBSCANFinal,{'m','x'});
%     title(['DBSCAN Clustering ($\epsilon$ = ' num2str(epsilon_DBSCAN_start) ...
%         ' samples, MinPts = ' num2str(MinPts_DBSCAN) ')']);
%     xlabel(feat_names(1))
%     ylabel(feat_names(2))
%     zlabel(feat_names(3))
%     axis tight
    
end
count_dbscan = [];
stop = [];
epsilon_DBSCAN_start = [];
dbscan_name = [];
clusteringSolutionDBSCAN = [];
clusteringSolutionDBSCANFinal = [];
centroids = [];


end