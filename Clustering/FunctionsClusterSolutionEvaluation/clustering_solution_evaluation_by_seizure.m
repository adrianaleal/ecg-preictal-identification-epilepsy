function [] = clustering_solution_evaluation_by_seizure(...
    feat_names2analyse, feat_comb, pats_name, clustering_folder_path)

% Inputs:
% - feat_names2analyse: cell array vector containing the names of the
%   features in the feature dataset
% - feat_comb: double indicating the feature combination under analysis,
%   e.g., if feat_comb = 3 a three-by-three feature combination is being
%   analysed
% - pats_name: cell array containing the indexes of the patients under
%   analysis
% - clustering_folder_path: char path to the folder containing the result
%   of the clustering task


% Outputs:
% - the values of the clustering solution evaluation step are stored in a
%   structure for each seizure


plotComb = 0;
feat_combs = nchoosek(1:numel(feat_names2analyse),feat_comb);
n_comb = size(feat_combs,1);

% ************************************************************************
k_cluster2test_vec = 4;
[clustering_methods, ind_method_dbscan, n_clustering_methods, ...
    k_cluster2test_vec, n_epsilon4dbscan] = initializeClusteringMethods( ...
    k_cluster2test_vec);
% ************************************************************************

methods_vec = [1 4 7:10 11];
n_clustering_methods = numel(methods_vec)


for pp = 1:numel(pats_name)
    
    disp(['Patient' pats_name{pp}])
    
    patFolder2save = ['pat_' pats_name{pp}];
    folderPath = fullfile(clustering_folder_path, patFolder2save);
    filesInFolder = dir(folderPath);
    filesInFolder(ismember({filesInFolder.name}, {'.', '..'})) = [];
    n_seizures_pat = numel(filesInFolder);
    % n_seizures_pat2 = sum(strcmp(names_struct_separated(:,2),file_name));
    
    
    % figure(3)
    for ss = 1:n_seizures_pat % goes through each seizures in each patient
        
        disp(['Seizure ' num2str(ss)])
        
        seizFolder2save = ['Seizure' num2str(ss)];
        
        seizFolderPath = fullfile(folderPath, seizFolder2save);
        
        %% get data2cluster
        % load(fullfile(seizFolderPath,'save_clustering_data'))
        save_clustering_data = parload(seizFolderPath, 'save_clustering_data');
        
        
        % filesInseizFolderPath = dir(seizFolderPath);
        % filesInseizFolderPathName = vertcat({filesInseizFolderPath.name})';
        % if any(strcmp(filesInseizFolderPathName,'results_feat_eval.mat'))
        %   results_feat_eval = parload(seizFolderPath, 'results_feat_eval');
        % else
        clustering_evaluation_measures = {'OD', 'ICV', 'C', 'DI', ...
            'CS', 'SI', 'continuous_cluster_mat', ...
            'duration_cluster_near_seiz_mat', 'nclusters', ...
            'noisy_clusters', 'samples_smaller_cluster'};
        for tt = 1:n_clustering_methods
            for ee = 1:numel(clustering_evaluation_measures)
                results_feat_eval.(clustering_evaluation_measures{ee}).(clustering_methods{methods_vec(tt)}) = NaN(n_comb,1);
            end
        end
        %  end
        
        seiz_data = save_clustering_data.data2cluster;
        % if plotComb
        n_wins = size(seiz_data,2);
        % end
        
        for ff = 1:n_comb % goes through the number of feature combinations
            % compTimePerFeat = tic;
            % disp(['Feature combination ' num2str(ff) ': ' strjoin(feat_names2analyse(feat_combs(ff,:)),' ')])
            % cluster validity indexes:
            
            
            %%
            
            
            seiz_data_feat_comb = squeeze(seiz_data(ff,:,:));
            % normalize data
            % The reason is that normalization gives the same importance to all
            % the variables. The standard example is considering age (in year)
            % and height (in cm). The age may range in [18 50], while the
            % height may range in [130 180]. If you use the classical Euclidean
            % distance, the height will have dispoportionately more importance
            % in its computation with respect to the age.
            feat_data_norm = (seiz_data_feat_comb-min(seiz_data_feat_comb))./(max(seiz_data_feat_comb)-min(seiz_data_feat_comb));
            ind_NaN = any(isnan(feat_data_norm),2);
            feat_data_norm_no_nan = feat_data_norm(~ind_NaN,:);
            
            % compute distances between feature values:
            distances_matrix_feat_data = pdist2(feat_data_norm_no_nan,feat_data_norm_no_nan);
            [nn_matrix,~] = nearest_neighbours(distances_matrix_feat_data);
            
            
            for vv = 5:n_clustering_methods
                % disp(['Clustering Method: ' clustering_methods{vv}])
                
                
                clusteringSolution = save_clustering_data.clustering_solutions.(clustering_methods{methods_vec(vv)})(ff,:);
                % centroids = save_clustering_data.clustering_prototypes.(clustering_methods{methods_vec(vv)}){ff,:};
                
                clusteringSolutionNoNaN = clusteringSolution(~isnan(clusteringSolution));
                % n_cluster_values = unique(clusteringSolutionNoNaN);
                
                
                [centroids, n_cluster_values, ~] = getClusterCentroids(feat_data_norm, ...
                    clusteringSolution);
                
                % first condition: do not analyse clusters with noise
                % (e.g. n_cluster_values = [0 1 2])
                % second condition: do not analyse clusters with one value
                % and NaN values
                
                % BEFORE: size(centroids,1)>1 --> analyse only clustering
                % solutions with more than one cluster
                % AFTER: size(centroids,1)==2 only for clustering solutions
                % with two clusters
                
                if 1 %~any(n_cluster_values==0) && size(centroids,1)>1 % there were no DBSCAN solutions (????)
                    % if the clustering solutions had zeros or didn't have
                    % more than one cluster, no analysis was performed
                    
                    results_feat_eval.nclusters.(clustering_methods{methods_vec(vv)})(ff) = numel(n_cluster_values);
                    results_feat_eval.noisy_clusters.(clustering_methods{methods_vec(vv)})(ff) = any(n_cluster_values==0);
                    frequency_clusters = histc(clusteringSolution, n_cluster_values);
                    results_feat_eval.samples_smaller_cluster.(clustering_methods{methods_vec(vv)})(ff) = min(frequency_clusters);
                    
                    % count_eval_methods = 0;
                    
                    %% Silhouette index
                    
                    % tic
                    silh = silhouette(feat_data_norm_no_nan, clusteringSolutionNoNaN);
                    results_feat_eval.SI.(clustering_methods{methods_vec(vv)})(ff) = mean(silh);
                    % count_eval_methods = count_eval_methods+2;
                    % get_evaluation_seizure(ss,count_eval_methods-1:count_eval_methods) = silh;
                    silh = [];
                    % disp(['Time to compute Silhouette Index: ' num2str(toc)])
                    
                    %% OverallDeviation
                    % tic
                    [overallDeviation, intraclusterVariance] = compactness(feat_data_norm_no_nan, clusteringSolutionNoNaN, []);
                    results_feat_eval.OD.(clustering_methods{methods_vec(vv)})(ff) = overallDeviation;
                    % count_eval_methods = count_eval_methods+2;
                    % get_evaluation_seizure(ss,count_eval_methods-1:count_eval_methods) = [overallDeviation, intraclusterVariance];
                    overallDeviation = [];
                    results_feat_eval.ICV.(clustering_methods{methods_vec(vv)})(ff) = intraclusterVariance;
                    intraclusterVariance = [];
                    % disp(['Time to compute Overall Deviation and ICV: ' num2str(toc)])
                    
                    %% Connectivity
                    % tic
                    C = connectivity(clusteringSolutionNoNaN, nn_matrix, 7, []);
                    results_feat_eval.C.(clustering_methods{methods_vec(vv)})(ff) = C;
                    nn_matrix = [];
                    C = [];
                    % disp(['Time to compute Connectivity: ' num2str(toc)])
                    % count_eval_methods = count_eval_methods+1;
                    % get_evaluation_seizure(ss,count_eval_methods) = C;
                    
                    %% Dunn Index
                    % tic
                    DI = dunns(clusteringSolutionNoNaN, distances_matrix_feat_data, []);
                    results_feat_eval.DI.(clustering_methods{methods_vec(vv)})(ff) = DI;
                    % if results_feat_eval.DI.(clustering_methods{vv})(ff)>=0.13
                    %     pause
                    % end
                    % disp(['Time to compute Dunn Index: ' num2str(toc)])
                    % count_eval_methods = count_eval_methods+1;
                    % get_evaluation_seizure(ss,count_eval_methods) = DI;
                    DI = [];
                    
                    %% Cluster Separation
                    % tic
                    CS = clusterSeparation(feat_data_norm_no_nan, ...
                        clusteringSolutionNoNaN, {centroids});
                    results_feat_eval.CS.(clustering_methods{methods_vec(vv)})(ff) = CS;
                    
                    % disp(['Time to compute Cluster Separation: ' num2str(toc)])
                    % count_eval_methods = count_eval_methods+1;
                    % get_evaluation_seizure(ss,count_eval_methods) = CS;
                    CS = [];
                    
                    %% evaluate time continuity
                    % tic
                    [continuous_cluster_mat, ...
                        duration_cluster_near_seiz_mat, ...
                        distances_matrix_cluster_solution] = ...
                        findContinuousClusters(clusteringSolutionNoNaN, ...
                        clusteringSolution, n_wins, n_cluster_values);
                    results_feat_eval.continuous_cluster_mat.(clustering_methods{vv})(ff) = continuous_cluster_mat;
                    results_feat_eval.duration_cluster_near_seiz_mat.(clustering_methods{vv})(ff) = duration_cluster_near_seiz_mat;
                    continuous_cluster_mat = [];
                    duration_cluster_near_seiz_mat = [];
                    continuous_clusters = [];
                    % disp(['Time to compute Continuous Clusters: ' num2str(toc)])
                    %% plot Clustering
                    
                    if plotComb
                        time = linspace(240,5,n_wins);
                        timeNoNaN = time(~ind_NaN);
                        plotClusterSolution(feat_data_norm_no_nan,clusteringSolutionNoNaN, ...
                            centroids, distances_matrix_cluster_solution, ...
                            n_cluster_values, numel(n_cluster_values), 3, ...
                            0, feat_names2analyse(feat_combs(ff,:)), timeNoNaN);
                        
                        % export_fig(['seiz_' num2str(ss) '_features_' ...
                        % strjoin(feat_names2analyse(comb(ff,:)),'_') ...
                        % '_clustering_method_' clustering_methods{vv} '_CC_' ...
                        % strrep(num2str(continuous_cluster_mat(ss,vv)),'.','') '.pdf'], '-pdf','-transparent')
                    end
                    % distances_matrix_cluster_solution = [];
                    % clusteringSolutionNoNaN = [];
                    
                end
            end
            
            distances_matrix_feat_data = [];
            feat_data_norm = [];
            feat_data_norm_no_nan = [];
            seiz_data_feat_comb = [];
            
            % disp(['Time to compute for one feature and all clustering methods: ' num2str(toc(compTimePerFeat))])
        end
        
        
        % save_clustering_data = [];
        % feat_data = [];
        
        
        %% see if a folder with the previous name exists and if not build
        % that folder in the specified path
        
        if ~exist(seizFolderPath, 'dir')
            mkdir(seizFolderPath)
        end
        
        parsave(seizFolderPath, {results_feat_eval}, 'results_feat_eval')
        
        % if the parsave is not being used the following line of code can
        % replace parsave
        % save(fullfile(seizFolderPath,'results_feat_comb.mat'),'results_feat_comb')
        results_feat_eval = [];
        
        
        
    end
end

end
