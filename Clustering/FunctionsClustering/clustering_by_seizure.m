function [] = clustering_by_seizure(feat_names2analyse, feat_comb, ...
    pats_name, names_struct_separated, clustering_folder_path, ...
    outer_folder_path) % , time


plotComb = 0;
feat_combs = nchoosek(1:numel(feat_names2analyse),feat_comb);
n_comb = size(feat_combs,1);
% considering only time feature:
% ind_time = strcmp(feat_names2analyse(comb(:,end)),'time');
% comb = comb(ind_time,:);

% *************************************************************************
k_cluster2test_vec = 4;
[clustering_methods, ind_method_dbscan, n_clustering_methods, ...
    k_cluster2test_vec, n_epsilon4dbscan] = initializeClusteringMethods( ...
    k_cluster2test_vec);

% clustering_methods = clustering_methods(11);
% n_clustering_methods = numel(clustering_methods);
% the following code can be uncommented when all methods are considered:
methods_vec = [1 4 7:10 11];
% methods_vec = 1; % to compute the GMM results
k_cluster2test_vec = k_cluster2test_vec(1); % kmeans and AH only for k = 2;
n_clustering_methods = numel(methods_vec);
% *************************************************************************

% see if a folder with the previous name exists and if not build
% that folder in the specified path
if ~exist(clustering_folder_path, 'dir')
    mkdir(clustering_folder_path)
end

% *************************************************************************

% for seizure original data loading:
feat_folder = 'ExtractedFeaturesNoEctopicBeatsOverlap';

current_folder = cd;

plotClustering = 0;
opts_KMEANS = statset('Display','final');

% parfor (ll = 1:size(comb,1), 'debug')
parfor pp = 1:numel(pats_name)
    
    
    % for loading feature original data for each seizure:
    file_name = pats_name{pp};
    disp(['Patient' file_name])
    patFolder2save = ['pat_' file_name];
    seizFolders = dir(fullfile(outer_folder_path,'ResultsEPILEPSIAE', ...
        feat_folder, patFolder2save));
    seizFolders = seizFolders(~ismember({seizFolders.name},{'.','..'}));
    
    
    patFolderPath = fullfile(clustering_folder_path, patFolder2save);
    if ~exist(patFolderPath, 'dir')
        mkdir(patFolderPath)
    end
    
    n_seizures_pat = sum(strcmp(names_struct_separated(:,2),file_name));
    % ind_seizures_pat = find(strcmp(file_name,C2(:,2)));
    
    for ss = 1:n_seizures_pat % goes through the number of seizures
        
        disp(['Seizure ' num2str(ss)])
        
        %% see if a folder with the previous name exists and if not build
        % that folder in the specified path and initialize
        % save_clustering_data structure
        seizFolderPath = fullfile(patFolderPath, ['Seizure' num2str(ss)]);
        if exist(seizFolderPath, 'dir')
            filesInseizFolderPath = dir(seizFolderPath);
            filesInSeizFolderPathName = vertcat({filesInseizFolderPath.name})';
            if any(strcmp(filesInSeizFolderPathName,'save_clustering_data.mat'))
                save_clustering_data = parload(seizFolderPath, 'save_clustering_data');
                % if there's no need to compute then:
                % continue
            else
                save_clustering_data = struct();
            end
        else
            mkdir(seizFolderPath)
            save_clustering_data = struct();
        end
        
        %% load of the original feature data:
        hrv_features = parload(fullfile(seizFolders(ss).folder, seizFolders(ss).name), ...
            [patFolder2save '_seizure' num2str(ss)], 'hrv_features');
        
        % the following load is not allowed in parfor
        % load(fullfile(seizFolders(ss).folder,seizFolders(ss).name, ...
        %     [patFolder2save '_seizure' num2str(ss)]), 'hrv_features')
        
        n_wins = length(hrv_features);
        
        %% initialize variables to save results
        temp_solutions = cell(n_clustering_methods,1);
        temp_prototypes = cell(n_clustering_methods,1);
        for cc = 1:n_clustering_methods
            temp_solutions(cc) = {zeros(n_comb,n_wins)};
            temp_prototypes(cc) = {cell(n_comb,1)};
        end
        
        % the previous instead of the following (not accepted by parfor)
        % for cc = 1:n_clustering_methods
        %     clustering_solutions.(clustering_methods{cc}) = zeros(n_seizures,n_wins);
        %     clustering_prototypes.(clustering_methods{cc}) = cell(n_seizures,1);
        % end
        
        
        data2cluster = zeros(n_comb,n_wins,feat_comb);
        
        %%
        % goes through the number of feature combinations
        for ff = 1:n_comb
            % disp(['Feature combination ' num2str(ff)])
            
            feat_names = feat_names2analyse(feat_combs(ff,:));
            
            
            % save_clustering_solutions = zeros(n_wins,n_clustering_methods);
            % save_clustering_prototypes = cell(1,n_clustering_methods);
            
            
            % load of the original data:
            
            feature_dataset = editInvalidWindows(feat_names, hrv_features);
            feat_data = [vertcat(feature_dataset.(feat_names{1})), ...
                vertcat(feature_dataset.(feat_names{2})), ...
                vertcat(feature_dataset.(feat_names{3}))];
            
            data2cluster(ff,:,:) = feat_data;
            dimension = size(data2cluster,3);
            
            %% normalize data
            % The reason is that normalization gives the same importance to all
            % the variables. The standard example is considering age (in year)
            % and height (in cm). The age may range in [18 50], while the
            % height may range in [130 180]. If you use the classical Euclidean
            % distance, the height will have dispoportionately more importance
            % in its computation with respect to the age.
            feat_data_norm = (feat_data-min(feat_data))./(max(feat_data)-min(feat_data));
            % feat_data_norm = bsxfun(@rdivide, bsxfun(@minus, feat_data, min(feat_data)), ...
            %     (max(feat_data)-min(feat_data)));
            
            
            ind_NaN = any(isnan(feat_data_norm),2);
            feat_data_norm_no_NaN = feat_data_norm(~ind_NaN,:);
            feat_data = [];
            % K-means clustering is "isotropic" in all directions of space and
            % therefore tends to produce more or less round (rather than
            % elongated) clusters. In this situation leaving variances unequal
            % is equivalent to putting more weight on variables with smaller
            % variance.
            % Source: https://stackoverflow.com/questions/15777201/why-vector-normalization-can-improve-the-accuracy-of-clustering-and-classificati
            
            %%
            
            time_min = linspace(5,240,n_wins);
            if plotComb
                %             time_data = time_dataset(ss,:);
                
                figure(1)
                %         set(gcf,'units','normalized','outerposition',[0 0 1 0.3])
                if feat_comb==2
                    scatter(feat_data_norm(:,1),feat_data_norm(:,2),[],time_min,'*')
                    
                elseif feat_comb==3
                    scatter3(feat_data_norm(:,1),feat_data_norm(:,2),feat_data_norm(:,3),[],time_min,'*')
                    zlabel(feat_names2analyse{feat_combs(ff,3)})
                end
                xlabel(feat_names2analyse{feat_combs(ff,1)})
                ylabel(feat_names2analyse{feat_combs(ff,2)})
                title(['Seizure ' num2str(ss)])
                c=colorbar;
                ylabel(c,'Time before seizure (min)','Interpreter','latex','Fontsize',14)
                pause
                % saveas(gcf,char(strcat(filename_dat{ii}(1:end-4), '_', name_feat2, '_', name_feat1, '.emf')))
                % close all
            end
            
            % perform clustering on feature combination
            
            
            %% Run K-means clustering algorithm for 2 clusters
            count_cluster_methods = 1;
            % ncluster2kmeans = 4;
            kmeans_name = 'kmeans_k';
            if any(~cellfun(@isempty,strfind(clustering_methods,kmeans_name)))
                for ncluster = k_cluster2test_vec
                    % disp([kmeans_name num2str(ncluster)])
                    
                    [clusteringSolutionKmeans, centroids] = kmeans(feat_data_norm_no_NaN, ...
                        ncluster,'Distance','sqeuclidean','Options',opts_KMEANS,'MaxIter',10);
                    if plotClustering==1
                        distances = [];
                        n_cluster_values = unique(clusteringSolutionKmeans);
                        % as unique function considers NaN values as unique the following code was
                        
                        % added to remove the NaN:
                        n_cluster_values = (n_cluster_values(~isnan(n_cluster_values)));
                        
                        plotClusterSolution(feat_data_norm_no_NaN, clusteringSolutionKmeans, ...
                            centroids, distances, n_cluster_values, ncluster, ...
                            dimension, 1, feat_names, time_min, [], []);
                    end
                    
                    % save clustering solutions and centroids *****************
                    
                    % save_clustering_solutions(:,count_cluster_methods) = clusteringSolutionKmeans;
                    
                    % if the input is feat_data_norm instead of feat_data_norm_no_NaN
                    % the following must be confirmed:
                    % ind_NaN_cluster_solution = any(isnan(clusteringSolutionKmeans),2);
                    % check_cluster_solution = isequal(ind_NaN, ind_NaN_cluster_solution)
                    
                    % add the cluster values to the proper non NaN positions:
                    clusteringSolutionKmeansFinal = NaN(size(ind_NaN));
                    clusteringSolutionKmeansFinal(~ind_NaN) = clusteringSolutionKmeans;
                    
                    
                    temp_solutions{count_cluster_methods}(ff,:) = clusteringSolutionKmeansFinal;
                    temp_prototypes{count_cluster_methods}(ff) = {centroids};
                    count_cluster_methods = count_cluster_methods+1;
                    
                    % the previous instead of the following (not accepted by parfor)
                    % clustering_solutions.(kmeans_name)(ss,:) = clusteringSolutionKmeans;
                    % clustering_prototypes.(kmeans_name)(ss) = {centroids};
                    
                end
                
            end
            kmeans_name = [];
            centroids = [];
            clusteringSolutionKmeans = [];
            clusteringSolutionKmeansFinal = [];
            
            %% Run agglomerative hierarchical clustering
            
            agglo_name = 'agglo_hier_k';
            if any(~cellfun(@isempty,strfind(clustering_methods,agglo_name)))
                % if ~any(find(isnan(data2cluster)))
                %   Y = pdist(feat_data_norm,'euclidean');
                % else
                % % Y = pdist(feat_data_norm,@naneucdist);
                % Y = pdist(feat_data_norm(~ind_NaN,:),'euclidean');
                % end
                % % Cluster the data with structure defined by neighbours
                % Y = pdist(adjacency_matrix,'euclidean');
                Z = linkage(feat_data_norm_no_NaN,'ward','euclidean');
                % Y = [];
                
                for ncluster = k_cluster2test_vec
                    % disp([agglo_name num2str(ncluster)])
                    
                    clusteringSolutionAgglo = cluster(Z,'maxclust',ncluster);
                    
                    % add the cluster values to the proper non NaN positions:
                    clusteringSolutionAggloFinal = NaN(size(ind_NaN));
                    clusteringSolutionAggloFinal(~ind_NaN) = clusteringSolutionAgglo;
                    
                    if plotClustering==1
                        [centroids, n_cluster_values, dimension] = ...
                            getClusterCentroids(feat_data_norm, clusteringSolutionAgglo);
                        distances = [];
                        plotClusterSolution(feat_data_norm, clusteringSolutionAggloFinal, ...
                            centroids, distances, n_cluster_values, ncluster, ...
                            dimension, 1, feat_names)
                    else
                        centroids = getClusterCentroids(feat_data_norm, clusteringSolutionAggloFinal);
                    end
                    
                    % if ~isempty(ind_NaN)
                    %     % when the vector of data with NaN
                    %     % is given to the cluster method
                    %     clusteringSolutionAgglo(ind_NaN) = NaN;
                    % end
                    
                    % save clustering solutions and centroids *****************
                    temp_solutions{count_cluster_methods}(ff,:) = clusteringSolutionAggloFinal;
                    temp_prototypes{count_cluster_methods}(ff) = {centroids};
                    count_cluster_methods = count_cluster_methods+1;
                end
                Z = [];
                clusteringSolutionAgglo = [];
                clusteringSolutionAggloFinal = [];
                centroids = [];
            end
            agglo_name = [];
            
            
            % [H,T] = dendrogram(Z,'colorthreshold',1.2,'orientation','right');
            
            %% Run DBSCAN Clustering Algorithm on the features intervals
            dbscan_name = 'dbscan_d';
            if any(~cellfun(@isempty,strfind(clustering_methods,dbscan_name)))
                [centroids, clusteringSolutionDBSCANFinal] = getDBSCANClustering(dimension, ...
                    feat_data_norm, feat_names, time_min, plotClustering, n_epsilon4dbscan);
                
                for tt = count_cluster_methods:count_cluster_methods+n_epsilon4dbscan-1
                    temp_prototypes{tt}(ff) = centroids(tt-count_cluster_methods+1);
                    temp_solutions{tt}(ff,:) = clusteringSolutionDBSCANFinal{tt-count_cluster_methods+1};
                end
                count_cluster_methods = count_cluster_methods+n_epsilon4dbscan;
            end
            
            
            
            %% Run gaussian mixture model clustering
            
            gmm_name = 'gmm_k';
            if any(~cellfun(@isempty,strfind(clustering_methods,gmm_name)))
                
                
                plotFigure = 0;
                if plotFigure
                    close all
                    figure()
                    count = 0;
                    SCtext = {'true','false'};
                end
                
                Sigma = {'diagonal','full'};
                nSigma = numel(Sigma);
                SharedCovariance = {true,false};
                
                nSC = numel(SharedCovariance);
                options = statset('MaxIter',4000); % Increase number of EM iterations
                
                distances_matrix = pdist2(feat_data_norm_no_NaN,feat_data_norm_no_NaN);
                % nn_matrix = nearest_neighbours(distances_matrix);
                save_cluster_solutions = zeros(nSigma,nSC,length(feat_data_norm_no_NaN));
                
                DI_mat = zeros(2,2);
                OD_mat = zeros(2,2);
                % ICV_mat = zeros(2,2);
                % C_mat = zeros(2,2);
                % CS_mat = zeros(2,2);
                
                for ii = 1:nSigma
                    for jj = 1:nSC
                        
                        % before performing clustering to force it to
                        % always use the same 'random' values and return
                        % the same result each time. It doesn't make it
                        % deterministic, but it makes it repeatable at
                        % least.
                        rng 'default'
                        
                        % we want solutions with two clusters: input of 2
                        gmfit = fitgmdist(feat_data_norm_no_NaN, 2, ...
                            'CovType', Sigma{ii}, 'SharedCov', ...
                            SharedCovariance{jj}, 'Options', options, ...
                            'Regularize',1e-5);
                        % Use 'Regularize' to add a very small positive
                        % number to the diagonal of every covariance matrix.
                        
                        clusterSolution = cluster(gmfit,feat_data_norm_no_NaN);
                        save_cluster_solutions(ii,jj,:) = clusterSolution;
                        % cluster_mean = unique(round(gmfit.mu*1e4)/1e4,'rows');
                        
                        % nsamples_Larger_cluster_mat(ii,jj) = {max(histc(clusterSolution,nClusterSolution))/numel(clusterSolution)};
                        
                        [centroids, n_cluster_values, dimension] = ...
                            getClusterCentroids(feat_data_norm_no_NaN, clusterSolution);
                        n_clusters = numel(n_cluster_values);
                        % Validity indices based on cluster labels
                        
                        DI_mat(ii,jj) = dunns(distances_matrix,clusterSolution, n_clusters);
                        
                        [overallDeviation,~] = compactness(feat_data_norm_no_NaN,clusterSolution);
                        OD_mat(ii,jj) = overallDeviation;
                        % ICV_mat(ii,jj) = intraclusterVariance;
                        % L = 10;
                        % connectedness = connectivity(nn_matrix,clusterSolution,L,ncluster);
                        % C_mat(ii,jj) = connectedness;
                        % CS_mat(ii,jj) = clusterSeparation(feat_data_norm_no_NaN,clusterSolution,{cluster_mean});
                        
                        if plotFigure
                            count = count+1;
                            subplot(2,2,count)
                            plotClusterSolution(feat_data_norm_no_NaN, clusterSolution, ...
                                centroids, [], n_cluster_values, n_clusters, ...
                                dimension, 1, feat_names, [], [], [])
                            title(sprintf('SigmaType = %s, SharedCov = %s',Sigma{ii},SCtext{jj}))
                        end
                    end
                end
                nSigma = [];
                nSC = [];
                gmfit = [];
                Sigma = [];
                SharedCovariance = [];
                options = [];
                distances_matrix = [];
                n_cluster_values = [];
                
                [ind_row_best_solution, ind_col_best_solution] = find(DI_mat>0.15);
                DI_mat = [];
                if ~isempty(ind_row_best_solution)
                    if numel(ind_row_best_solution)>1
                        OD_mat_best = zeros(numel(ind_row_best_solution),1);
                        for tt = 1:numel(ind_row_best_solution)
                            OD_mat_best(tt) = OD_mat(ind_row_best_solution(tt), ind_col_best_solution(tt));
                        end
                        
                        ind_best_final = find(OD_mat_best==max(OD_mat_best));
                        if numel(ind_best_final)==1
                            ind_row_best_solution = ind_row_best_solution(ind_best_final);
                            ind_col_best_solution = ind_col_best_solution(ind_best_final);
                        else
                            ind_row_best_solution = ind_row_best_solution(ind_best_final(1));
                            ind_col_best_solution = ind_col_best_solution(ind_best_final(1));
                        end
                    end
                    ind_best_final = [];
                else
                    ind_row_best_solution = 1;
                    ind_col_best_solution = 2;
                end
                
                clusteringSolutionGMM = squeeze(save_cluster_solutions(ind_row_best_solution, ...
                    ind_col_best_solution, :));
                
                ind_row_best_solution = [];
                ind_col_best_solution = [];
                save_cluster_solutions = [];
                
                [centroids, ~, ~] = getClusterCentroids(feat_data_norm_no_NaN, ...
                    clusteringSolutionGMM);
                
                % add the cluster values to the proper non NaN positions:
                clusteringSolutionGMMFinal = NaN(size(ind_NaN));
                clusteringSolutionGMMFinal(~ind_NaN) = clusteringSolutionGMM;
                
                temp_solutions{count_cluster_methods}(ff,:) = clusteringSolutionGMMFinal;
                temp_prototypes{count_cluster_methods}(ff) = {centroids};
                count_cluster_methods = count_cluster_methods+1;
                
                clusteringSolutionGMM = [];
                clusteringSolutionGMMFinal = [];
                centroids = [];
            end
            gmm_name = [];
            
            
            
            
            
        end
        
        
        %% Multiobjective algorithm applied to all seizures at once
        
        % if ~any(find(isnan(data2cluster))) % check if there are NaN values
        moea_name = 'moea';
        if any(strcmp(clustering_methods,moea_name))
            cd('C:\Users\Adriana\OneDrive - Universidade de Coimbra\analyse_features_SCE\All SCE\SCEpreictalStudy\SCE CODE')
            plotFigures = 0;
            feats2analyse = feat_names2analyse(feat_combs)';
            clusteringSolutionMOEA = MOEA_clustering(feats2analyse(:,2000), ...
                data2cluster(61,:,:), plotFigures);
            cd(current_folder)
            
            clusteringPrototypeMOEAFinal = cell(n_comb,1);
            clusteringSolutionMOEAFinal = clusteringPrototypeMOEAFinal;
            count_ff = 1;
            for ff = 1:n_comb
                feat_data = squeeze(data2cluster(ff,:,:));
                ind_NaN = any(isnan(feat_data),2);
                clusteringSolutionMOEAFinal_temp = NaN(size(ind_NaN));
                clusteringSolutionMOEAFinal_temp(~ind_NaN) = clusteringSolutionMOEA{ff};
                clusteringSolutionMOEAFinal(ff) = {clusteringSolutionMOEAFinal_temp};
                
                feat_data_norm = (feat_data-min(feat_data))./(max(feat_data)-min(feat_data));
                % feat_data_norm = bsxfun(@rdivide, bsxfun(@minus, feat_data, min(feat_data)), ...
                %   (max(feat_data)-min(feat_data)));
                
                feat_data_norm_no_NaN = feat_data_norm(~ind_NaN,:);
                
                [centroids, n_cluster_values, dimension] = ...
                    getClusterCentroids(feat_data_norm_no_NaN, clusteringSolutionMOEA{count_ff});
                
                clusteringPrototypeMOEAFinal(ff) = {centroids};
            end
            
            temp_solutions{count_cluster_methods} = cell2mat(clusteringSolutionMOEAFinal')';
            temp_prototypes{count_cluster_methods}(ff) = clusteringPrototypeMOEAFinal;
            
            
            % clustering_solutions.(moea_name) = cell2mat(clusteringSolutionMOEA')';
            
            
            filepath = fullfile(cd, results_folder, 'computed_moea.txt');
            if exist(filepath, 'file')==2
                fileID = fopen(filepath, 'a+');
                comp_date = datestr(now);
                fprintf(fileID,'%9s %4d %s %11s\n', 'Feat Comb', ff, 'computed at', comp_date);
                fclose(fileID);
            else
                fileID = fopen(filepath,'w');
                comp_date = datestr(now);
                fprintf(fileID,'%9s %4d %s %11s\n', 'Feat Comb', ff, 'computed at', comp_date);
                fclose(fileID);
            end
        end
        % end
        
        clc
        
        %% save results
        save_clustering_data.data2cluster = data2cluster;
        
        % save_clustering_data.clustering_solutions = temp_solutions;
        % save_clustering_data.clustering_prototypes = temp_prototypes;
        
        for tt = 1:numel(methods_vec)
            save_clustering_data.clustering_solutions.(clustering_methods{methods_vec(tt)}) = temp_solutions{tt};
            save_clustering_data.clustering_prototypes.(clustering_methods{methods_vec(tt)}) = temp_prototypes{tt};
        end
        
        parsave(seizFolderPath, {save_clustering_data}, 'save_clustering_data')

    end
    
    
end


end