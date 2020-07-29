function [clusteringSolutionGMM, centroids] = ...
    applyGaussianMixtureModelsClustering(data2cluster, plotFigure, ...
    feat_names)


if plotFigure
    figure()
    count = 0;
    SCtext = {'true','false'};
end

Sigma = {'diagonal','full'};
nSigma = numel(Sigma);
SharedCovariance = {true,false};

nSC = numel(SharedCovariance);
options = statset('MaxIter',10000); % Increase number of EM iterations

distances_matrix = pdist2(data2cluster,data2cluster);
% nn_matrix = nearest_neighbours(distances_matrix);
save_cluster_solutions = zeros(nSigma,nSC,length(data2cluster));

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
        gmfit = fitgmdist(data2cluster, 2, ...
            'CovType', Sigma{ii}, 'SharedCov', ...
            SharedCovariance{jj}, 'Options', options, ...
            'Regularize',1e-5);
        % Use 'Regularize' to add a very small positive
        % number to the diagonal of every covariance matrix.
        
        clusterSolution = cluster(gmfit,data2cluster);
        save_cluster_solutions(ii,jj,:) = clusterSolution;
        % cluster_mean = unique(round(gmfit.mu*1e4)/1e4,'rows');
        
        % nsamples_Larger_cluster_mat(ii,jj) = {max(histc(clusterSolution,nClusterSolution))/numel(clusterSolution)};
        
        [centroids, n_cluster_values, dimension] = ...
            getClusterCentroids(data2cluster, clusterSolution);
        n_clusters = numel(n_cluster_values);
        % Validity indices based on cluster labels
        
        DI_mat(ii,jj) = dunns(distances_matrix,clusterSolution, n_clusters);
        
        [overallDeviation,~] = compactness(data2cluster,clusterSolution);
        OD_mat(ii,jj) = overallDeviation;
        % ICV_mat(ii,jj) = intraclusterVariance;
        % L = 10;
        % connectedness = connectivity(nn_matrix,clusterSolution,L,ncluster);
        % C_mat(ii,jj) = connectedness;
        % CS_mat(ii,jj) = clusterSeparation(feat_data_norm_no_NaN,clusterSolution,{cluster_mean});
        
        if plotFigure
            count = count+1;
            subplot(2,2,count)
            plotClusterSolution(data2cluster, clusterSolution, ...
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

[centroids, ~, ~] = getClusterCentroids(data2cluster, ...
    clusteringSolutionGMM);


end