function clusteringSolutionGMM = getGMMClustering(data2cluster, ...
    ncluster, plotFigure)

% Inputs:
% - data2cluster (double): matrix whose columns contain information for
%                          each feature
% - plotFigure (double): flag to plot the clustering result


% Outputs:
% - clusteringSolutionGMM (double)

% NOTES:
% When not to use Gaussian Mixture Model (EM clustering):
% - Non-Gaussian dataset: as it is clear from the formulation, GMM assumes
% an underlying Gaussian generative distribution. However, many practical
% datasets do not satisfy this assumption.
% - Uneven cluster sizes. When clusters do not have even sizes there is a
% high chance that small cluster gets dominated by the large one.
% Source: http://hameddaily.blogspot.com/2015/03/when-not-to-use-gaussian-mixtures-model.html

%%

if plotFigure
    figure(12)
    count = 0;
    SCtext = {'true', 'false'};
end

Sigma = {'diagonal', 'full'};
nSigma = numel(Sigma);
SharedCovariance = {true, false};
nSC = numel(SharedCovariance);
options = statset('MaxIter', 4000); % Increase number of EM iterations

distances_matrix = pdist2(data2cluster, data2cluster);
% nn_matrix = nearest_neighbours(distances_matrix);

save_cluster_solutions = zeros(nSigma, nSC, length(data2cluster));

OD_mat = zeros(nSigma,nSC);
DI_mat = OD_mat;

for ii = 1:nSigma
    for jj = 1:nSC
        
        % before performing clustering to force it to
        % always use the same 'random' values and return
        % the same result each time. It doesn't make it
        % deterministic, but it makes it repeatable at
        % least.
        rng 'default'
        
        % we want solutions with two clusters: input of 2
        gmfit = fitgmdist(data2cluster, ncluster, 'CovType', Sigma{ii}, ...
            'SharedCov', SharedCovariance{jj}, 'Options', options, ...
            'RegularizationValue', 1e-5);
        % Use 'Regularize' to add a very small positive
        % number to the diagonal of every covariance matrix.
        
        clusterSolution = cluster(gmfit, data2cluster);
        save_cluster_solutions(ii,jj,:) = clusterSolution;
        % cluster_mean = unique(round(gmfit.mu*1e4)/1e4,'rows');
        
        % nsamples_Larger_cluster_mat(ii,jj) = {max(histc(clusterSolution,nClusterSolution))/numel(clusterSolution)};
        
        [centroids, n_cluster_values, dimension] = ...
            getClusterCentroids(data2cluster, clusterSolution);
        n_clusters = numel(n_cluster_values);
        
        % Validity indices based on cluster labels
        DI_mat(ii,jj) = dunns(clusterSolution, n_clusters, distances_matrix);
        
        % Validity indices based on cluster prototypes
        [overallDeviation,~] = compactness( ...
            clusterSolution, n_clusters, data2cluster);
        OD_mat(ii,jj) = overallDeviation;
        
        if plotFigure
            figure(12)
            count = count+1;
            subplot(2,2,count)
            distances = [];
            plotClusterSolution(data2cluster, clusterSolution, ...
                centroids, distances, dimension, 1, [], [], [], []);
            title(sprintf('SigmaType = %s, SharedCov = %s', Sigma{ii}, ...
                SCtext{jj}))
        end
    end
end


% Select the best solution among the four possible solutions
[ind_row_best_solution, ind_col_best_solution] = find(DI_mat>0.15);

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


end