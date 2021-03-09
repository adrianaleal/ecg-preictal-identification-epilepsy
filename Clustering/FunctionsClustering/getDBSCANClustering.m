function clusteringSolutionDBSCAN = getDBSCANClustering(dimension, ...
    data2cluster, epsilon)



% Inputs:
% - dimension (double): number of features to apply clustering
% - data2cluster (double): matrix whose columns contain information for
%                          each feature
% - epsilon: radius Euclidean distance that establishes a neighbourhood
%            among samples


% Outputs:
% - centroids (double)
% - clusteringSolutionDBSCAN (double)

% NOTES:
% the cluster solution of DBSCAN has zeros corresponding to noisy points

%% Run DBSCAN Clustering Algorithm

if dimension==2
    % For two-dimensional data: use default value of minPts=4
    % (Ester et al., 1996)
    MinPts_DBSCAN = 4;
elseif dimension>2
    % For more than 2 dimensions: minPts=2*dim
    % (Sander et al., 1998)
    MinPts_DBSCAN = 2*dimension;
end


% epsilon: Two points are considered neighbours if the distance between the
% two points is below the threshold epsilon.

if isempty(epsilon)
    getkPlot = 0;
    epsilon = getEpsilonDBSCAN(data2cluster, getkPlot);
end

clusteringSolutionDBSCAN = DBSCAN(data2cluster, epsilon, MinPts_DBSCAN);


end