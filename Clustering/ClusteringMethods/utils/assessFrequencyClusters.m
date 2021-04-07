function [clusteringSolutionChanged, frequency_clusters] = ...
    assessFrequencyClusters(clusteringSolution, frequency_clusters)

clusteringSolutionChanged = clusteringSolution;

if isequal(frequency_clusters(1,2), frequency_clusters(2,2))
    % when the two clusters have the same amount of samples
    % assign class 2 to the cluster with the beggining
    % sample that is closest to the seizure:
    
    indexes_cluster1 = find(clusteringSolution==frequency_clusters(1,1));
    indexes_cluster2 = find(clusteringSolution==frequency_clusters(2,1));
    
    if indexes_cluster1(1)>indexes_cluster2(1)
        clusteringSolutionChanged = clusteringSolution;
        clusteringSolutionChanged(clusteringSolution==frequency_clusters(1,1)) = frequency_clusters(2,1);
        clusteringSolutionChanged(clusteringSolution==frequency_clusters(2,1)) = frequency_clusters(1,1);
    end
else
    % assign class 1 with the highest amount of samples if
    % that is not already the case:
    if isequal(sort(frequency_clusters),frequency_clusters)
        % change the cluster classes
        frequency_clusters = [sort(frequency_clusters(:,1)) ...
            sort(frequency_clusters(:,2),'descend')];

        for ii = 1:size(frequency_clusters,1)
            clusteringSolutionChanged(clusteringSolution==frequency_clusters(ii,1)) = frequency_clusters(end-ii+1,1);
        end
        % see = [clusteringSolution, clusteringSolutionChanged];
    end
end

end