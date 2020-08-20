function clustering_methods = replaceClusteringMethodsNames(clustering_methods)

index = find(~cellfun(@isempty,strfind(clustering_methods,'kmeans_k2')));
if ~isempty(index)
    for ll = 1:numel(index)
        clustering_methods(index(ll)) = {'KMEANS'};
    end
end

index = find(strcmp(clustering_methods,'agglo_hier_k2'));
if ~isempty(index)
    for ll = 1:numel(index)
        clustering_methods(index(ll)) = {'AH'};
    end
end

IndexC = strfind(clustering_methods,'dbscan_d');
index = find(not(cellfun('isempty',IndexC)));

if ~isempty(index)
    for ll = 1:numel(index)
        clustering_methods(index(ll)) = {[strrep(clustering_methods{index(ll)},'dbscan_d','DBSCAN$_') '$']};
    end
end


IndexC = strfind(clustering_methods,'gmm_k2');
index = find(not(cellfun('isempty',IndexC)));

if ~isempty(index)
    for ll = 1:numel(index)
        clustering_methods(index(ll)) = {'GMM'};
    end
end

end

