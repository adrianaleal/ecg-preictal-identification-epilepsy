function [indexes_sorted_distances, sorted_distances] = nearest_neighbours(distances_matrix)


% sort each distance column vector from the lowest to highest distance
% in order words, the first point is separated from the other points by the
% value of each row in the first column in rightB, and so on
% in rightI it is possible to find the corresponding index of the point
% which distance is found in rightB



[sorted_distances,indexes_sorted_distances] = sort(distances_matrix);
% if nargin<2 
    n_points = length(sorted_distances);
% end



indexes_sorted_distances_final = indexes_sorted_distances;
sorted_distances_final = sorted_distances;
if ~isequal(indexes_sorted_distances(1,:),(1:n_points))
    % there are two samples with the same value (two samples with zero distance)
    indexes_cc = find(indexes_sorted_distances(1,:)~=(1:n_points));
    for dd = 1:length(indexes_cc)
        indexes_rr = find(indexes_sorted_distances(:,indexes_cc(dd))==indexes_cc(dd));
        indexes_sorted_distances_final(indexes_rr,indexes_cc(dd)) = indexes_sorted_distances(1,indexes_cc(dd));
        indexes_sorted_distances_final(1,indexes_cc(dd)) = indexes_sorted_distances(indexes_rr,indexes_cc(dd));
        sorted_distances_final(indexes_rr,indexes_cc(dd)) = sorted_distances(1,indexes_cc(dd));
        sorted_distances_final(1,indexes_cc(dd)) = sorted_distances(indexes_rr,indexes_cc(dd));
    end
end

% as the first row corresponds to the distance of the point to itself
% (which is zero), the first row of indexes_sorted_distances must 
% correspond to the index of that point in distance matrix column
if ~isequal(indexes_sorted_distances_final(1,:),(1:n_points))
    disp('#######################################################')
    disp('the matrix is NOT right')
    disp('#######################################################')
end

end