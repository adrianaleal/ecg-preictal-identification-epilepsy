function epsilon_dist = getEpsilonDBSCAN(X, plotFigure)

% Find epsilon parameter of DBSCAN:
% compute a k-distance plot of your dataset (compute the k-nearest
% neighbors (k-NN) for each data point to understand what is the
% density distribution of your data, for different k. the KNN is handy
% because it is a non-parametric method. Once you choose a minPTS
% (which strongly depends on your data), you fix k to that value. Then
% you use as epsilon the k-distance corresponding to the area of the
% k-distance plot (for your fixed k) with a low slope.
% Source: https://stackoverflow.com/questions/12893492/choosing-eps-and-minpts-for-dbscan-r

distances_matrix = pdist2(X,X);
upper_values = triu(ones(size(X,1)),1);
[ind_row_upper_values, ind_col_upper_values] = find(upper_values);
getDist = distances_matrix(logical(upper_values));

% get the sorted distances between each point and its neighbours
[nn_indexes, nn_sorted_distances] = nearest_neighbours(distances_matrix);

% for each neighbour (each row), sort the distances
sorted_dist_for_nth_neighbour = sort(nn_sorted_distances,2)';

mean_nn_sorted_distances = mean(sorted_dist_for_nth_neighbour);

n_points = size(X,1);

% nnn = 7; % number of nearest neighbours
% figure()
% plot(1:n_points,nn_sorted_distances(nnn,:),'*')
% ylabel([num2str(nnn) '-NN distance'])
% xlabel('Points (samples) sorted by distance')

% Having the distance between the points, the ideia is to find which
% point presents the nearest neighbour curve with the highest step
% since it means that there is a maximum separation between neighbours
% for that patient (I did that both using the derivative and the auc
% which is more straighforward even thoug it does not give the value of
% the step to assign no epsilon)
% Source: https://www.datanovia.com/en/lessons/dbscan-density-based-clustering-essentials/#method-for-determining-the-optimal-eps-value

% % using the derivative:
% sorted_dist_for_nth_neighbour_diff = diff(sorted_dist_for_nth_neighbour);
% figure()
% plot(1:n_points-1, sorted_dist_for_nth_neighbour_diff(:,2))
% xlabel('Nearest neighbours')
% ylabel('Derivative of the nearest neighbours matrix')
% axis tight
% [maximum_for_each_point,I] = max(sorted_dist_for_nth_neighbour_diff);
% index_point_max_step_der = find(maximum_for_each_point==max(maximum_for_each_point));
% index_distance_max_step_der = I(index_point_max_step_der)+1;

% using the auc

y = 1-sorted_dist_for_nth_neighbour(:,2);
x = 1-(1:n_points)';
geometric_mean = sqrt(x.^2 + y.^2);
[~,index_distance_max_step_auc] = max(geometric_mean,[],1)

% using the auc on the mean
y = 1-mean_nn_sorted_distances';
y = (y-min(y))/(max(y)-min(y));
x = 1-(1:n_points)';
x = (x-min(x))/(max(x)-min(x));
geometric_mean = sqrt(x.^2 + y.^2);
[~,index_distance_max_step_auc_mean] = max(geometric_mean,[],1)

figure()
plot(x, y), hold on

if plotFigure
    figure()
    plot(1:n_points,sorted_dist_for_nth_neighbour), hold on
    h1 = plot(1:n_points,mean_nn_sorted_distances,'b--');
    %h2 = plot(index_distance_max_step_auc_mean, mean_nn_sorted_distances(index_distance_max_step_auc_mean),'ro');
    hold off
    %     ylabel(['k-NN distance from point' num2str(pp)])
    ylabel('k-NN distance from each point')
    xlabel('Nearest neighbour')
    %legend([h1 h2], [{'mean curve found with auc'},{'mean point found with auc'}])
end


epsilon_dist = abs(mean(getDist)+std(getDist));


end


