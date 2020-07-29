%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function [IDX, isnoise] = DBSCAN(X,epsilon,MinPts, plotFigure)

    if isempty(epsilon)
        epsilon = getEpsilon(X, plotFigure)
    end
    
    C=0;
    
    n=size(X,1);
    IDX=zeros(n,1);
    
    D=pdist2(X,X);
    
    visited=false(n,1);
    isnoise=false(n,1);
    
    for ii=1:n
        if ~visited(ii)
            visited(ii)=true;
            
            Neighbors=RegionQuery(ii);
            if numel(Neighbors)<MinPts
                % X(i,:) is NOISE
                isnoise(ii)=true;
            else
                C=C+1;
                ExpandCluster(ii,Neighbors,C);
            end
            
        end
    end
    
    function ExpandCluster(ii,Neighbors,C)
        IDX(ii)=C;
        
        k = 1;
        while true
            jj = Neighbors(k);
            
            if ~visited(jj)
                visited(jj)=true;
                Neighbors2=RegionQuery(jj);
                if numel(Neighbors2)>=MinPts
                    Neighbors=[Neighbors Neighbors2];   %#ok
                end
            end
            if IDX(jj)==0
                IDX(jj)=C;
            end
            
            k = k + 1;
            if k > numel(Neighbors)
                break;
            end
        end
    end
    
    function Neighbors=RegionQuery(ii)
        Neighbors=find(D(ii,:)<=epsilon);
    end

end


function epsilon_dist = getEpsilon(X, plotFigure)

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

[nn_indexes, nn_sorted_distances] = nearest_neighbours(distances_matrix);

mean_nn_sorted_distances = mean(nn_sorted_distances,2);

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
% nn_sorted_distances_diff = diff(nn_sorted_distances);
% figure()
% plot(repmat((1:n_points-1)',1,n_points),nn_sorted_distances_diff)
% xlabel('Nearest neighbours')
% ylabel('Derivative of the nearest neighbours matrix')
% axis tight
% [maximum_for_each_point,I] = max(nn_sorted_distances_diff);
% index_point_max_step_der = find(maximum_for_each_point==max(maximum_for_each_point));
% index_distance_max_step_der = I(index_point_max_step_der)+1;

% using the auc
auc_under_each_curve = trapz(nn_sorted_distances);
index_point_max_step_auc = find(auc_under_each_curve==max(auc_under_each_curve));
y = nn_sorted_distances(:,index_point_max_step_auc);
y = (y-min(y))/(max(y)-min(y));
x = linspace(0,1,48)';
geometric_mean = sqrt((1-x).^2 + y.^2);
[~,index_distance_max_step_auc] = max(geometric_mean,[],1)

if plotFigure
    figure()
    plot(X(:,1), X(:,2), 'b*'), hold on
    plot(X(index_point_max_step_auc,1), X(index_point_max_step_auc,2), 'ro'), hold off
end

% using the auc on the mean
y = mean_nn_sorted_distances;
y = (y-min(y))/(max(y)-min(y));
geometric_mean = sqrt((1-x).^2 + y.^2);
[~,index_distance_max_step_auc_mean] = min(geometric_mean,[],1)


if plotFigure
    figure()
    plot(repmat((1:n_points)',1,n_points),nn_sorted_distances), hold on
    h1 = plot((1:n_points)',mean_nn_sorted_distances,'b--');
    % h1 = plot(repmat((1:n_points)',1,numel(max_step_curve_point_index)), ...
    %     nn_sorted_distances(:,max_step_curve_point_index),'k--');
    h2 = plot(repmat((1:n_points)',1,numel(index_point_max_step_auc)), ...
        nn_sorted_distances(:,index_point_max_step_auc),'r*');
    h3 = plot(repmat(index_distance_max_step_auc,1,numel(index_point_max_step_auc)),...
        nn_sorted_distances(index_distance_max_step_auc, index_point_max_step_auc),'bo');
    h4 = plot(index_distance_max_step_auc_mean, mean_nn_sorted_distances(index_distance_max_step_auc_mean),'ro');
    hold off
    %     ylabel(['k-NN distance from point' num2str(pp)])
    ylabel('k-NN distance from each point')
    xlabel('Nearest neighbour')
    legend([h1 h2, h3, h4], [{'mean curve found with auc'},{'curve found with auc'}, ...
        {'point found with auc'}, {'mean curve found with auc'}, {'mean point found with auc'}])
end

if numel(index_distance_max_step_auc_mean)==1
%     if ismember(index_point_max_step_auc,max_step_curve_point_index)
%         final_index_point_max_step = index_point_max_step_auc;
%         final_index_distance_max_step = max_step_curve_distance_index(max_step_curve_point_index==index_point_max_step_auc);
%     end
else
    disp('AUC has found more than one curve to define threshold --> COMPLETE CODE')
end


epsilon_maximum_auc_curve = nn_sorted_distances(index_distance_max_step_auc, index_point_max_step_auc);
epsilon_mean = mean_nn_sorted_distances(index_distance_max_step_auc_mean);

epsilon_dist = abs(mean(getDist)-std(getDist));


end










