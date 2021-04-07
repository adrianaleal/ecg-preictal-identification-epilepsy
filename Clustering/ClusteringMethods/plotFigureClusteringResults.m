function [frequency_smaller_cluster] = plotFigureClusteringResults( ...
    ss, clusteringSolution, n_cluster_values, n_seizures_pat, ...
    feat_names, feat2plot, title_name)

nclusters = numel(n_cluster_values);

[row, col] = size(feat2plot);
if col>row
    feat2plot = feat2plot';
end


[row1, col1] = size(n_cluster_values);
if col1>row1
    n_cluster_values = n_cluster_values';
end

frequency_clusters = [n_cluster_values histc(clusteringSolution, n_cluster_values)];

[clusteringSolution, frequency_clusters] = assessFrequencyClusters( ...
    clusteringSolution, frequency_clusters);

frequency_smaller_cluster = frequency_clusters(end,2);


[p,n] = numSubplots(n_seizures_pat);


% time from 240 min before until the seizure onset:
time_vec = linspace(0,240,size(feat2plot,1));

ind_NaN = any(isnan(feat2plot),2);
all_data = feat2plot(~ind_NaN,:);
clustering_solution_no_nan = clusteringSolution(~ind_NaN);
colors = time_vec(~ind_NaN)';
markers = '*oxds+ph.';


figure(67)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
subplot(p(1),p(2),ss)

count = 1;
for idx = n_cluster_values' % because DBSCAN has clusters with zeros (noise)
    ind_cluster_solution = clustering_solution_no_nan==idx;
    data = all_data(ind_cluster_solution,:);
    scatter3(data(:,1), data(:,2), data(:,3), [], ...
        colors(ind_cluster_solution), markers(count));
    count = count+1;
    hold on;
end
hold off


if ss==1 || ss==((p(1)*p(2))/2+1)
    ax = gca;
    pos = ax.Position;
    cb = colorbar();
    set(cb, 'Position',  [pos(1)-0.05 pos(2) 0.01 pos(4)])
    ylabel(cb,'Time before seizure (min)')
    n = 6;
    vec = 0:240/n:240+2;
    
    cb_tick_labels = cellstr(num2str((240:-40:0)'));
    set(cb,'XTickLabel',cb_tick_labels, 'XTick', vec);
    title(cb,'\bf{SEIZURE}');
end

set(gca,'FontSize',14)

title(title_name)

feat_labels = replaceFeatureNames(feat_names);

xlabel(feat_labels(1))
ylabel(feat_labels(2))
zlabel(feat_labels(3))
legend(num2str(n_cluster_values),'Location','best')

distances_matrix = pdist2(feat2plot, feat2plot);
DI = dunns(clusteringSolution, nclusters, distances_matrix)

if all(n_cluster_values>0) && nclusters==2 && DI>=0.15
    legend({[num2str(frequency_clusters(1,1)) ' - Interictal'], ...
        [num2str(frequency_clusters(2,1)) ' - Preictal']},'Location','best')
end

end