function h1 = plotClusterSolution(data2cluster, clusterSolution, ...
    centroids, distances, n_cluster_values, ncluster, dimension, ...
    plotOnlyClustering, feat_names, time, colors, markers)



% colors = ['r', 'b', 'g','m','k','c'];
h = zeros(1,ncluster);
h2 = zeros(1,ncluster);
marker = ['o','d','s','+','v','x','*','o','d','s','+','v','x','*','o','d','s','+','v','x','*'];
marker_size = [14 30 40 26 30 20 22 24 14 30 40 26 30 20 22 24 14 30 40 26 30 20 22 24];

if plotOnlyClustering
    if dimension==3
        for ll = 1:ncluster
            
            if isempty(markers)
                if isempty(colors)
                    colors = jet(ncluster);
                end
                h1(ll) = scatter3(data2cluster(clusterSolution==n_cluster_values(ll),1), ...
                    data2cluster(clusterSolution==n_cluster_values(ll),2),...
                    data2cluster(clusterSolution==n_cluster_values(ll),3),[],colors(ll,:));
                %         h(ll) = scatter3(data2cluster(clusterSolution==ll,1),data2cluster(clusterSolution==ll,2),...
                %             data2cluster(clusterSolution==ll,3),marker_size(ll),colors(ll,:),marker(ll));
                hold on
                if ~isempty(centroids)
                    h2(ll) = plot3(centroids(ll,1),centroids(ll,2),centroids(ll,3),'kx','MarkerSize',15,'LineWidth',3);
                end
            else
                h1(ll) = scatter3(data2cluster(clusterSolution==n_cluster_values(ll),1), ...
                    data2cluster(clusterSolution==n_cluster_values(ll),2),...
                    data2cluster(clusterSolution==n_cluster_values(ll),3),[],colors(ll,:), ...
                    markers(ll));
                
                %         h(ll) = scatter3(data2cluster(clusterSolution==ll,1),data2cluster(clusterSolution==ll,2),...
                %             data2cluster(clusterSolution==ll,3),marker_size(ll),colors(ll,:),marker(ll));
                hold on
                if ~isempty(centroids)
                    h2(ll) = plot3(centroids(ll,1),centroids(ll,2),centroids(ll,3), ...
                        'k','Marker',markers(ll), 'MarkerSize',15,'LineWidth',3);
                end
            
            end
            

        end
        xlabel(regexprep(feat_names{1},'_',' '))
        ylabel(regexprep(feat_names{2},'_',' '))
        zlabel(regexprep(feat_names{3},'_',' '))
        hold off
    else
        for ll = 1:ncluster
            colors = jet(ncluster);
            h1(ll) = scatter(data2cluster(clusterSolution==n_cluster_values(ll),1), ...
                data2cluster(clusterSolution==n_cluster_values(ll),2),[],colors(ll,:));
            %         h(ll) = scatter3(data2cluster(clusterSolution==ll,1),data2cluster(clusterSolution==ll,2),...
            %             data2cluster(clusterSolution==ll,3),marker_size(ll),colors(ll,:),marker(ll));
            hold on
            if ~isempty(centroids)
                h2(ll) = plot(centroids(ll,1),centroids(ll,2),'kx','MarkerSize',15,'LineWidth',3);
            end
            xlabel(regexprep(feat_names{1},'_',' '))
            ylabel(regexprep(feat_names{2},'_',' '))
        end
        hold off
        
    end
    
else
    
    if dimension==3 %&& ncluster<19
        figure()
        set(gcf,'units','normalized','outerposition',[0 0 1 0.6])
%         time = 5:5:240;% linspace(240,0,size(data2cluster,1));
        ax1 = subplot(131);
        colormap(ax1,'default')
        for ll = 1:ncluster
            h(ll) = scatter3(data2cluster(:,1),data2cluster(:,2),data2cluster(:,3),[],time,'filled');
            %         h(ll) = scatter3(data2cluster(clusterSolution==ll,1),data2cluster(clusterSolution==ll,2),...
            %             data2cluster(clusterSolution==ll,3),marker_size(ll),c(clusterSolution==ll),marker(ll));
            %         hold on
            %         if ~isempty(centroids)
            %             h2(ll) = plot3(centroids(ll,1),centroids(ll,2),centroids(ll,3),'kx','MarkerSize',15,'LineWidth',3);
            %         end
            %         holf off
        end
        ax1Pos = get(ax1,'Position');
        cb1 = colorbar('westoutside');
        set(ax1,'Position',ax1Pos);
        ylabel(cb1, 'Time before seizure (min)','Interpreter','Latex','FontSize',14)
        cb1Pos = get(cb1,'Position');
        cb1Pos(1) = cb1Pos(1)-0.03;
        set(cb1,'Position',cb1Pos);
        cb1.Label.Interpreter = 'Latex';
        index = find(~cellfun(@isempty,strfind(feat_names,'alpha')));
        if ~isempty(index)
            for oo = 1:numel(index)
                feat_names(index(oo)) = {strrep(feat_names{index(oo)}, ...
                    'alpha','$\alpha$')};
            end
        end
        
        xlabel(regexprep(feat_names{1},'_',' '))
        ylabel(regexprep(feat_names{2},'_',' '))
        zlabel(regexprep(feat_names{3},'_',' '))
        title('\textbf{(a)}')
        
        
        ax2 = subplot(132);
        %      colormap(lines(ncluster))
        
        map = brewermap(ncluster,'Spectral');
        for ll = 1:ncluster
            colors = jet(ncluster);
            h1(ll) = scatter3(data2cluster(clusterSolution==n_cluster_values(ll),1), ...
                data2cluster(clusterSolution==n_cluster_values(ll),2),...
                data2cluster(clusterSolution==n_cluster_values(ll),3),[],map(ll,:),'filled');
            %         h(ll) = scatter3(data2cluster(clusterSolution==ll,1),data2cluster(clusterSolution==ll,2),...
            %             data2cluster(clusterSolution==ll,3),marker_size(ll),colors(ll,:),marker(ll));
            hold on
            if ~isempty(centroids)
                h2(ll) = plot3(centroids(ll,1),centroids(ll,2),centroids(ll,3),'kx','MarkerSize',15,'LineWidth',3);
            end
        end
        hold off
        %     colormap(ax2,jet(ncluster))
        
        %         legend([h(1), h(2), h2(2)],{'Cluster 1','Cluster 2','Centroids'},'Location','NW')
        %     title (['Cluster Assignments and Centroids with ncluster = ' num2str(ncluster)])
        
        xlabel(regexprep(feat_names{1},'_',' '))
        ylabel(regexprep(feat_names{2},'_',' '))
        zlabel(regexprep(feat_names{3},'_',' '))
        title('\textbf{(b)}')
        
        
        ax3 = subplot(133);
        imagesc(distances,'CDataMapping','scaled')
        set(gca,'XDir','reverse')
        time_axis = linspace(240,0,7);
        set(gca,'XTickLabel',time_axis(2:end-1))
        set(gca,'YTickLabel',time_axis(2:end-1))
        colormap(ax3,map)
        if ncluster>1
%             caxis([min(n_cluster_values) max(n_cluster_values)]);
            cbh = colorbar('YTick',1+0.5*(ncluster-1)/ncluster:(ncluster-1)/ncluster:max(n_cluster_values)+2,...
            'YTickLabel',int2str(n_cluster_values'), 'YLim', [1 ...
            max(n_cluster_values)+1],'Location','EastOutside');
        else
            cbh = colorbar();
        end
        
        cbh.Label.Interpreter = 'latex';
        ylabel(cbh, 'Cluster','interpreter','latex','FontSize',14)
        title('\textbf{(c)}')
    elseif dimension==2 && ncluster<7
        %     figure()
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        time = linspace(240,0,size(data2cluster,1));
        %%
        ax1 = subplot(131);
        for ll = 1:ncluster
            h(ll) = scatter(data2cluster(clusterSolution==n_cluster_values(ll),1), ...
                data2cluster(clusterSolution==n_cluster_values(ll),2),...
                marker_size(ll),time(clusterSolution==n_cluster_values(ll)),marker(ll)); hold on
            %         if ~isempty(centroids)
            %             h2(ll) = plot(centroids(ll,1),centroids(ll,2),'kx','MarkerSize',15,'LineWidth',3);
            %         end
        end
        hold off
        colormap(ax1);
        cbh = colorbar;
        ylabel(cbh, 'Time before seizure (min)','interpreter','latex','FontSize',14)
        %%
        colors = lines(ncluster);
        subplot(132)
        for ll = 1:ncluster
            h(ll) = scatter(data2cluster(clusterSolution==n_cluster_values(ll),1), ...
                data2cluster(clusterSolution==n_cluster_values(ll),2),...
                marker_size(ll),colors(ll,:),marker(ll)); hold on
            if ~isempty(centroids)
                h2(ll) = plot(centroids(ll,1),centroids(ll,2),'kx','MarkerSize',15,'LineWidth',3);
            end
        end
        %         legend([h(1), h(2), h2(2)],{'Cluster 1','Cluster 2','Centroids'},'Location','NW')
        title(['Cluster Assignments and Centroids with ncluster = ' num2str(ncluster)])
        hold off
        
        ax3 = subplot(133);
        
        distances = pdist2(clusterSolution,clusterSolution)+1;
        imagesc(distances,'CDataMapping','scaled');
        colormap(ax3,colors)
        caxis([1 ncluster]);
        time_axis = linspace(240,0,7);
        set(gca,'XTickLabel',time_axis(2:end-1))
        set(gca,'YTickLabel',time_axis(2:end-1))
        cbh = colorbar('YTick',[1+0.5*(ncluster-1)/ncluster:(ncluster-1)/ncluster:ncluster],...
            'YTickLabel',int2str([1:ncluster]'), 'YLim', [1 ncluster]);
        ylabel(cbh, 'Cluster','interpreter','latex','FontSize',14)
        
    end
end


end