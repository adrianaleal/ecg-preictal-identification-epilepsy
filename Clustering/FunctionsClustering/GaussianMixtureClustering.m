function [validity_indices_out, chosenClusteringSolution] = GaussianMixtureClustering( ...
    data2cluster, plotFigures, ncluster, feat_names)

% When not to use Gaussian Mixture Model (EM clustering):
% - Non-Gaussian dataset: as it is clear from the formulation, GMM assumes 
% an underlying Gaussian generative distribution. However, many practical 
% datasets do not satisfy this assumption. 
% - Uneven cluster sizes. When clusters do not have even sizes there is a 
% high chance that small cluster gets dominated by the large one. 
% Source: http://hameddaily.blogspot.com/2015/03/when-not-to-use-gaussian-mixtures-model.html


SCtext = {'true', 'false'};


% ncluster = 4; % number of clusters
Sigma = {'diagonal', 'full'};
nSigma = numel(Sigma);
SharedCovariance = {true, false};
nSC = numel(SharedCovariance);
options = statset('MaxIter', 4000); % Increase number of EM iterations

if plotFigures
    d = 50;
    x1 = linspace(min(data2cluster(:,1)) - 2,max(data2cluster(:,1)) + 2,d);
    x2 = linspace(min(data2cluster(:,2)) - 2,max(data2cluster(:,2)) + 2,d);
    [x1grid, x2grid] = meshgrid(x1, x2);
    
    X0 = [x1grid(:) x2grid(:)];
    threshold = sqrt(chi2inv(0.99,2));
    
    
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    c = 1;
end

OD_mat = zeros(nSigma,nSC);
ICV_mat = OD_mat;
C_mat = OD_mat;
DI_mat = OD_mat;
CS_mat = OD_mat;
nclusters_mat = OD_mat;

options_names = cell(2,2);
nsamples_Larger_cluster_mat = options_names;

save_h1 = cell(2,2);
distances_matrix = pdist2(data2cluster,data2cluster);
nn_matrix = nearest_neighbours(distances_matrix);

save_cluster_solutions = zeros(nSigma, nSC, length(data2cluster));
save_gmfit = cell(nSigma, nSC);

for ii = 1:nSigma
    for jj = 1:nSC
        
        
        % contrarily to the clustering solutions in the preictal study, I
        % don't know how many clusters would be correct in the the RR
        % interval identification
        % see getGMMClustering
        options_names(ii,jj) = {['CovType ',Sigma{ii}, ' SharedCov ',SCtext{jj}]};
        stop = 0;
        try
            while stop==0
                gmfit = fitgmdist(data2cluster, ncluster, 'CovType', ...
                    Sigma{ii}, 'SharedCov', SharedCovariance{jj}, ...
                    'Options', options, 'RegularizationValue',0.01);
                
                clusterSolution = cluster(gmfit, data2cluster);
                nClusterSolution = unique(clusterSolution);
                if size(gmfit.mu,1)>numel(nClusterSolution)
                    ncluster = numel(nClusterSolution);
                else
                    stop = 1;
                end
            end
        catch
            continue
        end
        
        save_gmfit(ii,jj) = {gmfit};
        save_cluster_solutions(ii,jj,:) = clusterSolution;
        cluster_mean = unique(round(gmfit.mu*1e4)/1e4,'rows')
        
        nclusters_mat(ii,jj,:) = size(cluster_mean,1);
        
        %% Validity indices based on cluster labels
        [overallDeviation,intraclusterVariance] = compactness(data2cluster,clusterSolution);
        OD_mat(ii,jj) = overallDeviation;
        ICV_mat(ii,jj) = intraclusterVariance;
        
        L = 10;
        connectedness = connectivity(nn_matrix,clusterSolution,L,ncluster);
        C_mat(ii,jj) = connectedness;
        
        DI_mat(ii,jj) = dunns(clusterSolution, distances_matrix, ncluster);
 
        CS_mat(ii,jj) = clusterSeparation(data2cluster,clusterSolution,{cluster_mean});

        
        %%
        
        nsamples_Larger_cluster_mat(ii,jj) = {max(histc(clusterSolution,nClusterSolution))/numel(clusterSolution)};
        if plotFigures
            subplot(2,2,c);
            h1 = gscatter(data2cluster(:,1),data2cluster(:,2),clusterSolution);
            
            save_h1(ii,jj) = {h1};
            % D = mahal(obj,X) computes the Mahalanobis distance (in
            % squared units) of each observation in X to the mean of each
            % of the k components of the Gaussian mixture distribution
            % defined by obj
            
            mahalDist = mahal(gmfit,X0);
            hold on
            Legends = {};
            
            
            % if numel(nClusterSolution)<ncluster
            %    disp('problem')
            % end
            
            for m = 1:numel(nClusterSolution)
                idx = mahalDist(:,nClusterSolution(m))<=threshold;
                Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
                h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
                uistack(h2,'bottom');
                Legends{end+1} = num2str(nClusterSolution(m));
            end
            
            
            if size(cluster_mean,1)>numel(nClusterSolution)
                disp('problem')
            end
            plot(cluster_mean(:,1),cluster_mean(:,2),'kx','LineWidth',2,'MarkerSize',10)
            title(sprintf('SigmaType = %s, SharedCov = %s',Sigma{ii},SCtext{jj}))
            %             legend(h1,{'1','2','3'});
            legend(h1,Legends);
            if c==1 || c==3
                ylabel('$RR_i$')
            end
            
            if c==3 || c==4
                xlabel('$RR_{i-1}$')
            end
            hold off
            c = c + 1;
            
            
        end
    end
end

%% choose best clustering solution
[OD_mat_sorted, ind_sorted] = sort(OD_mat(:))

% [M,I] = min(OD_mat(:)); % best clustering solution

I = ind_sorted(2)
[ind_sigma,ind_sharedCov] = ind2sub(size(OD_mat),I);


chosenGMdist = save_gmfit(ind_sigma,ind_sharedCov);

cluster_means = chosenGMdist{1}.mu;
nclusters = size(cluster_means,1);


% I do this here, because if I put directly in the table the names of the
% variables woudn't show up
options_names = options_names(:);
OD_mat = OD_mat(:);
ICV_mat = ICV_mat(:);
C_mat = C_mat(:);
DI_mat = DI_mat(:);
CS_mat = CS_mat(:);
nclusters_mat = nclusters_mat(:);
nsamples_Larger_cluster_mat = nsamples_Larger_cluster_mat(:);


validity_indices = table(options_names, OD_mat, ICV_mat, C_mat, ...
    DI_mat, CS_mat, nsamples_Larger_cluster_mat, nclusters_mat)

validity_indices_out = validity_indices(I,:)
chosenClusteringSolution = squeeze(save_cluster_solutions(ind_sigma,ind_sharedCov,:));


end