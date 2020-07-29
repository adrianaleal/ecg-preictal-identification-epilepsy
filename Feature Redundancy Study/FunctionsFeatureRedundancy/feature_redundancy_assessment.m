function [results_out, normalityAssessment] = ...
    feature_redundancy_assessment(feature_dataset, ...
    feat_names2analyse, folder2save, names_struct, removeSeizures, ...
    chosenDim)

close all

if removeSeizures==1 && strcmp(chosenDim, 'windows')
    % this is only necessary when analysing the correlation between vectors
    % with n_seiz length and not for the n_win length
    disp('Removed seizures: ')
    
    indFeat2remove = [71, 191];
    for ii = indFeat2remove
        disp(names_struct{ii})
    end
    feature_dataset(indFeat2remove,:,:) = [];
end

n_original_feat = numel(feat_names2analyse);
comb2 = nchoosek(1:n_original_feat, 2);
n_comb2 = size(comb2,1);


[n_seiz, n_feat, n_wins] = size(feature_dataset);

if strcmp(chosenDim, 'seizures')
    dim = n_seiz;
    dimension_string = 'over_seizures';
elseif strcmp(chosenDim, 'windows')
    dim = n_wins;
    dimension_string = 'over_time_windows';
end


variables2compute = {'correlation_between_features_pearson', ...
    'correlation_between_features_spearman', ...
    'correlation_between_features_kendall', ...
    'mutual_information_between_features_AMI', ...
    'mutual_information_between_features_MIGrinstead', ...
    'mutual_information_between_features_MIMatlab', ...
    'mutual_information_between_features_MI', ...
    'mutual_information_between_features_NMI', 'normalityAssessment'}';
variables2load = strcat(variables2compute, ['_' dimension_string '.mat']);
n_variables2compute = numel(variables2compute);

computeFromScratch = 0;

if computeFromScratch==1
    %save_correlated_features = [];
    normalityAssessment = zeros(n_feat, dim);
    pNormalityAssessment = normalityAssessment;
    
    
    % % test if correlation
    % load('NaN_seizures_HF_POWER.mat')
    % load('NaN_seizures_LF_POWER.mat')
    % load('less_than_300s_seizures.mat')
    
    % ver = [NaN_seizures_LF_POWER, NaN_seizures_HF_POWER];
    % NaN_vec = (NaN_seizures_HF_POWER<1 & NaN_seizures_LF_POWER<1);
    %
    % indexes_discarded_features = load('indexes_seizures_discarded_by_MVs.mat')
    % NaN_vec([171,175]) = [];
    
    % less_than_300s_vec = less_than_300s_seizures<1;
    % less_than_300s_vec([171,175]) = [];
    
    
    % initialize variables2compute
    eval([variables2compute{1} ' = zeros(n_comb2, dim);'])
    for ll = 2:n_variables2compute
        eval([variables2compute{ll} ' = ' variables2compute{1} ';'])
    end
    
    for ii = 1:dim

        disp([dimension_string ' ' num2str(ii)])

        %% Compute correlation between features (Pearson's)
        
        
        if strcmp(chosenDim, 'seizures')
            data2compare = squeeze(feature_dataset(ii,:,:))';
        elseif strcmp(chosenDim, 'windows')
            data2compare = squeeze(feature_dataset(:,:,ii))';
        end
        
        % data_each_win = data_each_win(less_than_300s_vec,:);
        if any(strcmp(variables2compute, 'correlation_between_features_pearson'))
            correlation_between_features_pearson(:,ii) = corr_features(data2compare);
        end
        
        %% Check feature normality
        for ff = 1:n_feat
            [normalityAssessment(ff,ii), pNormalityAssessment(ff,ii)] = kstest(data2compare(:,ff));
        end
        % returns a test decision for the null hypothesis that the data in vector
        % x comes from a standard normal distribution, against the alternative that
        % it does not come from such a distribution, using the one-sample
        % Kolmogorov-Smirnov test. The result h is 1 if the test rejects the null
        % hypothesis at the 5% significance level, or 0 otherwise.
        
        
        %% Compute correlation and mutual information between features
        
        for jj = 1:n_comb2
            getMat2correlate = data2compare(:,comb2(jj,:));
            getMat2correlateNoNaN = getMat2correlate(~sum(isnan(getMat2correlate),2)>0,:);
            
            % correlation *************************************************
            if any(strcmp(variables2compute, 'correlation_between_features_spearman'))
                [RHO,PVAL] = corr(getMat2correlateNoNaN(:,1), getMat2correlateNoNaN(:,2), ...
                    'Type','Spearman');
                correlation_between_features_spearman(jj,ii) = RHO;
            end
            
            if any(strcmp(variables2compute, 'correlation_between_features_kendall'))
                [RHO,PVAL] = corr(getMat2correlateNoNaN(:,1), getMat2correlateNoNaN(:,2), ...
                    'Type','Kendall');
                correlation_between_features_kendall(jj,ii) = RHO;
            end
            % *************************************************************
            
            
            % mutual information ******************************************
            if any(strcmp(variables2compute, 'mutual_information_between_features_MI')) || ...
                    any(strcmp(variables2compute, 'mutual_information_between_features_NMI'))
                [MI, NMI, HXY] = MutualInformation(getMat2correlateNoNaN(:,1), ...
                    getMat2correlateNoNaN(:,2));
                
                if any(strcmp(variables2compute, 'mutual_information_between_features_NMI'))
                    mutual_information_between_features_NMI(jj,ii) = NMI;
                end
                
                if any(strcmp(variables2compute, 'mutual_information_between_features_MI'))
                    mutual_information_between_features_MI(jj,ii) = MI;
                end
            end
            
            if any(strcmp(variables2compute, 'mutual_information_between_features_MIMatlab'))
                % MI matlab
                arrLag = 1:10;
                minfo = zeros(size(arrLag));
                for kk = 1:length(arrLag)
                    tmp = predmaint.internal.NonlinearFeatures.getPhaseSpace(getMat2correlateNoNaN, arrLag(kk), 2);
                    minfo(kk) = predmaint.internal.NonlinearFeatures.getMutualInfo(tmp(:,1), tmp(:,2), 10);
                end
                MIMatlab = minfo(find(minfo<(minfo(1)/exp(1)), 1));
                if isempty(MIMatlab)
                    MIMatlab = min(minfo);
                end
                mutual_information_between_features_MIMatlab(jj,ii) = MIMatlab;
            end
            
            
            if any(strcmp(variables2compute, 'mutual_information_between_features_MIMatlab'))
                % MI Grinstead
                [MIGrinstead,lag] = amiGrinstead(getMat2correlateNoNaN(:,1), ...
                    getMat2correlateNoNaN(:,2));
                mutual_information_between_features_MIGrinstead(jj,ii) = MIGrinstead;
            end
            
            if any(strcmp(variables2compute, 'mutual_information_between_features_AMI'))
                plotFigure = 0;
                mutual_information_between_features_AMI(jj,ii) = average_mutual_information( ...
                    getMat2correlateNoNaN, plotFigure);
            end
            
        end
        
    end
    
    for ll = 1:n_variables2compute
        eval(['save(fullfile(cd, folder2save, ''' variables2load{ll} ...
            '''), ''' variables2compute{ll} ''')'])
    end

else
    
    for ll = 1:n_variables2compute
        eval(['load(fullfile(cd, folder2save, ''' variables2load{ll} ...
            '''), ''' variables2compute{ll} ''')'])
    end
    
end


if removeSeizures==1 || strcmp(chosenDim, 'seizures')
    % this is only necessary when analysing the correlation between vectors
    % with n_seiz length and not for the n_win length
    disp('Removed seizures: ')
    
    indFeat2remove = [71, 191];
    feat2remove = zeros(1,dim);
    feat2remove(indFeat2remove) = 1; 
    logicFeat2remove = logical(feat2remove);
    for ll = 1:n_variables2compute-1
        eval([variables2compute{ll} ' = ' variables2compute{ll} '(:,~logicFeat2remove);'])
    end
    dim = dim-numel(indFeat2remove);
end

%% observe the matrices obtained for correlation and the different MI
% methods

figure()
subplot(221)
imagesc(correlation_between_features_pearson')
cb = colorbar();
ylabel(cb, 'Pearson''s correlation coefficient')
ylabel('Number of 5-min overlapped windows')
xlabel('Index of two-by-two feature combination')
subplot(222)
imagesc(correlation_between_features_spearman')
cb = colorbar();
ylabel(cb, 'Spearman''s correlation coefficient')
xlabel('Index of two-by-two feature combination')
subplot(223)
imagesc(mutual_information_between_features_MIMatlab')
cb = colorbar();
ylabel(cb, 'Kendall''s correlation coefficient')
ylabel('Number of 5-min overlapped windows')
xlabel('Index of two-by-two feature combination')
subplot(224)
imagesc(mutual_information_between_features_AMI')
cb = colorbar();
ylabel(cb, 'AMI (bits)')
xlabel('Index of two-by-two feature combination')


%% observe the matrices obtained for correlation and the different MI
% methods



figure()
subplot(221)
imagesc(mutual_information_between_features_AMI')
cb = colorbar();
ylabel(cb, 'AMI (bits)')
ylabel('Number of 5-min overlapped windows')
xlabel('Index of two-by-two feature combination')
subplot(222)
imagesc(correlation_between_features_pearson')
cb = colorbar();
ylabel(cb, 'Correlation')
xlabel('Index of two-by-two feature combination')
subplot(223)
imagesc(mutual_information_between_features_MIMatlab')
cb = colorbar();
ylabel(cb, 'MI matlab (bits)')
ylabel('Number of 5-min overlapped windows')
xlabel('Index of two-by-two feature combination')
subplot(224)
imagesc(mutual_information_between_features_MIGrinstead')
cb = colorbar();
ylabel(cb, 'MI Grinstead (bits)')
xlabel('Index of two-by-two feature combination')

%% define which correlation to use:
correlation_between_features = correlation_between_features_spearman;

% define which mutual information to use:
mutual_information_between_features = mutual_information_between_features_AMI;


%% get the LDA classification in order to identify the best threshold

computeLDAclassification = 1;
variablesLDA = {'SE_mat_corr', 'SP_mat_corr', 'corr_threshold_vec', ...
    'SE_mat_mi', 'SP_mat_mi', 'mi_threshold_vec'};
variablesLDA2load = strcat(variablesLDA, ['_' dimension_string '.mat']);
n_variablesLDA = numel(variablesLDA);
    
if computeLDAclassification
    
    % correlation
    [SE_mat_corr, SP_mat_corr, corr_threshold_vec] = ...
        getLDAclassification(correlation_between_features);
    
    % mutual information
    [SE_mat_mi, SP_mat_mi, mi_threshold_vec] = ...
        getLDAclassification(mutual_information_between_features);
    
    for ll = 1:n_variablesLDA
        eval(['save(fullfile(cd, folder2save, ''' variablesLDA2load{ll} ...
            '''), ''' variablesLDA{ll} ''')'])
    end
    
else
    
    for ll = 1:n_variablesLDA
        eval(['load(fullfile(cd, folder2save, ''' variablesLDA2load{ll} ...
            '''), ''' variablesLDA{ll} ''')'])
    end
end



%% Select the best threshold on the ROC curve
method_info = {'PCC', ''};
corr_threshold = selectMeasureThreshold(correlation_between_features, ...
    SE_mat_corr, SP_mat_corr, corr_threshold_vec, n_comb2, method_info);
% export_fig(fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_threshold.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_threshold.fig']))
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_threshold.png']))

method_info = {'AMI', ' (bits)'};
mi_threshold = selectMeasureThreshold(mutual_information_between_features, ...
    SE_mat_mi, SP_mat_mi, mi_threshold_vec, n_comb2, method_info);
% mi_threshold = 1;
% export_fig(fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_threshold.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_threshold.fig']))
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_threshold.png']))

%% get the features for which the threshold was overcome

getMeasuresValues = 1;
if getMeasuresValues
    save_correlated_features = [];
    save_dependent_features = [];
    
    for ii = 1:dim
        wins_above_threshold = correlation_between_features(:,ii)>corr_threshold;
        
        ind_features = num2cell(comb2(wins_above_threshold,:));
        features = feat_names2analyse(comb2(wins_above_threshold,:));
        if sum(wins_above_threshold)==1
            features = features';
        end
        
        save_correlated_features = [save_correlated_features; features, ...
            ind_features];
        
        wins_above_threshold = mutual_information_between_features(:,ii)>mi_threshold;
        ind_features = num2cell(comb2(wins_above_threshold,:));
        features = feat_names2analyse(comb2(wins_above_threshold,:));
        if sum(wins_above_threshold)==1
            features = features';
        end
        
        save_dependent_features = [save_dependent_features; features, ...
            ind_features];
        
    end
    
    save(fullfile(cd, folder2save, ['save_correlated_features_' dimension_string '.mat']), ...
        'save_correlated_features')
    save(fullfile(cd, folder2save, ['save_dependent_features_' dimension_string '.mat']), ...
        'save_dependent_features')
    
else
    load(fullfile(cd, folder2save, ['save_correlated_features_' dimension_string '.mat']), ...
        'save_correlated_features')
    load(fullfile(cd, folder2save, ['save_dependent_features_' dimension_string '.mat']), ...
        'save_dependent_features')
end


%% get features for high correlation
featureSelection = 0

close all
figure()
set(gcf,'units','normalized','outerposition',[0 0 0.51 1])
ax(1) = subplot_tight(2,1,1, [0.04, 0.04]);
method_info = {'PCC', num2str(corr_threshold)};
[correlatedFeatures2remove, selected_features_namesCorr, ...
    selected_features_indexesCorr] = getFeatureSet(save_correlated_features, ...
    method_info, feat_names2analyse, chosenDim, ax(1), featureSelection);

results_out.correlatedFeatures2remove = correlatedFeatures2remove;
results_out.selected_features_namesCorr = selected_features_namesCorr;
results_out.selected_features_indexesCorr = selected_features_indexesCorr;

% export_fig(fullfile(cd, folder2save, [method_info{1}, ...
%     '_feature_selection_overlap_no_interpolation.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation.fig']))
%
% export_fig(fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.fig']))

xlabel('\textbf{(a)}')

%% get features for high mutual information


ax(2) = subplot_tight(2,1,2, [0.04, 0.04]);
method_info = {'AMI', num2str(mi_threshold)};
[features2removeMI, selected_features_namesMI, ...
    selected_features_indexes_MI] = getFeatureSet(save_dependent_features, ...
    method_info, feat_names2analyse, chosenDim, ax(2), featureSelection);

results_out.features2removeMI = features2removeMI;
results_out.selected_features_namesMI = selected_features_namesMI;
results_out.selected_features_indexes_MI = selected_features_indexes_MI;

% the features here were selected by the user, so that they are in
% accordance with the features selected in correlation


indexes_all_feat = 1:n_original_feat;
selected_features_total_ind = unique([selected_features_indexes_MI; selected_features_indexesCorr]);
logical_indexes = ismember((1:n_original_feat),selected_features_total_ind);

results_out.indexesFeat2Remove = indexes_all_feat(~logical_indexes);


% export_fig(fullfile(cd, folder2save, [method_info{1}, ...
%     '_feature_selection_overlap_no_interpolation.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation.fig']))

% export_fig(fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.fig']))

xlabel('\textbf{(b)}')
% export_fig(fullfile(cd, folder2save, ...
%     ['feature_selection_overlap_no_interpolation_' dimension_string '.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, ...
%     ['feature_selection_overlap_no_interpolation_' dimension_string '.fig']))

%% evaluate behaviour of correlated features over time

% n_seizures = size(feature_dataset,1);
%
% for ii = 1:n_seizures
%
%     for jj = 1:size(correlatedFeatures2remove,1)
%
%         % visualy inspect the relationship among the features
%
%         c = time_dataset(ii,:);
%         feat1_data = feature_dataset(ii,cell2mat(correlatedFeatures2remove(jj,2)),:);
%         feat2_data = feature_dataset(ii,cell2mat(correlatedFeatures2remove(jj,2)),:);
%         figure(1)
%         scatter(feat1_data,feat2_data,[],c,'*') % só o canal 1
%         xlabel(correlatedFeatures2remove{jj,1})
%         ylabel(correlatedFeatures2remove{jj,3})
%         title(['Seizure ' num2str(ii)])
%         pause
%     end
%
% end


end

function [SE_mat, SP_mat, measure_level_vec] = getLDAclassification(measure)

max_feature = max(max(measure));
step2increment = 0.01;
measure_level_vec = step2increment:step2increment:max_feature-step2increment;

SP_mat = zeros(length(measure_level_vec),1);
SE_mat = SP_mat;

for ii = 1:length(measure_level_vec)
    mutual_information_group1 = measure>measure_level_vec(ii);
    Xtrain = measure(:);
    Ytrain = mutual_information_group1(:);
    MdlLinearLDA = fitcdiscr(Xtrain, Ytrain);
    [ypredTrain, scores] = predict(MdlLinearLDA,Xtrain);
    [SE_mat(ii), SP_mat(ii)] = Validation(Ytrain, ypredTrain);
end



end

function [measure_threshold] = selectMeasureThreshold(measure, SE_mat, ...
    SP_mat, measure_threshold_vec, n_comb2, method_info)


% roc curve
FPR = 1-SP_mat;
geometric_mean = sqrt(SE_mat.^2 + SP_mat.^2);





see = [SE_mat SP_mat geometric_mean measure_threshold_vec'];

save_ratio = zeros(numel(measure_threshold_vec),1);
for ii = 1:numel(measure_threshold_vec)
    mutual_information_group1 = measure>=measure_threshold_vec(ii);
    save_ratio(ii) = sum(sum(mutual_information_group1))/numel(measure);
end

indexes_R = (save_ratio>=0.069 & (measure_threshold_vec>=0.4)');

geometric_mean_indexes_R = geometric_mean(indexes_R);
[~,I_GM] = max(geometric_mean_indexes_R,[],1);

measure_threshold_vec_indexsR = measure_threshold_vec(indexes_R);
measure_threshold = measure_threshold_vec_indexsR(I_GM);

figure()
set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
subplot(321)
plot(FPR, SE_mat, '*')
xlabel('1-SP')
ylabel('SE')
axis tight
subplot(322)
yyaxis left
plot(measure_threshold_vec, geometric_mean), hold on
h2 = plot(measure_threshold_vec_indexsR(I_GM), geometric_mean_indexes_R(I_GM),'ok');
ylabel('GM')
hold off
ylim([min(geometric_mean) max(geometric_mean)+0.01])
yyaxis right
plot(measure_threshold_vec, save_ratio), hold on
save_ratio_indexes_R = save_ratio(indexes_R);
h2 = plot(measure_threshold_vec_indexsR(I_GM), save_ratio_indexes_R(I_GM),'ok');
hold off
ylabel('R')
legend(h2, [method_info{1} '$_{threshold}$'],'Location', 'Best')
xlabel('th')

axis tight

ax2 = gca;
subplot(3,2,3:4)
imagesc(measure')
% xlabel('Index of two-by-two feature combination')
ylabel('5-min window')
cb = colorbar;
ylabel(cb, [method_info{1} method_info{2} ], 'Interpreter', 'Latex')
pos = get(cb,'Position');
pos(1) = ax2.Position(1)+ax2.Position(3)+0.01;
set(cb,'Position',pos);
set(cb, 'TickLabelInterpreter', 'latex')

subplot(3,2,5:6)

plot(1:n_comb2, measure);
axis tight
xlabel('Index of 2-combinations of features')
ylabel(method_info{1})

hold on
% h1 = plot(1:n_comb2,mean(correlation_between_features,2),'k', ...
%     'LineWidth',2);

h2 = plot(1:n_comb2,ones(n_comb2,1)*measure_threshold,'k--', ...
    'LineWidth',2);
legend(h2, [method_info{1} '$_{threshold}$'])


% measure_group1 = measure(measure>measure_threshold);
% measure_group2 = measure(measure<=measure_threshold);
% plot(1:numel(measure_group1), measure_group1(:))
% hold on
% plot(1:numel(measure_group2), measure_group2(:), 'r')
hold off
axis tight


end