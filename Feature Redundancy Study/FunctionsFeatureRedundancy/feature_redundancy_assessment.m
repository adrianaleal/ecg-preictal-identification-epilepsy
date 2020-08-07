function results_out = feature_redundancy_assessment(feature_dataset, ...
    feat_names2analyse, folder2save, seizure_names, removeSeizures, ...
    chosenDim)

% Inputs:
% - feature_dataset (double): 3D feature dataset of size [n_seiz x n_feat x
%                             n_win]
% - feat_names2analyse (cell): names of the features to analyse
% - folder2save (char): folder to save the results
% - removeSeizures (double): flag to remove seizures that present a high 
%                            number of missing values
% - chosenDim (char): variable that indicates the dimension of the 3D 
%                     feature dataset ('windows' or 'seizures') along which 
%                     to compute redundancy  

% Outputs:
% - results_out (struct): information about the correlated features and the
%                         number of windows/seizures for which redundancy
%                         was verified for each feature combination




close all

if removeSeizures==1 && strcmp(chosenDim, 'windows')
    % this is only necessary when analysing the correlation between vectors
    % with n_seiz length and not for the n_win length
    disp('Removed seizures: ')
    
    indFeat2remove = [71, 191];
    for ii = indFeat2remove
        disp(seizure_names{ii})
    end
    feature_dataset(indFeat2remove,:,:) = [];
end

n_original_feat = numel(feat_names2analyse);
comb = nchoosek(1:n_original_feat, 2);
n_comb = size(comb,1);


[n_seiz, n_feat, n_wins] = size(feature_dataset);

if strcmp(chosenDim, 'seizures')
    dim = n_seiz;
    dimension_string = 'over_seizures';
elseif strcmp(chosenDim, 'windows')
    dim = n_wins;
    dimension_string = 'over_time_windows';
end


%% compute the redundancy measures

variables2compute = {'correlation_between_features_pearson', ...
    'mutual_information_between_features_AMI', 'normalityAssessment'}';
% 'correlation_between_features_spearman', ...
% 'correlation_between_features_kendall', ...
    
n_variables2compute = numel(variables2compute);

computeRedundancy = 1;

subfolder2save = fullfile(cd, folder2save, dimension_string);
if ~exist(subfolder2save, 'dir')
    mkdir(subfolder2save)
end

if computeRedundancy==1
    computeRedundancy(n_feat, n_comb, dim, variables2compute, ...
        dimension_string, feature_dataset, comb, chosenDim, ...
        n_variables2compute, subfolder2save);
else
    for ll = 1:n_variables2compute
        eval(['load(fullfile(subfolder2save, ''' variables2compute{ll} ...
            '.mat''), ''' variables2compute{ll} ''')'])
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
        eval([variables2compute{ll} ' = ' variables2compute{ll} ...
            '(:,~logicFeat2remove);'])
    end
    dim = dim-numel(indFeat2remove);
end

%% observe the redundancy measures matrices

figure()
subplot(211)
imagesc(correlation_between_features_pearson')
cb = colorbar();
ylabel(cb, 'Pearson''s correlation coefficient')
ylabel('Number of 5-min overlapped windows')
xlabel('Index of two-by-two feature combination')
subplot(212)
imagesc(mutual_information_between_features_AMI')
cb = colorbar();
ylabel(cb, 'AMI (bits)')
xlabel('Index of two-by-two feature combination')


% define which correlation to use if more than one is tested:
correlation_between_features = correlation_between_features_pearson;

% define which mutual information to use if more than one is tested:
mutual_information_between_features = mutual_information_between_features_AMI;


%% use the LDA classification to identify the best threshold for each 
% measure

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
        eval(['save(fullfile(cd, subfolder2save, ''' variablesLDA2load{ll} ...
            '''), ''' variablesLDA{ll} ''')'])
    end    
else    
    for ll = 1:n_variablesLDA
        eval(['load(fullfile(cd, subfolder2save, ''' variablesLDA2load{ll} ...
            '''), ''' variablesLDA{ll} ''')'])
    end    
end


% Select the best threshold for each measure on the ROC curve

method_info = {'PCC', ''};
corr_threshold = selectMeasureThreshold(correlation_between_features, ...
    SE_mat_corr, SP_mat_corr, corr_threshold_vec, n_comb, method_info);

% SAVE FIGURE:
% name2save = [method_info{1} '_feature_redundancy_threshold'];
% export_fig(fullfile(cd, subfolder2save, [name2save '.pdf']), '-painters', ...
%     '-transparent')
% saveas(gcf,fullfile(cd, subfolder2save, [name2save '.fig']))
% saveas(gcf,fullfile(cd, subfolder2save, [name2save '.png']))


method_info = {'AMI', ' (bits)'};
mi_threshold = selectMeasureThreshold(mutual_information_between_features, ...
    SE_mat_mi, SP_mat_mi, mi_threshold_vec, n_comb, method_info);

% SAVE FIGURE:
% name2save = [method_info{1} '_feature_redundancy_threshold'];
% export_fig(fullfile(cd, subfolder2save, [name2save '.pdf']), '-painters', ...
%     '-transparent')
% saveas(gcf,fullfile(cd, subfolder2save, [name2save '.fig']))
% saveas(gcf,fullfile(cd, subfolder2save, [name2save '.png']))


%% get the features for which the threshold was overcome

getMeasuresValues = 1;
if getMeasuresValues
    save_correlated_features = [];
    save_dependent_features = [];
    
    for ii = 1:dim
        wins_above_threshold = correlation_between_features(:,ii)>corr_threshold;
        
        ind_features = num2cell(comb(wins_above_threshold,:));
        features = feat_names2analyse(comb(wins_above_threshold,:));
        if sum(wins_above_threshold)==1
            features = features';
        end
        
        save_correlated_features = [save_correlated_features; features, ...
            ind_features];
        
        wins_above_threshold = mutual_information_between_features(:,ii)>mi_threshold;
        ind_features = num2cell(comb(wins_above_threshold,:));
        features = feat_names2analyse(comb(wins_above_threshold,:));
        if sum(wins_above_threshold)==1
            features = features';
        end
        
        save_dependent_features = [save_dependent_features; features, ...
            ind_features];
        
    end
    
    save(fullfile(cd, subfolder2save, ['save_correlated_features_' ...
        dimension_string '.mat']), 'save_correlated_features')
    save(fullfile(cd, subfolder2save, ['save_dependent_features_' ...
        dimension_string '.mat']), 'save_dependent_features')
    
else
    load(fullfile(cd, subfolder2save, ['save_correlated_features_' dimension_string '.mat']), ...
        'save_correlated_features')
    load(fullfile(cd, subfolder2save, ['save_dependent_features_' dimension_string '.mat']), ...
        'save_dependent_features')
end


%% get graph for correlation
featureSelection = 0;

close all
figure()
set(gcf,'units','normalized','outerposition',[0 0 0.51 1])
ax(1) = subplot_tight(2,1,1, [0.04, 0.04]);
method_info = {'PCC', num2str(corr_threshold)};
restructureGraph = 0;
% if we are plotting the graph for the first time, set restructureGraph = 0
% afterwards we can manually change the layout of the nodes and edges

[redundadantFeatures2removePCC, selected_features_namesPCC, ...
    selected_features_indexesPCC] = getFeatureSet(save_correlated_features, ...
    method_info, feat_names2analyse, chosenDim, ax(1), featureSelection, ...
    restructureGraph);

results_out.redundadantFeatures2removePCC = redundadantFeatures2removePCC;
results_out.selected_features_namesPCC = selected_features_namesPCC;
results_out.selected_features_indexesPCC = selected_features_indexesPCC;


% SAVE GRAPH FIGURES:
% name2save = [method_info{1} '_feature_redundancy_graph'];
% export_fig(fullfile(cd, subfolder2save, [name2save '.pdf']), '-painters', ...
%     '-transparent')
% saveas(gcf,fullfile(cd, subfolder2save, [name2save '.fig']))

% SAVE PLOTS CONTAINING THE THRESHOLD FOR EITHER THE NUMBER OF WINDOWS OR
% THE NUMBER OF SEIZURES:
% name2save = [method_info{1} '_feature_redundancy_graph'];
% export_fig(fullfile(cd, subfolder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, subfolder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.fig']))


%% get graph for mutual information


ax(2) = subplot_tight(2,1,2, [0.04, 0.04]);
method_info = {'AMI', num2str(mi_threshold)};
[selected_features_indexesAMI, selected_features_namesAMI, ...
    selected_features_indexes_AMI] = getFeatureSet(save_dependent_features, ...
    method_info, feat_names2analyse, chosenDim, ax(2), featureSelection, ...
    restructureGraph);

results_out.selected_features_indexesAMI = selected_features_indexesAMI;
results_out.selected_features_namesAMI = selected_features_namesAMI;
results_out.selected_features_indexes_AMI = selected_features_indexes_AMI;

% the features here were selected by the user, so that they are in
% accordance with the features selected in correlation

if featureSelection
    indexes_all_feat = 1:n_original_feat;
    selected_features_total_ind = unique([selected_features_indexes_AMI; selected_features_indexesPCC]);
    logical_indexes = ismember((1:n_original_feat),selected_features_total_ind);

    results_out.indexesFeat2Remove = indexes_all_feat(~logical_indexes);
end

% SAVE GRAPH FIGURES:
% name2save = [method_info{1} '_feature_redundancy_graph'];
% export_fig(fullfile(cd, folder2save, [name2save '.pdf']), '-painters', ...
%     '-transparent')
% saveas(gcf,fullfile(cd, folder2save, [name2save '.fig']))

% SAVE PLOTS CONTAINING THE THRESHOLD FOR EITHER THE NUMBER OF WINDOWS OR
% THE NUMBER OF SEIZURES:
% name2save = [method_info{1} '_feature_redundancy_graph'];
% export_fig(fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.pdf']), ...
%     '-painters','-transparent')
% saveas(gcf,fullfile(cd, folder2save, [method_info{1} ...
%     '_feature_selection_overlap_no_interpolation_windows.fig']))



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

function measure_threshold = selectMeasureThreshold(measure, SE_mat, ...
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

% The threshold for each measure corresponds to the first maximum of GM for
% which a ratio, R (variable save_ratio), between the number of samples in 
% the group of data crossing the tested thresholds, th = 0.01, 0.02, 0.03, 
% ..., 1, and the total number of samples, N = 496 x 2768, was higher than 
% 5%.

indexes_R = (save_ratio>=0.05 & (measure_threshold_vec>=0.4)');

geometric_mean_indexes_R = geometric_mean(indexes_R);
[~, I_GM] = max(geometric_mean_indexes_R,[],1);

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