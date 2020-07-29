fclose all; clear; close all; clc;
set(groot, 'DefaultAxesFontSize',11,'defaultTextFontSize',11)
set(groot, 'DefaultAxesTickLabelInterpreter','tex');
set(groot, 'DefaultLegendInterpreter','tex');
set(groot, 'DefaultTextInterpreter','tex')


cd ..
outer_folder_path = cd;
load('patients_info_final_preictal_study_240min_before_seiz')
cd PrepareFeatureDataset2Clustering
    
results_folder_path = fullfile(cd,'ResultsOverlap'); % Annotated
if ~exist(results_folder_path, 'dir')
    mkdir(results_folder_path)
end

if ~isempty(strfind(results_folder_path,'Annotated'))
    feat_folder = 'ExtractedFeaturesNoEctopicBeatsOverlapAnnotated';
    windows_folder = 'OverlappedWindowsAnnotated';
else
    feat_folder = 'ExtractedFeaturesNoEctopicBeatsOverlap';
    windows_folder = 'OverlappedWindows';
end

computeFromScratch = 0;

if computeFromScratch==1
    
    patients_folder_path = fullfile(outer_folder_path, ...
        'ResultsEPILEPSIAE', feat_folder ,'pat_*');
    pats_name = getPatientsNames(patients_folder_path);
    
    mat_number_seizures = [pats_name, cell(numel(pats_name),1)];
    clear patients split1

    names_struct = [];
    time_vec_original = [];
    percent_noise_original = [];
    comp_time_total = [];
    n_5min_wins_per_seizure = [];
    count_seizures = 0;
    
    for pp = 1:length(pats_name)
        file_name = pats_name{pp}
        seiz_folders = dir(fullfile(outer_folder_path,'ResultsEPILEPSIAE',feat_folder, file_name));
        seiz_folders = seiz_folders(~ismember({seiz_folders.name},{'.','..'}));
        num_seizures_pat = numel(seiz_folders);
        mat_number_seizures(pp,2) = num2cell(num_seizures_pat);
        for ss = 1:numel(seiz_folders)
            count_seizures = count_seizures+1;
            
            comp_time = load(fullfile(seiz_folders(ss).folder, ...
                seiz_folders(ss).name, [file_name '_seizure' num2str(ss) '_comp_time']));
            comp = fieldnames(comp_time);
            comp_time_total = horzcat(comp_time_total, comp_time.(comp{1}));
            
            
            % load hrv features
            load(fullfile(seiz_folders(ss).folder,seiz_folders(ss).name, ...
                [file_name '_seizure' num2str(ss)]))
            feat_names = fieldnames(hrv_features);
            
            n_5min_wins_per_seizure = [n_5min_wins_per_seizure; numel(hrv_features)];
            time_vec_original = [time_vec_original; {vertcat(hrv_features.time)}];
            percent_noise_original = [percent_noise_original; {vertcat(hrv_features.percent_detected_noise)}];
            
            if pp==1 && ss==1
                for ff = 1:numel(feat_names)
                    feature_dataset_240min_before_seizure.(feat_names{ff}) = vertcat(hrv_features.(feat_names{ff}));
                end
                min_samples_before = numel(hrv_features);
            else
                min_samples_after = numel(hrv_features);
                if min_samples_before>min_samples_after
                    min_samples = min_samples_after;
                else
                    min_samples = min_samples_before;
                end
                min_samples_before = min_samples;
                
                for ff = 1:numel(feat_names)
                    new_feat_values = vertcat(hrv_features.(feat_names{ff}))';
                    old_feat_values = feature_dataset_240min_before_seizure.(feat_names{ff});
                    if size(old_feat_values,2)==1
                        old_feat_values = old_feat_values';
                    end
                    % the values of the features are removed from the
                    % beggining and not from seizure neighbourhood
                    feature_dataset_240min_before_seizure.(feat_names{ff}) ...
                        = [old_feat_values(:,end-min_samples+1:end); ...
                        new_feat_values(end-min_samples+1:end)];
                end
            end
        end
        
        names_struct = [names_struct; strcat([file_name '_' num2str(pp) ...
            '_seizure_'], cellstr(num2str((1:num_seizures_pat)')))];
    end
    clear num_seizures_pat hrv_features count_seizures new_feat_values ...
        min_samples old_feat_values seiz_folders file_name ff ss pp ...
        min_samples_after min_samples_before comp comp_time
    
    feature_dataset_240min_before_seizure.names = names_struct;
    
    n_seizures = vertcat(patients_info_final.n_seizures)- ...
        vertcat(patients_info_final.n_discarded_seizures_by_240min)- ...
        vertcat(patients_info_final.n_discarded_seizures_by_no_file);
    check_n_seizures_patients_info_final = sum(n_seizures)
    check_n_seizures_each_seizure = sum(cell2mat(mat_number_seizures(:,2)))
    clear mat_number_seizures
    
    % put the invalid windows to NaN *************************************
    % there are already some segment_size assigned to NaN...
    
    feature_dataset_240min_before_seizure = editInvalidWindows(feat_names, ...
        feature_dataset_240min_before_seizure);
    
    % get the 3D matrix to perform feature selection and clustering ******
    [n_seizures, n_wins] = size(feature_dataset_240min_before_seizure.nHB);
    get_feat_names = fieldnames(feature_dataset_240min_before_seizure);
    featsNot2analyse = {'long_window', 'time', 'segment_size_seconds', ...
        'n_clean_segments', 'percent_detected_noise', 'nHB', 'VLF_NORM', ...
        'tau_rec_phase_space', 'embDim_rec_phase_space', 'names', 'abnormalRR'};
    indFeatsNot2analyse = ismember(get_feat_names,featsNot2analyse);
    feat_names2analyse = get_feat_names(~indFeatsNot2analyse);
    
    
    % compute 3D matrix:
    n_features = numel(feat_names2analyse);
    feature_dataset_240min_before_seizure_3D = zeros(n_seizures, n_wins, n_features);
    for ff = 1:n_features
        feature_dataset_240min_before_seizure_3D(:,:,ff) = ...
            feature_dataset_240min_before_seizure.(feat_names2analyse{ff});
    end
    clear ff
    
    % change matrix dimensions to [n_seizures, n_feat, n_wins]
    feature_dataset_240min_before_seizure_3D = permute(feature_dataset_240min_before_seizure_3D,[1,3,2]);
    
    % save 3D matrix
    save(fullfile(results_folder_path, 'feature_dataset_240min_before_seizure_3D.mat'),...
        'feature_dataset_240min_before_seizure_3D')
    
    % save original structure with the features
    feature_dataset_240min_before_seizure_all_feat = feature_dataset_240min_before_seizure;
    save(fullfile(results_folder_path, 'feature_dataset_240min_before_seizure_all_feat.mat'), ...
        'feature_dataset_240min_before_seizure_all_feat')
    

    % remove features that are not informative:
    feature_dataset_240min_before_seizure = rmfield(feature_dataset_240min_before_seizure, ...
        get_feat_names(indFeatsNot2analyse));
    
    % get feature extraction computational time 
    all_values = squeeze(cell2mat(struct2cell(comp_time_total)))';
    feats = [''; fieldnames(comp_time_total)];
    ct_mean = nanmean(all_values)';
    ct_std = nanstd(all_values)';
    compTimeTable = table(feats, ct_mean, ct_std);
    clear comp_time_total ct_mean ct_std feats all_values
    
    n_5min_wins_per_seizure = [names_struct, num2cell(n_5min_wins_per_seizure)];
    
    structure2selectFeatPreIctalStudy.n_5min_wins_per_seizure = n_5min_wins_per_seizure;
    structure2selectFeatPreIctalStudy.percent_noise_original = percent_noise_original;
    structure2selectFeatPreIctalStudy.time_vec_original = time_vec_original;
    structure2selectFeatPreIctalStudy.nwins = n_wins;
    structure2selectFeatPreIctalStudy.names_struct = names_struct;
    structure2selectFeatPreIctalStudy.feature_dataset_240min_before_seizure_3D = feature_dataset_240min_before_seizure_3D;
    structure2selectFeatPreIctalStudy.feat_names2analyse = feat_names2analyse;
    structure2selectFeatPreIctalStudy.feature_dataset_240min_before_seizure = feature_dataset_240min_before_seizure;
    structure2selectFeatPreIctalStudy.n_seizures = n_seizures;
    structure2selectFeatPreIctalStudy.compTimeTable = compTimeTable;
    
    % time_dataset = feature_dataset_240min_before_seizure.time;
    % structure2selectFeatPreIctalStudy.time_dataset = time_dataset;
    save(fullfile(results_folder_path, 'structure2selectFeatPreIctalStudy.mat'), ...
        'structure2selectFeatPreIctalStudy')
    
    clear get_feat_names indFeatsNot2analyse feat_names featsNot2analyse I
    
else
    load(fullfile(results_folder_path, 'feature_dataset_240min_before_seizure_all_feat.mat'))
    load(fullfile(results_folder_path, 'structure2selectFeatPreIctalStudy.mat'))
    feature_dataset_240min_before_seizure_3D = structure2selectFeatPreIctalStudy.feature_dataset_240min_before_seizure_3D;
    feature_dataset_240min_before_seizure = structure2selectFeatPreIctalStudy.feature_dataset_240min_before_seizure;
    feat_names2analyse = structure2selectFeatPreIctalStudy.feat_names2analyse;
    names_struct = structure2selectFeatPreIctalStudy.names_struct;
    n_wins = structure2selectFeatPreIctalStudy.nwins;
    time_vec_original = structure2selectFeatPreIctalStudy.time_vec_original;
    percent_noise_original = structure2selectFeatPreIctalStudy.percent_noise_original;
    n_seizures = structure2selectFeatPreIctalStudy.n_seizures;
    compTimeTable = structure2selectFeatPreIctalStudy.compTimeTable;
    n_5min_wins_per_seizure = structure2selectFeatPreIctalStudy.n_5min_wins_per_seizure;
    % time_dataset = structure2clusterPreIctalStudy.time_dataset;
    clear structure2selectFeatPreIctalStudy
end


% get the patients' names to initialize the structure to save the results
C = cellfun(@(x)strsplit(x, '_' ), names_struct(:,1), 'UniformOutput', false);
names_struct_separated = vertcat(C{:});
clear C
pats_name = unique(names_struct_separated(:,2),'stable');
n_pat = numel(pats_name);

% to perform clustering, clustering solution evaluation and the other
% analysis: ***************************************************************
feat_comb = 3;
if ~isempty(strfind(results_folder_path,'Annotated'))
    clustering_folder_path = fullfile(cd, ...
        'ResultsPreIctalStudyClusteringOverlapAnnotated', ...
        ['ResultsPreICtalStudyClusteringFeatComb' num2str(feat_comb) ...
        'ByPatientSeizure']);
else
    clustering_folder_path = fullfile(cd, ...
        'ResultsPreIctalStudyClusteringOverlap', ...
        ['ResultsPreICtalStudyClusteringFeatComb' num2str(feat_comb) ...
        'ByPatientSeizure']);
end
% ************************************************************************  



% get patient metadata ***************************************************
% run the following function in order to run function getSaveInfoPatients.m

[var2correlate, matrix_patient_info_classes, ...
    n_seiz, matrix_patient_info, mat4AssociationLearning, metaData, ...
    store_seizures_start_date] = ...
    getMetaDataPreIctal(results_folder_path, patients_info_final);

% get time intervals between seizures
getBetweenSeizureIntervals(matrix_patient_info)

% In MEGAsync\Doutoramento\Tese\Papers Literatura\EPILEPSIAE
% from PatientInformationWithLocalisation 
n_patients_with_localization = 263;
% from PatientInformationLocalisationtemporal
percent_patients_with_temporal_lobe = round(215/n_patients_with_localization*100*100)/100
% from PatientInformationLocalisationOnlytemporal
percent_patients_only_temporal_lobe = round(174/n_patients_with_localization*100*100)/100




%% check for segments less than 300 seconds
AbnormalRR_seizures_stat_cell = checkLessThan300secsSegments(n_seizures, ...
    patients_info_final, feature_dataset_240min_before_seizure_all_feat, ...
    results_folder_path, pats_name, names_struct_separated);

%% Interpolate NaNs in matrix in order to perform feature selection
% Interpolate over time
perform_interpolation = 1;

if perform_interpolation==1
    time_dataset = [];
    [feature_dataset_240min_before_seizure_3D_interpolated, ...
        feature_dataset_240min_before_seizure_3D_NaN, ...
        feature_dataset_240min_before_seizure_interpolated, time_dataset, ...
        tab_NaN_stat] = interpolate_over_time(feature_dataset_240min_before_seizure_3D, ...
        feature_dataset_240min_before_seizure, feat_names2analyse, ...
        time_dataset, results_folder_path);
    
    structureInterpolatedPreIctalStudy.feature_dataset_240min_before_seizure_3D = feature_dataset_240min_before_seizure_3D_interpolated;
    structureInterpolatedPreIctalStudy.feat_names2analyse = feat_names2analyse;
    structureInterpolatedPreIctalStudy.feature_dataset_240min_before_seizure = feature_dataset_240min_before_seizure_interpolated;
    save(fullfile(results_folder_path, 'structureInterpolatedPreIctalStudy.mat'), ...
        'structureInterpolatedPreIctalStudy')
else
    load(fullfile(results_folder_path, 'structureInterpolatedPreIctalStudy.mat'))
    feature_dataset_240min_before_seizure_3D_interpolated = structureInterpolatedPreIctalStudy.feature_dataset_240min_before_seizure_3D;
    feature_dataset_240min_before_seizure_interpolated = structureInterpolatedPreIctalStudy.feature_dataset_240min_before_seizure;
    feat_names2analyse = structureInterpolatedPreIctalStudy.feat_names2analyse;
    clear structureInterpolatedPreIctalStudy
end

close all

% see clustergram by seizure (?)
clustergram(squeeze(feature_dataset_240min_before_seizure_3D_interpolated(:,1,:)))

%% Perform feature selection before applying clustering methods
% By feature selection one means assess which features are correlated and
% also the mutual information
% does not work when data contains missig values

removeSeizures = 1;
if removeSeizures==1
    folder2save = 'ResultsFeatureSelectionRemovedSeizures';
else
    folder2save = 'ResultsFeatureSelection';
end

folder2savePath = fullfile(cd, folder2save);
if ~exist(folder2savePath, 'dir')
    mkdir(folder2savePath)
end

functionsFolder = 'FunctionsFeatureSelection';
folder2savePath = fullfile(cd, functionsFolder);
if exist(folder2savePath, 'dir')
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder2savePath));
end

% [results_feat_selection, normalityAssessment] = feature_redundancy_assessment_over_ime(...
%     feature_dataset_240min_before_seizure_3D, feat_names2analyse, ...
%     folder2save, names_struct_separated(:,2), removeSeizures);

chosenDim = 'windows' % seizures
[results_feat_selection, normalityAssessment] = feature_redundancy_assessment(...
    feature_dataset_240min_before_seizure_3D, feat_names2analyse, ...
    folder2save, names_struct_separated(:,2), removeSeizures, chosenDim);

% By removing seizures 4 and 2 of patients 53402 and 104602, respectively,
% no change was seen in the feature selection results, except for the
% number of windows verifying a given value of PCC and AMI

if exist(folder2savePath, 'dir')
    % Remove that folder plus all subfolders to the path.
    rmpath(genpath(folder2savePath))
end


%%


close all
folder2save = 'ResultsFeatureSelection';
matrix_patient_info_VS = matrix_patient_info(:,3);
matrix_patient_info_ST = matrix_patient_info(:,7);
newStr = split(matrix_patient_info(:,2), ' ');
matrix_patient_info_Hour = newStr(:,2);

variables_names = {'VS'; 'ST'; 'Hour'};
matrix_patient_info_final = [matrix_patient_info_VS, matrix_patient_info_ST, ...
    matrix_patient_info(:,8)];

for vv = 1:3

    [matrix_patient_info_sorted, indexes] = sort(matrix_patient_info_final(:,vv));
    
    computeStats = 0;
    if computeStats
        feature_dataset_240min_before_seizure_3D_sorted = feature_dataset_240min_before_seizure_3D;
        for ff = 1:numel(feat_names2analyse)
            feature_dataset_240min_before_seizure_3D_sorted(:,ff,:) = ...
                feature_dataset_240min_before_seizure_3D(indexes,ff,:);
        end
        
        [statisticalDifferences_sorted_by.(variables_names{vv}), ...
            pStatisticalDifferences_sorted_by.(variables_names{vv})] = ...
            assess_seizure_differences(feature_dataset_240min_before_seizure_3D_sorted, ...
            feat_names2analyse);
        save(fullfile(folder2save,'statisticalDifferences_sorted_by.mat'), ...
            'statisticalDifferences_sorted_by')
        save(fullfile(folder2save,'pStatisticalDifferences_sorted_by.mat'), ...
            'pStatisticalDifferences_sorted_by')
    else
        load(fullfile(folder2save,'statisticalDifferences_sorted_by.mat'), ...
            'statisticalDifferences_sorted_by')
        load(fullfile(folder2save,'pStatisticalDifferences_sorted_by.mat'), ...
            'pStatisticalDifferences_sorted_by')
    end
    
    
    for ff = 1:numel(feat_names2analyse)
        statisticalDifferences_sorted = statisticalDifferences_sorted_by.(variables_names{vv}).(feat_names2analyse{ff});
        pStatisticalDifferences_sorted = pStatisticalDifferences_sorted_by.(variables_names{vv}).(feat_names2analyse{ff});
        
        % sort by vigilance state
        figure(22)
        statisticalDifferences_sorted_inverted = flipud(fliplr(statisticalDifferences_sorted));
        imagesc(statisticalDifferences_sorted)
        colormap(gca, pink(2));
        ncolors = 2;
        cb = colorbar('YTick',0.5*(ncolors-1)/ncolors:(ncolors-1)/ncolors:ncolors,...
            'YTickLabel',{'h = 0 (same population)', 'h = 1 (not same population)'}, ...
            'TickLabelInterpreter','latex');
        
        [var_classes, ~, ind2] = unique(matrix_patient_info_sorted);
        n_classes = numel(var_classes);
        ax1 = gca;
        xticks_vec = zeros(n_classes,1);
        for pp = 1:n_classes
            ind_seizures_pat = find(strcmp(var_classes{pp},matrix_patient_info_sorted));
            n_seizures_pat = numel(ind_seizures_pat);
            ind = round(n_seizures_pat/2);
            if bitget(n_seizures_pat,1) %odd
                val = ind_seizures_pat(ind);
            else %even
                val = ind_seizures_pat(ind)+0.5;
            end
            
            xticks_vec(pp) = val;
            
            line([ind_seizures_pat(end)+0.5,ind_seizures_pat(end)+0.5], ...
                ax1.YLim, 'Color', 'k');
            
            line(ax1.XLim, [ind_seizures_pat(end)+0.5,ind_seizures_pat(end)+0.5], ...
                'Color', 'k');
        end
        xticks(xticks_vec)
        xticklabels(var_classes)
        yticks(xticks_vec)
        yticklabels(var_classes)
        title(['statistical differences for ' ...
            regexprep(feat_names2analyse{ff},'_',' ') ', sorted by ' ...
            variables_names{vv}])
    end
    % figure()
    % clustergram(pStatisticalDifferences.NN50)

end

%% Perform feature selection using auc for each feature individually
% RES_KFOLD = auc_feat_selection(feature_dataset_240min_before_seizure_3D_interpolated, ...
%     feat_names2analyse);


%% Perform feature selection using a linear classifier for features
% combinatin 3 by 3, using the same model for all training data

% feat_comb = 3;
% [results_final_mean_kfold, results_final_std_kfold] = ...
%     linear_classification_feat_selection(feature_dataset_240min_before_seizure_3D, ...
%     feat_names2analyse, feat_comb);



% % logistic regression
% results_SE_logistic_test = results_final_mean_kfold(2:end,5);
% results_SP_logistic_test = results_final_mean_kfold(2:end,6);
% results_AUC_logistic = results_final_mean_kfold(2:end,2);
% % geometric mean
% results_GM_logistic = sqrt(cell2mat(results_SE_logistic_test).*cell2mat(results_SP_logistic_test));
% sorted_features_logistic_3D = sort_feature_value(results_AUC_logistic,[feat_combos,indexes', ...
%     results_SE_logistic_test, results_SP_logistic_test, results_AUC_logistic, num2cell(results_GM_logistic)]);

% % lda
% results_AUC_lda = results_final_mean_kfold(2:end,7);
% results_SE_lda_test = results_final_mean_kfold(2:end,10);
% results_SP_lda_test = results_final_mean_kfold(2:end,11);

% % geometric mean
% geometric_mean = num2cell(sqrt(cell2mat(results_SE_lda_test).*cell2mat(results_SP_lda_test)));
% feat_combos = results_final_mean_kfold(2:end,1);
% indexes = num2cell(1:length(feat_combos));
% sorted_features_lda_3D = sort_feature_value(results_AUC_lda,[feat_combos, indexes', ...
%     results_SE_lda_test, results_SP_lda_test, geometric_mean]); %


% % linear svm
% results_SE_linear_svm_test = results_final_mean_each_seizure_3D(2:end,16);
% results_SP_linear_svm_test = results_final_mean_each_seizure_3D(2:end,17);
% % geometric mean
% geometric_mean = num2cell(sqrt(cell2mat(results_SE_linear_svm_test).*cell2mat(results_SP_linear_svm_test)));
% sorted_features_linear_svm_3D = sort_feature_value(num2cell(geometric_mean),[feat_combos,indexes', ...
%     results_SE_linear_svm_test, results_SP_linear_svm_test, results_AUC_linear_svm]);

%% Perform feature selection using a linear classifier for features
% combinatin 2 by 2 and 3 by 3, using a given model for a each seizure

classifiers_name = {'Logistic', 'LDA','Linear_SVM'}
feat_comb = 2;
plotFigure = 0;
samples_each_class = 20;

% linear classification in 2D
feature_relevance_2D = getLinearClassificationEachDimension(...
    feature_dataset_240min_before_seizure_3D_interpolated, feat_names2analyse, ...
    feat_comb, samples_each_class, classifiers_name, plotFigure);

% linear classification in 3D
feature_relevance_3D = getLinearClassificationEachDimension(...
    feature_dataset_240min_before_seizure, feat_names2analyse, feat_comb, ...
    samples_each_class, classifiers_name, plotFigure);
% the best results were found for the lda classifier
% the results for this classifier were


%%
saveSelectedFeatures = 1;

if saveSelectedFeatures==1
    clear feature_dataset_240min_before_seizure_clustering
    clear feature_dataset_240min_before_seizure_3D_clustering
    indexes_selected = [2,4,6:9,11,16,19,21:26,28,31:32];
    indexes_discarded = [1,3,5,10,12:15,17:18,20,27,29:30];
    fields = fieldnames(feature_dataset_240min_before_seizure_interpolated);
    feature_dataset_240min_before_seizure_clustering = rmfield(feature_dataset_240min_before_seizure_interpolated,feat_names2analyse(indexes_discarded)');
    feature_dataset_240min_before_seizure_3D_clustering = feature_dataset_240min_before_seizure_3D_interpolated;
    feature_dataset_240min_before_seizure_3D_clustering(:,indexes_discarded,:) = [];
    feat_names2analyse(indexes_discarded) = [];
    
    structure2clusterPreIctalStudy.feature_dataset_240min_before_seizure_3D = feature_dataset_240min_before_seizure_3D_clustering;
    structure2clusterPreIctalStudy.feat_names2analyse = feat_names2analyse;
    structure2clusterPreIctalStudy.feature_dataset_240min_before_seizure = feature_dataset_240min_before_seizure_clustering;
    save(fullfile(results_folder_path, 'structure2clusterPreIctalStudy.mat'), ...
        'structure2clusterPreIctalStudy')
else
    load(fullfile(results_folder_path, 'structure2clusterPreIctalStudy.mat'))
    feature_dataset_240min_before_seizure_3D_clustering = structure2clusterPreIctalStudy.feature_dataset_240min_before_seizure_3D_clustering;
    feature_dataset_240min_before_seizure_clustering = structure2clusterPreIctalStudy.feature_dataset_240min_before_seizure_clustering;
    feat_names2analyse = structure2clusterPreIctalStudy.feat_names2analyse;
    clear structure2clusterPreIctalStudy
end




% include_time_feat = 0;
%
% if include_time_feat
%     n_features = n_features+1;
%     feat_names2analyse = [feat_names2analyse; {'time'}];
%     time_dataset_3D = reshape(time_dataset,size(time_dataset,1), 1, n_wins);
%     feature_dataset_240min_before_seizure_3D_interpolated2 = cat(2, ...
%         feature_dataset_240min_before_seizure_3D_interpolated,time_dataset_3D);
% end


%% PERFORM CLUSTERING

functionsFolder = 'FunctionsClustering';
folder2savePath = fullfile(cd, functionsFolder);
if exist(folder2savePath, 'dir')
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder2savePath));
end

% tic
clustering_by_seizure(feat_names2analyse, feat_comb, pats_name, ...
    names_struct_separated, clustering_folder_path, outer_folder_path);
% time_clustering = toc;
% disp(['Elapsed time for clustering: ' num2str(time_clustering/60/60) ' hours'])
% save([feat_folder '_clustering.mat', time_clustering])

if exist(folder2savePath, 'dir')
    % Remove that folder plus all subfolders to the path.
    rmpath(genpath(folder2savePath))
end


% DONE ONLY WITH DBSCAN:
% clustering_by_seizure_over_time(feat_names2analyse, feat_comb, ...
%     names_struct, clustering_folder_path, outer_folder_path);


% to run the following code, it is require to use
% feature_dataset_240min_before_seizure_3D_clustering which was built by
% cropping the seizures clustering solutions:
% feature_dataset_240min_before_seizure_3D_clustering = feature_dataset_240min_before_seizure_3D;
% results_clustering_comb = clustering_by_feature_combination( ...
%     feature_dataset_240min_before_seizure_3D_clustering, ...
%     feat_names2analyse, feat_comb, clustering_folder_path, k_cluster2test_vec, ...
%     n_epsilon4dbscan, computeFromScratch, init_seiz);




%% PERFORM CLUSTERING EVALUATION *******************************************


% tic
% results = clustering_solution_evaluation_by_feature_combination(feat_names2analyse, ...
%     feat_comb, names_struct, clustering_folder_path, n_seizures, n_wins);
% time_cluster_evaluation = toc;
% disp(['Elapsed time for cluster solution evaluation: ' ...
%     num2str(time_cluster_evaluation/60/60) ' hours'])
% save([feat_folder '_cluster_evaluation.mat', time_cluster_evaluation])

functionsFolder = 'FunctionsClusterSolutionEvaluation';
folder2savePath = fullfile(cd, functionsFolder);
if exist(folder2savePath, 'dir')
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder2savePath));
end

clustering_solution_evaluation_by_seizure(feat_names2analyse, feat_comb, ...
    pats_name, clustering_folder_path)


if exist(folder2savePath, 'dir')
    % Remove that folder plus all subfolders to the path.
    rmpath(genpath(folder2savePath))
end

%% PERFORM PROTOTYPE CLUSTERING AMONG PATIENTS *****************************
% understand if all seizures clustering solutions comprised the same number
% of clusters
prototype_patient_clustering(feat_names2analyse, feat_comb, names_struct, ...
    results_folder_path, clustering_folder_path);

seizure_similarity_features(feat_names2analyse, ...
    feat_comb, names_struct, results_folder_path, clustering_folder_path, ...
    feature_dataset_240min_before_seizure_3D);


% first seizure of pat_402, kmeans_k2
ind1 = results_feat_eval.DI.kmeans_k2>=0.09;
ind2 = results_feat_eval.DI.kmeans_k2<0.2;
ind3 = find(all([ind1'; ind2']));
save('ind3.mat', 'ind3')


%% ANALYSE CLUSTERING RESULTS

functionsFolder = 'FunctionsResultsClusteringAnalysis';
folder2savePath = fullfile(cd, functionsFolder);
if exist(folder2savePath, 'dir')
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder2savePath));
end


prototype_patient_clustering_by_seizure(feat_names2analyse, ...
    feat_comb, pats_name, names_struct_separated, results_folder_path, ...
    clustering_folder_path, time_vec_original, percent_noise_original)

% run the following functions after running 
% prototype_patient_clustering_by_seizure.m
save_info_patients = getSaveInfoPatients(feat_names2analyse, feat_comb, ...
    pats_name, names_struct_separated, results_folder_path, ...
    patients_info_final, matrix_patient_info, n_wins, time_vec_original, ...
    percent_noise_original, clustering_folder_path);


getSmallerClusterStartingTime(feat_names2analyse, feat_comb, ...
    pats_name, names_struct_separated, results_folder_path, n_wins, ...
    clustering_folder_path, matrix_patient_info, time_vec_original, ...
    percent_noise_original);

% run the following function after running getSaveInfoPatients.m
metadataSeizure = plotMetadataQ1Q2(metaData, pats_name, results_folder_path);

[mat_all_seizures_all_methods, mat_all_seizures_all_methods_selected, ...
    mat_all_seizures_all_methods_selected_final, mat_all_seizures_all_methods_no_intervals] ...
    = analyseSmallerClusterStartingTime(feat_names2analyse, feat_comb, ...
    pats_name, names_struct, names_struct_separated, results_folder_path, ...
    matrix_patient_info, metadataSeizure, mat4AssociationLearning);


indexes_seiz_empty = cellfun(@isempty,mat_all_seizures_all_methods_selected_final(:,2));

seizures_preictal_info = mat_all_seizures_all_methods_selected_final(~indexes_seiz_empty, :);





observeSignal2Features(feat_names2analyse, feat_comb, feat_folder, ...
    windows_folder, outer_folder_path, patients_info_final)

getClusteringSolutions(feat_names2analyse, feat_comb, pats_name, ...
    names_struct_separated, results_folder_path, clustering_folder_path, ...
    matrix_patient_info, time_vec_original, percent_noise_original)

% get the duration of the cluster in min:
% (1) duration = end_time_cluster-start_time_cluster for continuous clusters
% (2) duration = (samples_cluster*5-5)/60
% get the samples of the cluster in min:
% (3) samples = (duration*60+5)/5 

if exist(folder2savePath, 'dir')
    % Remove that folder plus all subfolders to the path.
    rmpath(genpath(folder2savePath))
end


n_pat = numel(pat_names);


%% analyse clustering evaluation measures with the goal of defining DI threshold
compute_clust_indexes = 1;
clust_eval = {'DI', 'OD', 'ICV', 'C', 'CS', 'SI', ...
    'samples_smaller_cluster', 'samples_bigger_cluster', ...
    'NaN_samples', 'number_solutions', 'count_DI_solutions'};
folder2save = 'SelectionClusterSolutions';
folder2savePath = fullfile(cd, folder2save);
if ~exist(folder2savePath, 'dir')
    mkdir(folder2savePath)
end

save_number_solutions = cell(n_seiz+1,7);
save_number_solutions(2:end,1) = names_struct;
save_number_solutions(1,1:6) ={'seiz', 'total solutions', 'two clusters', ...
    'solutions no noise', 'min samples', 'solutions 2 clust non noise min samples'};

if compute_clust_indexes
    [save_DI, save_OD, save_ICV, save_C, save_CS, save_SI,...
        save_samples_smaller_cluster, save_samples_bigger_cluster, ...
        save_count_DI_solutions, save_number_solutions, save_NaN_samples] = ...
        analyse_DI_results(feat_names2analyse, feat_comb, pats_name, ...
        clustering_folder_path, save_number_solutions);
    for ll = 1:numel(clust_eval)
        eval(['save(fullfile(cd, folder2save, ''save_' clust_eval{ll} ...
            '''), ''save_' clust_eval{ll} ''')'])
    end
else
    for ll = 1:numel(clust_eval)
        eval(['load(fullfile(cd, folder2save, ''save_' clust_eval{ll} ...
            '''), ''save_' clust_eval{ll} ''')'])
    end
end

% *************************************************************************
disp('Total number of solutions:')
n_total_solutions = 4960*n_seiz*7
disp('Total number of solutions with 2 clusters:')
n_2cluster_solutions = sum(cell2mat(save_number_solutions(2:end,3)))
disp(['Total number of solutions with 2 clusters, no noisy samples, ' ...
    'minimum of 20 samples and DI>=15: ' num2str(save_count_DI_solutions)])
disp('Percentage of solutions obtained from Algorithm 1: ')
algorithm1_percent = save_count_DI_solutions/save_count_solutions*100;

% *************************************************************************

data2cluster = [save_DI save_OD save_ICV save_C save_CS save_SI ...
    save_samples_smaller_cluster save_samples_bigger_cluster];
[B,I] = sort(data2cluster(:,1));
data2cluster = data2cluster(I,:);
data2cluster_norm = (data2cluster-min(data2cluster))./(max(data2cluster)-min(data2cluster));


figure()
plot(1:size(data2cluster,1), data2cluster_norm)






dim = 3;
feat_combs = nchoosek(1:numel(clust_eval(1:end-2)),dim);
n_comb = size(feat_combs,1);

figure(10)
ind_DI = save_DI>=0.15;
for ii = 1:n_comb

    if dim==2
        eval(['scatter(save_' clust_eval{feat_combs(ii,1)} ...
            ', save_' clust_eval{feat_combs(ii,2)} ', [], ind_DI)'])
    else
        eval(['scatter3(save_' clust_eval{feat_combs(ii,1)} ...
            ', save_' clust_eval{feat_combs(ii,2)} ...
            ', save_' clust_eval{feat_combs(ii,3)} ', [], ind_DI)'])
        zlabel(clust_eval{feat_combs(ii,3)})
    end
    
    % yline(0.15)
    xlabel(clust_eval{feat_combs(ii,1)})
    ylabel(clust_eval{feat_combs(ii,2)}) 

    axis tight
    pause
end


ind_cluster_OD_SI = save_SI<=0.92 & save_OD<=30;

figure()
scatter(save_OD, save_SI, '*'), hold on
scatter(save_OD(ind_cluster_OD_SI), save_SI(ind_cluster_OD_SI), 'ro'), hold off

figure()
scatter3(save_OD, save_SI, save_DI, '*'), hold on
scatter3(save_OD(ind_cluster_OD_SI), save_SI(ind_cluster_OD_SI), ...
    save_DI(ind_cluster_OD_SI),'ro'), hold off


n_cluster_OD_SI = sum(ind_cluster_OD_SI)/numel(save_SI)*100

DI_in_cluster_OD_SI = save_DI(ind_cluster_OD_SI);


variables2inspect = {'OD', 'CS', 'SI'};

% apply a clustering method to CS vs SI vs OD




data2cluster = [save_OD save_CS save_SI];
indexes_less_samples = save_OD<=100;
data2cluster = data2cluster(indexes_less_samples,:);
data2cluster_norm = (data2cluster-min(data2cluster))./(max(data2cluster)-min(data2cluster));

dimension = size(data2cluster,2);

opts_KMEANS = statset('Display','final');
[clusteringSolution, centroids] = kmeans(data2cluster_norm, 2, ...
    'Distance', 'sqeuclidean', 'Options', opts_KMEANS, 'MaxIter', 1000);


Z = linkage(data2cluster_norm,'ward','spearman');
clusteringSolution = cluster(Z,'maxclust',2);




plotFigure = 1
[clusteringSolution, centroids] = ...
    applyGaussianMixtureModelsClustering(data2cluster_norm, plotFigure, ...
    variables2inspect);


clusteringSolution = spectralcluster(data2cluster_norm, 2, 'Distance', ...
    'cityblock');


[centroids, clusteringSolutionDBSCANFinal] = getDBSCANClustering(dimension, ...
    data2cluster_norm, variables2inspect, [], 1, 4);
                
                
                
figure(20)
plotClusterSolution(data2cluster_norm, clusteringSolution, ...
    centroids, [], unique(clusteringSolution), 2, ...
    dimension, 1, variables2inspect, [], [], []);    

yline(0.93)

smaller_cluster_percentage = count_DI_solutions/count_solutions*100;

figure()
histogram(save_DI)

%%
% find accepted solutions with ApEn and SampEn
close all
clc
feat_combs = nchoosek(1:numel(feat_names2analyse),feat_comb);
indexes_ApEn_SampEn = find(feat_combs(:,2)==21 & feat_combs(:,3)==21); % &feat_combs(:,3)==21

indexes_empty = ~cellfun(@isempty,mat_all_seizures_all_methods_selected(:,8));
mat_all_seizures_all_methods_selected2 = mat_all_seizures_all_methods_selected(indexes_empty,:);

for ss = 1:length(mat_all_seizures_all_methods_selected2)
    get_feat = mat_all_seizures_all_methods_selected2{ss,8}{3};
    
    [LIA,LOCB] = ismember(indexes_ApEn_SampEn,get_feat);
    
    if any(LIA)
        disp('Found')
        indexes_ApEn_SampEn(LIA)
        [mat_all_seizures_all_methods_selected2{ss,1} ...
            mat_all_seizures_all_methods_selected2{ss,8}]
        mat_all_seizures_all_methods_selected2{ss,8}{3}(LOCB(LIA))
        feat_combs(indexes_ApEn_SampEn(LIA),:)
    end
    
end

%%
analyseCentroidTrajectory(feat_names2analyse, feat_comb, names_struct, ...
    results_folder_path, clustering_folder_path, T)

%%

%% how many save_clustering_data.mat files are missing

check = checkExistingFiles(pats_name, patients_info_final);

%% Plot the percent_noise_detection for patient 402 with and whitout
% annotated noisy intervals and R peaks

clc
clear all
close all

% get noise percentage for overlaped non annotated patient 402
results_folder_path = 'ResultsOverlap'
load(fullfile(cd, results_folder_path, 'feature_dataset_240min_before_seizure_all_feat.mat'))
percent_noise = feature_dataset_240min_before_seizure_all_feat.percent_detected_noise;
[n_seiz, n_wins1] = size(percent_noise)


% get noise percentage for overlaped annotated patient 402
results_folder_path = 'ResultsOverlapAnnotated'
load(fullfile(cd, results_folder_path, 'feature_dataset_240min_before_seizure_all_feat.mat'))
percent_noise_annotated = feature_dataset_240min_before_seizure_all_feat.percent_detected_noise;
[n_seiz, n_wins2] = size(percent_noise_annotated)


figure()
subplot(511)
plot(1:n_wins1,percent_noise(1,:)), hold on
plot(1:n_wins2,percent_noise_annotated(1,:),'--'),
axis tight
ylim([0 100])
subplot(512)
plot(1:n_wins1,percent_noise(2,:)), hold on
plot(1:n_wins2,percent_noise_annotated(2,:),'--'),
axis tight
ylim([0 100])
subplot(513)
plot(1:n_wins1,percent_noise(3,:)), hold on
plot(1:n_wins2,percent_noise_annotated(3,:),'--'),
axis tight
ylim([0 100])
subplot(514)
plot(1:n_wins1,percent_noise(4,:)), hold on
plot(1:n_wins2,percent_noise_annotated(4,:),'--'),
axis tight
ylim([0 100])
subplot(515)
plot(1:n_wins1,percent_noise(5,:)), hold on
plot(1:n_wins2,percent_noise_annotated(5,:),'--'),
axis tight
ylim([0 100])

