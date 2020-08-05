function [] = loadFeatureDataset(results_folder_path, outer_folder_path, ...
    variables2load)


feat_folder = 'ExtractedFeaturesNoEctopicBeatsOverlap';
if ~isempty(strfind(results_folder_path,'Annotated'))
    feat_folder = [feat_folder 'Annotated'];
end

patients_folder_path = fullfile(outer_folder_path, ...
    'ResultsEPILEPSIAE', feat_folder ,'pat_*');
pats_name = getPatientsNames(patients_folder_path);

mat_number_seizures = [pats_name, cell(numel(pats_name),1)];

names_struct = [];
time_vec_original = [];
percent_noise_original = [];
comp_time_total = [];
n_5min_wins_per_seizure = [];
count_seizures = 0;

for pp = 1:length(pats_name)
    file_name = pats_name{pp}
    
    seiz_folders = dir(fullfile(outer_folder_path,'ResultsEPILEPSIAE', ...
        feat_folder, file_name));
    seiz_folders = seiz_folders(~ismember({seiz_folders.name},{'.','..'}));
    num_seizures_pat = numel(seiz_folders);
    mat_number_seizures(pp,2) = num2cell(num_seizures_pat);
    
    for ss = 1:numel(seiz_folders)
        count_seizures = count_seizures+1;
        
        % get information regarding the computational time required for to
        % obtained each feature:
        comp_time = load(fullfile(seiz_folders(ss).folder, ...
            seiz_folders(ss).name, [file_name '_seizure' num2str(ss) '_comp_time']));
        comp = fieldnames(comp_time);
        comp_time_total = horzcat(comp_time_total, comp_time.(comp{1}));
        
        
        % load HRV features:
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

feature_dataset_240min_before_seizure.names = names_struct;

n_seizures = vertcat(patients_info_final.n_seizures)- ...
    vertcat(patients_info_final.n_discarded_seizures_by_240min)- ...
    vertcat(patients_info_final.n_discarded_seizures_by_no_file);
check_n_seizures_patients_info_final = sum(n_seizures)
check_n_seizures_each_seizure = sum(cell2mat(mat_number_seizures(:,2)))
clear mat_number_seizures

% put the invalid windows to NaN ******************************************
% there are already some samples in variable segment_size assigned to NaN

feature_dataset_240min_before_seizure = editInvalidWindows(feat_names, ...
    feature_dataset_240min_before_seizure);

% get the 3D matrix to perform feature redundancy analysis ****************
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

% change matrix dimensions to [n_seizures, n_feat, n_wins]
feature_dataset_240min_before_seizure_3D = permute(feature_dataset_240min_before_seizure_3D,[1,3,2]);

% save 3D matrix with the feature information
save(fullfile(results_folder_path, 'feature_dataset_240min_before_seizure_3D.mat'),...
    'feature_dataset_240min_before_seizure_3D')

% save original structure with all the information from Feature Extraction:
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

for ii = 1:numel(variables2load)
    eval(['structure2selectFeatPreIctalStudy.' variables2load{ii} ' = ' ...
        variables2load{ii} ';'])
end
    
save(fullfile(results_folder_path, 'structure2selectFeatPreIctalStudy.mat'), ...
    'structure2selectFeatPreIctalStudy')

end