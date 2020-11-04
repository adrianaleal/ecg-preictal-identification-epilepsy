function [] = loadFeatureDataset(results_folder_path, features_folder_path, ...
    variables2load, patients_info)


patients_folder_path = fullfile(features_folder_path, 'pat_*');
pats_name = getPatientsNames(patients_folder_path);

% get the number of seizures for each patient:
n_seizures_each_pat = vertcat(patients_info.n_seizures)- ...
    vertcat(patients_info.n_discarded_seizures_by_240min)- ...
    vertcat(patients_info.n_discarded_seizures_by_no_file);
n_seizures = sum(n_seizures_each_pat)

% initialize a structure to save information for each seizure:
c = cell(n_seizures,1);
[seizure_struct(1:n_seizures).seizure_name] = deal(c{:});
[seizure_struct(1:n_seizures).time_vec] = deal(c{:});
[seizure_struct(1:n_seizures).percent_noise] = deal(c{:});
[seizure_struct(1:n_seizures).n_5min_wins] = deal(c{:});

comp_time_total = [];
count_seizures = 0;
start_seiz_index = 1;

for pp = 1:length(pats_name)
    file_name = pats_name{pp}
    
    seiz_folders = dir(fullfile(features_folder_path, file_name));
    seiz_folders = seiz_folders(~ismember({seiz_folders.name},{'.','..'}));
    n_seizures_pat = numel(seiz_folders);
    
    % get seizures indexes in the total number of seizures;
    indexes_seizures = start_seiz_index:start_seiz_index+n_seizures_pat-1;
    start_seiz_index = start_seiz_index+n_seizures_pat;
    seizure_name = strcat([file_name '_' num2str(pp) ...
        '_seizure_'], cellstr(num2str((1:n_seizures_pat)')));
    [seizure_struct(indexes_seizures).seizure_name] = deal(seizure_name{:});
    
    for ss = 1:numel(seiz_folders)
        count_seizures = count_seizures+1;
        
        % get information regarding the computational time required to
        % obtain each feature:
        comp_time = load(fullfile(seiz_folders(ss).folder, ...
            seiz_folders(ss).name, [file_name '_seizure' num2str(ss) '_comp_time']));
        comp = fieldnames(comp_time);
        comp_time_total = horzcat(comp_time_total, comp_time.(comp{1}));
        
        
        % load HRV features:
        load(fullfile(seiz_folders(ss).folder,seiz_folders(ss).name, ...
            [file_name '_seizure' num2str(ss)]))
        feat_names = fieldnames(hrv_features);
        
        seizure_struct(indexes_seizures(ss)).n_5min_wins = numel(hrv_features);
        seizure_struct(indexes_seizures(ss)).time_vec = {vertcat(hrv_features.time)};
        seizure_struct(indexes_seizures(ss)).percent_noise = {vertcat(hrv_features.percent_detected_noise)};
        
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
end

feature_dataset_240min_before_seizure.names = seizure_struct.seizure_name;


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
save(fullfile(results_folder_path, 'feature_dataset_240min_before_seizure.mat'), ...
    'feature_dataset_240min_before_seizure')

% get feature extraction computational time
all_values = squeeze(cell2mat(struct2cell(comp_time_total)))';
feats = [''; fieldnames(comp_time_total)];
ct_mean = nanmean(all_values)';
ct_std = nanstd(all_values)';
compTimeTable = table(feats, ct_mean, ct_std);

for ii = 1:numel(variables2load)
    eval(['structureData.' variables2load{ii} ' = ' variables2load{ii} ';'])
end

save(fullfile(results_folder_path, 'structureData.mat'), 'structureData')

end