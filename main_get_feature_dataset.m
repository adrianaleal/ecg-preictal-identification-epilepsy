fclose all; clear; close all; clc;
set(groot, 'DefaultAxesFontSize',11,'defaultTextFontSize',11)
set(groot, 'DefaultAxesTickLabelInterpreter','tex');
set(groot, 'DefaultLegendInterpreter','tex');
set(groot, 'DefaultTextInterpreter','tex')


cd ..
outer_folder_path = cd;
% load variable patients_info_final:
load('patients_info_final_preictal_study_240min_before_seiz')
cd('Feature Redundancy Study')

% open folder containing the feature dataset:
dataset_folder_path = fullfile(cd,'FeatureDataset');
if ~exist(dataset_folder_path, 'dir')
    mkdir(dataset_folder_path)
end

% get the variable to 
computeFromScratch = 0;
variables2load = {'feature_dataset_240min_before_seizure_3D', ...
        'feature_dataset_240min_before_seizure', 'feat_names2analyse', ...
        'n_wins', 'n_seizures', 'compTimeTable', 'seizure_struct'}';
    
if computeFromScratch==1
    loadFeatureDataset(dataset_folder_path, outer_folder_path, ...
        variables2load, patients_info_final)
end
load(fullfile(dataset_folder_path, 'structureData.mat'))
for ii = 1:numel(variables2load)
    eval([variables2load{ii} ' = structureData.' variables2load{ii} ';'])
end
clear structureData

load(fullfile(dataset_folder_path, 'feature_dataset_240min_before_seizure_3D.mat'), ...
    'feature_dataset_240min_before_seizure_3D')
load(fullfile(dataset_folder_path, 'feature_dataset_240min_before_seizure_all_feat.mat'), ...
    'feature_dataset_240min_before_seizure_all_feat')
load(fullfile(dataset_folder_path, 'feature_dataset_240min_before_seizure.mat'), ...
    'feature_dataset_240min_before_seizure')

seizure_names = {seizure_struct(:).seizure_name};

% get the patients' names to initialize the structure to save the results
C = cellfun(@(x)strsplit(x, '_' ), seizure_names, 'UniformOutput', false);
seizure_names_separated = vertcat(C{:});
pats_name = unique(seizure_names_separated(:,2),'stable');
n_pat = numel(pats_name);
clear C pats_name 



