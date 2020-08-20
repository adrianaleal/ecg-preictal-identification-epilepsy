fclose all; clear; close all; clc;
set(groot, 'DefaultAxesFontSize',11,'defaultTextFontSize',11)
set(groot, 'DefaultAxesTickLabelInterpreter','tex');
set(groot, 'DefaultLegendInterpreter','tex');
set(groot, 'DefaultTextInterpreter','tex')


cd ..
outer_folder_path = cd;
% load variable patients_info_final:
load('patients_info_preictal_study_240min_before_seiz')
cd('Clustering')

% define the path to the folder containing the features extracted for each
% patient and seizure:
feature_folder_path = fullfile(outer_folder_path, 'ResultsFeatureExtraction');


dataset_folder_path = fullfile(outer_folder_path, 'FeatureDataset');
if ~exist(dataset_folder_path, 'dir')
    mkdir(dataset_folder_path)
end

% load feature dataset
load(fullfile(dataset_folder_path, 'structureData.mat'))
variables2load = {'feat_names2analyse', 'seizure_struct'}';
    
for ii = 1:numel(variables2load)
    eval([variables2load{ii} ' = structureData.' variables2load{ii} ';'])
end
clear structureData

seizure_names = {seizure_struct(:).seizure_name}';
index_nonempty = find(~cellfun(@isempty,seizure_names));
seizure_names = seizure_names(index_nonempty);

% get the patients' names to initialize the structure to save the results
C = cellfun(@(x)strsplit(x, '_' ), seizure_names, 'UniformOutput', false);
seizure_names_separated = vertcat(C{:});
patients_name = unique(seizure_names_separated(:,2),'stable');
n_pat = numel(patients_name);



%% Perform cluster solution evaluation


% defined number of features to combine(e.g. 3-by-3 or 2-by-2)
feat_comb = 3;

% define the folder path to save the results
clustering_folder_path = fullfile(cd, 'ResultsClustering', ...
    ['ResultsClusteringFeatComb' num2str(feat_comb)]);


% add the clustering evaluation functions to the path
functionsFolder = 'FunctionsClusterSolutionEvaluation';
folder2savePath = fullfile(cd, functionsFolder);
if exist(folder2savePath, 'dir')
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder2savePath));
end


functionsFolder = 'utils';
folder2savePath = fullfile(outer_folder_path, functionsFolder);
if exist(folder2savePath, 'dir')
    % Add that folder plus all subfolders to the path.
    addpath(genpath(folder2savePath));
end



% tic
clustering_solution_evaluation_by_seizure(feat_names2analyse, feat_comb, ...
    patients_name, clustering_folder_path)


if exist(folder2savePath, 'dir')
    % Remove that folder plus all subfolders to the path.
    rmpath(genpath(folder2savePath))
end

