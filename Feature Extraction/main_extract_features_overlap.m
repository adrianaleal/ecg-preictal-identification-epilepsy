fclose all; clear; close all; clc;
set(groot, 'defaultAxesFontSize',11,'defaultTextFontSize',11)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex')


% p = gcp;
% delete(p)


code_language = 'matlab'

if strcmp(code_language,'octave')
    outer_folder_path = pwd(); %cd();
else
    outer_folder_path = cd();
end
clear code_language

adrianaPC = 1;

if adrianaPC==1
    addpath(genpath(fullfile(outer_folder_path, 'utils')))
    addpath(genpath(fullfile(outer_folder_path, 'epi2pub_with_changes', 'Toolbox','CRPtool')))
end
clear adrianaPC
addpath(genpath(fullfile(outer_folder_path, 'HRVFeatures')))

%% TO DEFINE
file_type = 'resultsNoiseDet'% 'resultsNoiseDet' % 'hrv_features'%
plotFigures = 0;
computeFeatures = 1;
getOverlappedWindows = 1; % 1 --> compute, 0 --> load already stored data
correctRR = 1; % 1 --> compute, 0 --> load already stored data
plotRRCorrectionFigure = 0;


noise_folder = 'NoiseDetectionAnnotated'; % Annotated
feat_folder = 'ExtractedFeaturesNoEctopicBeatsOverlapAnnotated'; %Annotated
overlapped_windows_folder = 'OverlappedWindowsAnnotated'; %Annotated
RR_intervals_folder = 'RRIntervalsCorrectedAnnotated'; %Annotated


%%

figureFolder2save = fullfile(cd, 'PrepareFeatureDataset2Clustering', ...
    'FiguresAnalysedECGSignals');

if ~exist(figureFolder2save, 'dir')
    mkdir(figureFolder2save)
end
            
            
results_folder_path = fullfile(outer_folder_path, 'ResultsEPILEPSIAE');
patients = dir(fullfile(results_folder_path, noise_folder ,'pat_*'));
pat_names = {patients.name}';
split1 = regexp(pat_names, '_', 'split');
get_pat_ids = vertcat(split1{:});
[~,I] = sort(cellfun(@str2double,get_pat_ids(:,2)));
pat_names = pat_names(I);
clear I get_pat_ids split1 patients outer_folder_path

% patients_info_final = [];
load('patients_info_final_preictal_study_240min_before_seiz')
include_subdir = 1;



% checkFeatures = checkFeatureExtractionFiles(pat_names, patients_info_final, ...
%     results_folder_path);

feature_dataset_240min_before_seizure = struct;

% parpool('local',feature('numcores'),'IdleTimeout', Inf)
parfor pp = 1:6%:length(pat_names)
    
    tic
    pat_name = pat_names{pp}
    % adm_folders = dir(fullfile(dir_name,file_name,'adm_*'));
    
    patient_info = patients_info_final(pp)
    
    h = struct;
    
    count_files = 0;
    count_seizures = 0;
    h.removed_samples_seiz_file = [];
    h.removed_samples_first_file = [];
    save_last_recording_windows = [];
    [h, count_seizures, count_files, save_last_recording_windows] = ...
        GetAllFilesInDir4FeatureExtraction(fullfile(results_folder_path, ...
        noise_folder, pat_name), h, include_subdir, count_files, ...
        count_seizures, patient_info, file_type, save_last_recording_windows);
    
    h = rmfield(h,'noise_dataset');
    h = rmfield(h,'get_n_5min_until_seizure_file');
    
    
    % to compute the time gap vector **************************************
    removed_samples_seiz_file = h.removed_samples_seiz_file;
    removed_samples_first_file = h.removed_samples_first_file;
    h = rmfield(h,'removed_samples_seiz_file');
    h = rmfield(h,'removed_samples_first_file');
    ind_files_with_seizure = find(~cellfun(@isempty, ...
        {patients_info_final(pp).files.non_discarded_seizures_info}));
    % get the number of files to analyse before seizure, required to obtain
    % 240 min of data:
    n_files_before_seizure = vertcat(patients_info_final(pp).files(ind_files_with_seizure).n_files_before_seizure);
    
    ind_files2analyse = [ind_files_with_seizure'-n_files_before_seizure, ...
        ind_files_with_seizure'];
    n_samples_in_files = vertcat(patients_info_final(pp).files.n_samples);
    files_gap_secs = [0; vertcat(patients_info_final(pp).files.info_gap_secs)];
    
    
    % ver = [vertcat(patients_info_final(pp).files.n_seiz), n_samples_in_files, ...
    %     files_gap_secs];
    % *********************************************************************
    
    
    
    h_fields = fieldnames(h);
    h.names = strcat([pat_name '_seizure'],cellstr(num2str((1:count_seizures)')));
    plotFigure = zeros(count_seizures,1);
    % plotFigure(2) = 1;
    
    hrv_comp_time_pat = cell(count_seizures,2);
    hrv_all_feat_comp_time_pat = hrv_comp_time_pat;
    
    for ss = 1:count_seizures
        
        pat_seiz_name = [pat_name '_seizure' num2str(ss)];
        
        seiz_noise_det = h.(h_fields{ss});
        time_entire_seiz = horzcat(seiz_noise_det.time_entire_signal);
        
        inds_negative = find(time_entire_seiz(2:end)-time_entire_seiz(1:end-1)<0)';
        check_negative_time = [inds_negative, ...
            (time_entire_seiz(inds_negative+1)-time_entire_seiz(inds_negative))'];
        
        fs = patient_info.samp_rate;
        
        %% COMPUTE THE TIME GAP LOGICAL VECTOR ****************************
        [gap_time, gap_logical] = getFileGap(ind_files2analyse(ss,:), ...
            files_gap_secs, n_samples_in_files, seiz_noise_det, fs, ...
            removed_samples_first_file(ss), removed_samples_seiz_file(ss));
        % IMPORTANT: THE GAP VECTOR IS NOT ENTIRELY EQUAL TO THE DATA
        % VECTOR --> makes sense as the gap samples are not part of the
        % time_entire_seiz vector
        % length_time_gap_logical_final = length(gap_logical_final)
        % length_time_entire_seiz = length(time_entire_seiz)
        
        
        %% COMPUTE THE OVERLAPPED 5-MIN WINDOWS ***************************

        overlapped_windows_folder_path = fullfile(results_folder_path, ...
            overlapped_windows_folder, pat_name, ['seizure' num2str(ss)]);
        if ~exist(overlapped_windows_folder_path, 'dir')
            mkdir(overlapped_windows_folder_path)
        end
        
            
        if getOverlappedWindows==1
            [resultsOutput, save2FigurePlotECG] = getSignalWithOverlapSlidingWindow( ...
                seiz_noise_det, time_entire_seiz, gap_logical, gap_time, ...
                fs, 0, pat_seiz_name);
            % plotFigure(ss)
            resultsOutput = rmfield(resultsOutput,'sig_preprocess');
            parsave(overlapped_windows_folder_path, {resultsOutput}, ...
                'resultsOutput', pat_seiz_name)
            parsave(figureFolder2save, {save2FigurePlotECG}, ...
                'save2FigurePlotECG', [pat_seiz_name, '_save2FigurePlotECG'])
        else
            resultsOutput = parload(overlapped_windows_folder_path, ...
                pat_seiz_name, 'resultsOutput');
        end

        % check for repeated R peaks 
        % for qq = 1:length(resultsOutput)
        %    time_RR_intervals_before = resultsOutput(qq).QRSindxs_entire_sig_nan;
        %    if ~isequal(unique(time_RR_intervals_before), time_RR_intervals_before)
        %        qq
        %        disp('asdasd')
        %    end
        % end
        if plotRRCorrectionFigure
            figure()
            for ll = 1:size(save2FigurePlotECG.RRI_240min)
                plot(save2FigurePlotECG.time_Rpeaks_240min{ll}, ...
                    save2FigurePlotECG.RRI_240min{ll}, 'b')
                hold on
                plot(save2FigurePlotECG.time_Rpeaks_240min_corrected{ll}, ...
                    save2FigurePlotECG.RRI_240min_corrected{ll}, '--r')
            end
            hold off
            axis tight
            legend('RR intervals', 'RR intervals corrected')
            ylabel('ms')
            xlabel('Time (hour)')
            
            savefig(fullfile(figureFolder2save, [pat_seiz_name '_RRI_signal']))
            close all
            save(fullfile(cd, 'PrepareFeatureDataset2Clustering', ...
                'FiguresAnalysedECGSignals', 'save2FigurePlotECG.mat'), ...
                'save2FigurePlotECG')
        end
        

        %% get the RR intervals in each 5min window
        disp('>> RR intervals analysis <<')
        close all
        
        RR_intervals_folder_path = fullfile(results_folder_path, ...
            RR_intervals_folder, pat_name, ['seizure' num2str(ss)]);
        if ~exist(RR_intervals_folder_path, 'dir')
            mkdir(RR_intervals_folder_path)
        end
        
        if correctRR==1
            plotFigureRR = 1;
            sig_preprocess_entire_seiz = horzcat(seiz_noise_det.sig_preprocess);
            resultsRR = getRRInterval(resultsOutput, time_entire_seiz, ...
                sig_preprocess_entire_seiz, fs, plotFigureRR);
            close all
           
            parsave(RR_intervals_folder_path, {resultsRR}, 'resultsRR', ...
                pat_seiz_name)
        else
            %load resultsRR
            resultsRR = parload(RR_intervals_folder_path, pat_seiz_name, 'resultsRR');
        end
        
        
        
        %% UNIVARIATE ECG FEATURE EXTRACTION
        
        % OR **************************************************************
        %         QRS_ComplexDurSec = [];
        %         [v_RIndices, v_DataInFilt] = RWaveDet(ecg_sig, samp_freq,...
        %             QRS_ComplexDurSec, 1);
        % *****************************************************************
        
        % when I want to observe the features' figures for a given window:
        % example: 5-min window after 80 min before seizure and 5-min
        % window after 25 min before seizure onset

        if computeFeatures
           
            disp('>> UNIVARIATE ECG FEATURE EXTRACTION <<')
            close all
            
            % in order to know when the windows that are captured to the
            % paper occur before the seizure onset:
            % time_RR_secs = [0; cumsum(diff(vertcat(resultsRR.time)))];
            % time_RR = time_RR_secs/60; % in min
            % time_win_1800 = time_RR(end)-time_RR(1800)
            % time_win_2800 = time_RR(end)-time_RR(2800)
            
            plotFigureFeat = 0;
            [hrv_features, hrv_comp_time] = getHRVFeatures(resultsRR, ...
                resultsOutput, plotFigureFeat);
            
            
%             [hrv_features, results] = HRVFeatureExtraction(resultsOutput, ...
%                 sig_preprocess_entire_seiz, time_entire_seiz, fs, ...
%                 plotFigureFeat);

            %% Plot features
            
            if plotFigures==1
                ecg_sig_preprocess = vertcat(resultsOutput.sig_preprocess);
                saveFig = 0;
                dispFigures(time_entire_seiz, ecg_sig_preprocess, ...
                    [], results, fs, recording_name, hrv_features, saveFig)
                %         close all
            end
            
            
            %% Save feature extraction results
            extracted_features_folder_path = fullfile(results_folder_path, ...
                feat_folder, pat_name, ['seizure' num2str(ss)]);
            plotsRQATT_folder_path = fullfile(results_folder_path, ...
                'PlotsRQATT', pat_name, ['seizure' num2str(ss)]);
            
            
            % SAVE PLOTS
            %if ~exist(plotsRQATT_folder_path, 'dir')
            %    mkdir(plotsRQATT_folder_path)
            %end
            
            %count_fig = 1;
            %FigList = flipud(findobj(allchild(0), 'flat', 'Type', 'figure'));
            %for iFig = 1:2:length(FigList)
            %    FigHandle1 = FigList(iFig);
            %    FigHandle2 = FigList(iFig+1);
            %    set(0, 'CurrentFigure', FigHandle1);
            %    set(0, 'CurrentFigure', FigHandle2);
            %    FigName1 = ['distance_matrix_' num2str(count_fig)];
            %    FigName2 = ['recurrence_plot' num2str(count_fig)];
            %    savefig(fullfile(plotsRQATT_folder_path, [FigName1 '.fig']));
            %    savefig(fullfile(plotsRQATT_folder_path, [FigName2 '.fig']));
            %    export_fig(FigHandle1,fullfile(plotsRQATT_folder_path, [FigName1 '.png']),'-transparent')
            %    export_fig(FigHandle2,fullfile(plotsRQATT_folder_path, [FigName2 '.png']),'-transparent')
            %    count_fig = count_fig +1;
            %end
            %close all
            
            
            % SAVE FILES
            % see if a folder with the previous name exists and if not build
            % that folder in the specified path
            if ~exist(extracted_features_folder_path, 'dir')
                mkdir(extracted_features_folder_path)
            end
            
            stop_time = toc;
            
            parsave(extracted_features_folder_path, {hrv_features}, ...
                'hrv_features', pat_seiz_name)
            
            parsave(extracted_features_folder_path, {stop_time}, ...
                'stop_time', [pat_seiz_name '_comp_time_all_feat'])
            
            parsave(extracted_features_folder_path, {hrv_comp_time}, ...
                'hrv_comp_time', [pat_seiz_name '_comp_time'])

            pat_name
            fprintf('[%.3f sec] >> Elapsed total time for patient \n', stop_time);
            
            
            
        end
  
    end
    
end




