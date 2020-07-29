function [h, count_seizures, countFiles2analyse, ...
    save_last_recording_windows] = GetAllFilesInDir4FeatureExtraction(dir_path, ...
    h, include_sub, countFiles2analyse, count_seizures, patient_info, ...
    file_type, save_last_recording_windows)

all_files = dir(dir_path); % get all files in directory
%gap_secs=0;
all_files = all_files(~ismember({all_files.name},{'.','..'}));

analyseFileAgain = 0;

%% when a given recording has to be analysed for more than one seizure:

recordingsFoundInFolder = cellstr(vertcat(all_files.name))
[~,recordingsFoundInFolder,file_type2] = cellfun(@fileparts,recordingsFoundInFolder,'un',0);

if any(strcmp(file_type2,'.mat'))
    % check for repeated files in recordings2analyse:
    recordings2analyse = patient_info.files2analyse(:,1)
    [~,recordings2analyse,~] = cellfun(@fileparts,recordings2analyse,'un',0);
    a = unique(recordings2analyse);
    b = cell2mat(cellfun(@(x) sum(ismember(sort(recordings2analyse),x)),a,'un',0));

    if any(b>1) % there are repeated recordings
        names_repeated_files = a(b>1);
        for ll = 1:numel(names_repeated_files)
            name_repeated = names_repeated_files(ll);
            recordingsFoundInFolder = [recordingsFoundInFolder; name_repeated];
        end
        recordingsFoundInFolder = sort(recordingsFoundInFolder);
    end
    recordingsFoundInFolder = strcat(recordingsFoundInFolder,'.mat');
end

%%

for ii = 1:length(recordingsFoundInFolder)
    %% entire directory
    path = fullfile(dir_path, recordingsFoundInFolder{ii});
    if include_sub && exist(path, 'dir') == 7 % Check for subdirectories
        [h, count_seizures, countFiles2analyse, ...
            save_last_recording_windows] = GetAllFilesInDir4FeatureExtraction(path, ...
            h, include_sub, countFiles2analyse, count_seizures, patient_info, ...
            file_type, save_last_recording_windows);
    end
    
    
    if exist(path, 'file') ~= 2 || strcmp(recordingsFoundInFolder(ii), '.DS_Store')% Check for files only
        continue;
    end
    
    fprintf('\n')
    fprintf('\n')
    disp('****************************************************************')
    disp('                  NEW FILE ANALYZED                            ')
    disp('****************************************************************')
    path
    
    
    
    %% Load noise detection data
    file_name = recordingsFoundInFolder{ii}
    
    % check if the file belongs to the files of interest (e.g. 240 min
    % before the seizure)
    indFile2analyse = find(~cellfun(@isempty,strfind(patient_info.files2analyse(:,1),file_name(1:end-4))));
    % indFile2analyse may contain more than one index because the file may
    % be needed for the anlysis of two seizures
    
    ind_col_contains_seizure = find(strcmp(patient_info.files2analyse(1,:), 'contains_seizure'));
    indSeizure2 = strcmp(patient_info.files2analyse(indFile2analyse,ind_col_contains_seizure),'seiz')
    indFile2analyseSeiz = find(~cellfun(@isempty,strfind(patient_info.files2analyse(:,ind_col_contains_seizure),'seiz')));
    indSeizure = find(ismember(indFile2analyseSeiz,indFile2analyse));
    if numel(indFile2analyse)>1  && analyseFileAgain==0
        indFile2analyse = indFile2analyseSeiz(indSeizure);
        analyseFileAgain = 1;
    elseif numel(indFile2analyse)>1  && analyseFileAgain==1
        indFile2analyse = indFile2analyse(indFile2analyse~=indFile2analyseSeiz(indSeizure));
        analyseFileAgain = 1;
        indSeizure = [];
    end
    
    if ~isempty(indFile2analyse)
        S = load(path);
        var_names = fieldnames(S.(file_type));
        
        countFiles2analyse = countFiles2analyse+1;
        
        
        mat_data = S.(file_type);
        
        %%  COMPUTE THE TIME VECTOR FOR SIGNAL REPRESENTATION
        
        % compute the time interval, in seconds, since the first recording
        n_5min_win_in_file = cell2mat(patient_info.files2analyse(indFile2analyse,2))
        
        % starting time of first recording
        time_start_first_recording = patient_info.files(1).start_ts;
        
        % starting time of current recording
        time_start_current_recording = patient_info.files2analyse{indFile2analyse,3};
        
        % stopping time of current recording
        time_stop_current_recording = patient_info.files2analyse{indFile2analyse,4};
        
        % time interval from starting time of first recording to starting
        % time of current recording
        time_interval_since_first_rec = etime(datevec(time_start_current_recording, 'yyyy-mm-dd HH:MM:SS'), ...
            datevec(time_start_first_recording,'yyyy-mm-dd HH:MM:SS'));
        
        
        % just to check the duration of the time intervals of current
        % recording *******************************************************
        % time interval from starting time of first recording to sttoping
        % time of current recording
        time_stop_since_first_rec = etime(datevec(time_stop_current_recording, 'yyyy-mm-dd HH:MM:SS'), ...
            datevec(time_start_first_recording,'yyyy-mm-dd HH:MM:SS'));
        
        time_interval_current_recording = etime(datevec(time_stop_current_recording, 'yyyy-mm-dd HH:MM:SS'), ...
            datevec(time_start_current_recording,'yyyy-mm-dd HH:MM:SS'));
        
        check = time_interval_current_recording==(time_stop_since_first_rec-time_interval_since_first_rec);
        if check~=1
            pause
        end
        % *****************************************************************
        
        for pp = 1:length(mat_data)
            samples_win = numel(mat_data(pp).sig_preprocess);
            time_win_sec = samples_win/patient_info.samp_rate;
            time_entire_signal = linspace(time_interval_since_first_rec, ...
                time_interval_since_first_rec+time_win_sec,samples_win);
            mat_data(pp).time_entire_signal = time_entire_signal;
            time_interval_since_first_rec = time_entire_signal(end);
            % by adding the ts value as below, we are adding a ts value to
            % the time vector which will yield a time vector which is not
            % in accordance with the gap vector further on
            % ts = 1/patient_info.samp_rate; %time_entire_signal(2)-time_entire_signal(1)
            % time_interval_since_first_rec = time_entire_signal(end)+ts;
        end
        
        
        %% concatenate data in the current file to the remaining files data
        
        if ~isempty(mat_data)
            if any(ismember(fieldnames(mat_data),'detected'))
                mat_data = rmfield(mat_data, 'detected');
            end
        else
            disp('recording with no data computed')
            pause
        end
        
        if countFiles2analyse==1
            for ww = 1:length(mat_data)
                
                time_RR_intervals = horzcat(mat_data(ww).QRSindxs_recording);
                check_no_repeated_Rpeaks = isequal(numel(time_RR_intervals), ...
                    numel(unique(time_RR_intervals)));
                
                if check_no_repeated_Rpeaks==0
                    disp('line 162: problem: repeated peaks!!')
%                     pause
                end
                
                mat_data(ww).QRSindxs_240min = mat_data(ww).QRSindxs_recording;
                mat_data(ww).window_indexes_seiz = mat_data(ww).window_indexes;
                mat_data(ww).clean_segment_indexes_seiz = mat_data(ww).clean_segment_indexes;
                mat_data(ww).file_name = file_name;
            end
            
            h.noise_dataset = mat_data;
        else
            
            % the QRSindxs_recording contains the information about the
            % indexes of the R peaks since the beginning of the recording
            % In order to represent the R peak content for all recordings
            % of a given seizure it is necessary to update the indexes of
            % all but the first recording:
            
            last_window = h.noise_dataset(end).window_indexes(end);
            save_last_recording_windows = [save_last_recording_windows, last_window];
            
            for ww = 1:length(mat_data)
                vec_QRSindxs_recording = mat_data(ww).QRSindxs_recording;
                QRSindxs_240min = vec_QRSindxs_recording+sum(save_last_recording_windows);
                mat_data(ww).QRSindxs_240min = QRSindxs_240min;
                
                vec_window_indexes = mat_data(ww).window_indexes;
                window_indexes_seiz = vec_window_indexes+sum(save_last_recording_windows);
                mat_data(ww).window_indexes_seiz = window_indexes_seiz;
                
                clean_segment_indexes = mat_data(ww).clean_segment_indexes;
                if size(clean_segment_indexes,1)>1
                    a = unique(clean_segment_indexes,'rows');
                    if any(size(clean_segment_indexes)~=size(a))
%                         pause
                        clean_segment_indexes = a;
                        % b = cell2mat(cellfun(@(x) sum(ismember(clean_segment_indexes,x)),a,'un',0));
                    end
                end
                
                if ~isempty(clean_segment_indexes)
                    mat_data(ww).clean_segment_indexes_seiz = clean_segment_indexes+sum(save_last_recording_windows);
                else
                    mat_data(ww).clean_segment_indexes_seiz = [];
                end
                mat_data(ww).file_name = file_name;
            end
            
            
            mat_data_matrix = vertcat(mat_data.clean_segment_indexes_seiz);
            h_noise_dataset_matrix = vertcat(h.noise_dataset.clean_segment_indexes_seiz);
            if ~isempty(mat_data_matrix)
                if mat_data_matrix(1,1)<h_noise_dataset_matrix(end,end)
                    %mat_data(1).clean_segment_indexes_seiz(1)<h.noise_dataset(end).clean_segment_indexes_seiz(2)
                    % a linha anterior estava a dar erro quando o
                    % clean_segment_indexes_seiz era vazio na primeira janela
                    pause
                end
            end
            
            
            h.noise_dataset = horzcat(h.noise_dataset, mat_data);

        end

        % add each feature data to the output structure
        if ~isempty(indSeizure)
            
            %% Evaluate the seizure information
            
            if h.noise_dataset(end).clean_segment_indexes_seiz>numel(horzcat(h.noise_dataset.sig_preprocess))
                pause
            end
            
            
            n_5min_win_required = 48; % 240min/5min
            
            ind_n_box_5min_in_seiz_file = find(strcmp('n_box_5min_in_seiz_file', patient_info.seizures2analyseInfo(1,:)));
            n_5min_win_in_seiz_file = cell2mat(patient_info.seizures2analyseInfo(indSeizure,ind_n_box_5min_in_seiz_file));
            count_seizures = count_seizures+1;
            n_5min_win_in_seiz_file_effective = length(mat_data);
            n_5min_win_remove_end = n_5min_win_in_seiz_file_effective-n_5min_win_in_seiz_file;
            
            % because sometimes in files2analyse there are small files for
            % which 0 windows are computed for which in fact there is an
            % effective computation of signal features (even though the signal
            % in that file is smaller than 5 min --> see patient 402)
            
            n_5min_win_total = h.get_n_5min_until_seizure_file+n_5min_win_in_seiz_file;
            n_5min_win_remove_init = n_5min_win_total-n_5min_win_required;
            
            
            h.removed_samples_seiz_file = [h.removed_samples_seiz_file; ...
                numel([h.noise_dataset(end-n_5min_win_remove_end+1:end).sig_preprocess])];
            h.removed_samples_first_file = [h.removed_samples_first_file; ...
                numel([h.noise_dataset(1:n_5min_win_remove_init).sig_preprocess])];
            
            % remove the unnecessary windows starting from the fisrt
            % windows and also the last windows:
            if n_5min_win_total>=n_5min_win_required
                h.noise_dataset = h.noise_dataset(:,n_5min_win_remove_init+1:end-n_5min_win_remove_end);
            elseif n_5min_win_total<n_5min_win_required
                disp('problem: missing feature file in rec folder')
                pause
            end
            
            if any(any(vertcat(h.noise_dataset.clean_segment_indexes_seiz)<0))
                pause
            end
            
            % update clean_segment_indexes_seiz, the window_indexes_seiz
            % and the QRSindxs_240min:
            % find the first window index to subtract to the remaining
            % windows so that I can have the indexes of the clean segments
            % for the whole time of analysis (in this case for the 240 min)
            ind_first_win = h.noise_dataset(1).window_indexes_seiz(1);
            for ww = 1:length(h.noise_dataset)
                h.noise_dataset(ww).clean_segment_indexes_seiz = ...
                    h.noise_dataset(ww).clean_segment_indexes_seiz-ind_first_win+1;
                
                h.noise_dataset(ww).QRSindxs_240min = ...
                    h.noise_dataset(ww).QRSindxs_240min-ind_first_win+1;
                
                h.noise_dataset(ww).window_indexes_seiz = ...
                    h.noise_dataset(ww).window_indexes_seiz-ind_first_win+1;
            end
            % *************************************************************
            
            
            % the indexes of the clean segments must be lower than the
            % number of samples of the whole concatenated windows for each
            % seizure:
            if h.noise_dataset(end).clean_segment_indexes_seiz>numel(horzcat(h.noise_dataset.sig_preprocess))
                pause
            end
            
            if any(any(vertcat(h.noise_dataset.clean_segment_indexes_seiz)<0))
                pause
            end
            
            eval(['h.seiz' num2str(count_seizures) ' = h.noise_dataset;'])
            
            h.noise_dataset = [];
            h.get_n_5min_until_seizure_file = [];
            countFiles2analyse = 0;
            save_last_recording_windows = [];
        else
            h.get_n_5min_until_seizure_file = size(h.noise_dataset,2);
        end
        
    end
    
end


end
