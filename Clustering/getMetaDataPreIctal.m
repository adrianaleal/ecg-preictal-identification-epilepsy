function [var2correlate, matrix_patient_info_classes, n_seiz, ...
    matrix_patient_info_out, mat4AssociationLearning, metaData, ...
    store_seizures_start_date] = ...
    getMetaDataPreIctal(outer_folder, patients_info_final)

%% for the correlation between the results and the type of seizure and
% vigilance state

n_seiz_vec = vertcat(patients_info_final.n_seizures)- ...
    vertcat(patients_info_final.n_discarded_seizures_by_240min)- ...
    vertcat(patients_info_final.n_discarded_seizures_by_no_file);
n_seiz = sum(n_seiz_vec);



n_max_seizures = max(n_seiz_vec);
n_pat = size(patients_info_final,1);
seiz_type_mat = repmat({''},n_pat,n_max_seizures);% cell(n_pat,n_max_seizures);
vigilance_mat = seiz_type_mat;
hour_mat = seiz_type_mat;
pattern_mat = seiz_type_mat;
medication_mat = seiz_type_mat;
dosage_mat = seiz_type_mat;


%%
matrix_patient_info = cell(n_seiz,12);
start = 1;
store_seizures_start_date = cell(n_seiz,2);
for pp = 1:n_pat
    n_seiz_in_patient = size(patients_info_final(pp).seizures2analyseInfo,1)-1;
    stop = start+n_seiz_in_patient-1;
    matrix_patient_info(start:stop,1) = cellstr(repmat(patients_info_final(pp).name,n_seiz_in_patient,1));
    matrix_patient_info(start:stop,2:7) = patients_info_final(pp).seizures2analyseInfo(2:end,[1,3,5:7,end]);
    
    seiz_type_mat(pp,1:n_seiz_in_patient) = (patients_info_final(pp).seizures2analyseInfo(2:end,end))';
    vigilance_mat(pp,1:n_seiz_in_patient) = (patients_info_final(pp).seizures2analyseInfo(2:end,3))';
    pattern_mat(pp,1:n_seiz_in_patient) = (patients_info_final(pp).seizures2analyseInfo(2:end,5))';
    
    store_seizures_start_date(start:stop,1) = matrix_patient_info(start:stop,1);
    
    %% Medication and dosage
    meds = patients_info_final(pp).seizures2analyseInfo(2:end,6);
    newStr = cellfun(@(s) strsplit(s, ' '), meds, 'UniformOutput', false);
    try
        % if is numeric
        meds_numeric = str2double(vertcat(newStr{:}));
        if size(unique(meds_numeric,'rows'),1)==2 && numel(unique(sum(meds_numeric,2)))
            % disp('seizures with same medication but in different order')
            meds = repmat(meds(1,:),n_seiz_in_patient,1);
        end
    catch
        % it is No info
    end
    
    medication_mat(pp,1:n_seiz_in_patient) = meds';
    dosage_mat(pp,1:n_seiz_in_patient) = (patients_info_final(pp).seizures2analyseInfo(2:end,7))';
    
    %% Seizure onset hour
    seizures_start_date = vertcat(patients_info_final(pp).seizures2analyseInfo(2:end,1));
    C = cellfun(@(x)strsplit(x, ' ' ), seizures_start_date, 'UniformOutput', false);
    seiz_hour = vertcat(C{:});
    seiz_onset_hour = cellfun(@hour,seiz_hour(:,2),'uni',false);
    hour_mat(pp,1:n_seiz_in_patient) = seiz_onset_hour;
    matrix_patient_info(start:stop,8) = cellstr(num2str(cell2mat(seiz_onset_hour)));
    
    store_seizures_start_date(start:stop,2) = seizures_start_date;
    
    %% Information regarding the temporal distance between the seizures
    
    all_seizures_start_date = vertcat(patients_info_final(pp).AllSeizuresInfo(2:end,1));
    
    x_label = cellstr(num2str((1:numel(all_seizures_start_date))'));
    vec_names = ['start recording'; strcat('seiz',x_label)];
    vec_names_final = strcat(vec_names(1:end-1), {':'}, vec_names(2:end));
    
    first_recording_start_date = {patients_info_final(pp).files(1).start_ts};
    
    dates_vec = [first_recording_start_date; all_seizures_start_date];
    time_vec = datevec(datenum(dates_vec));
    
    all_seizure_difference_in_secs = etime(time_vec(2:end,:),time_vec(1:end-1,:));
    
    
    all_seizure_difference_in_date = cellstr(datestr(datenum(0,0,0,0,0,all_seizure_difference_in_secs), 'DD HH:MM:SS'));
    
    ind_analysed_seizures = ismember(all_seizures_start_date,seizures_start_date);
    all_seizure_difference_in_date = all_seizure_difference_in_date(ind_analysed_seizures);
    matrix_patient_info(start:stop,9) = num2cell(all_seizure_difference_in_secs(ind_analysed_seizures));
    matrix_patient_info(start:stop,10) = all_seizure_difference_in_date;
    
    matrix_patient_info(start:stop,11) = num2cell(vec_names_final(ind_analysed_seizures));
    
    %%
    
    matrix_patient_info(start:stop,end) = num2cell((1:n_seiz_in_patient)');
    
    start = stop+1;
    
    
end

metaData.seiz_type_mat = seiz_type_mat;
metaData.vigilance_mat = vigilance_mat;
metaData.hour_mat = hour_mat;
metaData.pattern_mat = pattern_mat;
metaData.medication_mat = medication_mat;
metaData.dosage_mat = dosage_mat;


%% in order to know information about the discarded seizures based on the
% missing values
filename = 'indexes_seizures_discarded_by_MVs';
if exist(fullfile(outer_folder,[filename '.mat']), 'file') == 2
    % File exists.
    load(fullfile(outer_folder,filename))
    disp('Seizures removed: ')
    ind_seiz_in_pat = matrix_patient_info(indexes_seizures_discarded_by_MVs,4);
    patient = matrix_patient_info(indexes_seizures_discarded_by_MVs,1);
    vigilance_state = matrix_patient_info(indexes_seizures_discarded_by_MVs,2);
    seizure_type = matrix_patient_info(indexes_seizures_discarded_by_MVs,3);
    
    table(ind_seiz_in_pat, patient, vigilance_state, seizure_type)
    matrix_patient_info(indexes_seizures_discarded_by_MVs,:) = [];
end


% all variables but the index of the seizures for each patient:
matrix_patient_info_out = matrix_patient_info(:,1:end-1);

% the categorical variables to turn into numeric ones:
matrix_patient_info = matrix_patient_info(:,[1,3,4,7,8]);

%% Information for all patients, seizure types and vigilance state
matrix_patient_info_classes = matrix_patient_info;
ind_vigilance = 2;
matrix_patient_info_classes(strcmp(matrix_patient_info_classes(:,ind_vigilance),'A'),ind_vigilance) = {'3'};
matrix_patient_info_classes(strcmp(matrix_patient_info_classes(:,ind_vigilance),'R'),ind_vigilance) = {'4'};
% matrix_patient_info_classes(strcmp(matrix_patient_info_classes(:,2),'?'),2) = {'6'};

ind_seiz_type = 4;
[seiz_type_names, ~, seiz_type_discrete] = unique(matrix_patient_info_classes(:,ind_seiz_type));
matrix_patient_info_classes(:,ind_seiz_type) = cellstr(num2str(seiz_type_discrete));

ind_pattern = 3;
[pattern_names, ~, pattern_discrete] = unique(matrix_patient_info_classes(:,ind_pattern));
matrix_patient_info_classes(:,ind_pattern) = cellstr(num2str(pattern_discrete));

% seizure onset hour
ind_onset_hour = 5;
hour_cell = str2double(cellstr(matrix_patient_info(:,ind_onset_hour)));
ind_day = hour_cell>=7|hour_cell<22;
ind_night1 = hour_cell<7;
ind_night2 = hour_cell>=22;
hour_day_night(ind_day) = 1;
hour_day_night(ind_night1) = 2;
hour_day_night(ind_night2) = 2;
matrix_patient_info_classes(:,ind_onset_hour) = cellstr(num2str(hour_day_night'));

% assign classes to the patients
split1 = regexp(matrix_patient_info(:,1), '_', 'split');
get_pat_ids = vertcat(split1{:});
[B,I] = sort(cellfun(@str2double,get_pat_ids(:,2)));
[patient_names, ~, patient_discrete] = unique(B);
matrix_patient_info_classes(:,1) = cellstr(num2str(patient_discrete));
matrix_patient_info_classes = cell2mat(cellfun(@str2num,matrix_patient_info_classes,'un',0));

ind_seiz_240_0min_before_seiz_preictal = matrix_patient_info;
save('ind_seiz_240_0min_before_seiz_preictal.mat','ind_seiz_240_0min_before_seiz_preictal')

% in order to inspect if correlation exists between the possible groups
% among the variables: vigilance_state w/ seiz_type, vigilance_state w/
% patient, patient w/ seiz_type and also the whole possible combination
% vigilance_state w/ seiz_type w/ patient
combine_patient_VS = strcat(matrix_patient_info(:,1), {'_'}, matrix_patient_info(:,2));
[patient_VS_names, ~, patient_VS_discrete] = unique(combine_patient_VS);

combine_patient_ST = strcat(matrix_patient_info(:,1), {'_'}, matrix_patient_info(:,3));
[patient_ST_names, ~, patient_ST_discrete] = unique(combine_patient_ST);

combine_VS_ST = strcat(matrix_patient_info(:,2), {'_'}, matrix_patient_info(:,3));
[VS_ST_names, ~, VS_ST_discrete] = unique(combine_VS_ST);

combine_patient_VS_ST = strcat(matrix_patient_info(:,1), {'_'}, matrix_patient_info(:,2), {'_'}, matrix_patient_info(:,3));
[patient_VS_ST_names, ~, patient_VS_ST_discrete] = unique(combine_patient_VS_ST);
matrix_patient_info_classes = [matrix_patient_info_classes, patient_VS_discrete, ...
    patient_ST_discrete, VS_ST_discrete, patient_VS_ST_discrete];

var2correlate = {'patient','vigilance_state','seiz_type', 'patient_VS', ...
    'patient_ST', 'VS_ST', 'patient_VS_ST'};


%% compute the binary matrix to give to the classifier

x_label = cellstr(num2str((1:n_max_seizures)'));


% number of seizures
seiz_number = num2str(unique(n_seiz_vec));
vec_names_NS = strcat({'NoSeiz_'}, seiz_number);

% vigilance state
VS_classes = unique(matrix_patient_info(:,ind_vigilance));
VS = repmat(VS_classes,1,n_max_seizures)';
vec_names_VS = strcat(repmat(strcat('seiz',x_label), numel(VS_classes),1),{'_VS_'}, VS(:));

% seizure type
ST_classes = unique(matrix_patient_info(:,ind_seiz_type));
ST = repmat(ST_classes,1,n_max_seizures)';
vec_names_ST = strcat(repmat(strcat('seiz',x_label), numel(ST_classes),1),{'_ST_'}, ST(:));


% seizure pattern
pattern_classes = unique(matrix_patient_info(:,ind_pattern));
pattern = repmat(pattern_classes,1,n_max_seizures)';
vec_names_pattern = strcat(repmat(strcat('seiz',x_label), ...
    numel(pattern_classes),1),{'_pattern_'}, pattern(:));


hour_classes_day_night = cellstr(num2str(unique(hour_day_night)'));
onset_hour = repmat(hour_classes_day_night,1,n_max_seizures)';
vec_names_hour = strcat(repmat(strcat('seiz',x_label), ...
    numel(hour_classes_day_night),1),{'_onset_hour_'}, onset_hour(:));


binary_matrix_names = [vec_names_NS' vec_names_VS' vec_names_ST' ...
    vec_names_pattern' vec_names_hour'];
binary_matrix = zeros(n_pat,numel(binary_matrix_names));

start = 1;
for pp = 1:n_pat
    n_seiz_in_patient = size(patients_info_final(pp).seizures2analyseInfo,1)-1;
    get_NS = strcat({'NoSeiz_'}, num2str(n_seiz_in_patient));
    ind_var = strcmp(get_NS, binary_matrix_names);
    binary_matrix(pp,find(ind_var)) = 1;
    stop = start+n_seiz_in_patient-1;
    ind_seiz = start:stop;
    for ss = 1:n_seiz_in_patient
        get_VS = matrix_patient_info(ind_seiz(ss), ind_vigilance);
        ind_var = strcmp(strcat('seiz', num2str(ss), '_VS_',get_VS), binary_matrix_names);
        binary_matrix(pp,find(ind_var)) = 1;
        
        get_ST = matrix_patient_info(ind_seiz(ss), ind_seiz_type);
        ind_var = strcmp(strcat('seiz', num2str(ss), '_ST_',get_ST), binary_matrix_names);
        binary_matrix(pp,find(ind_var)) = 1;
    end
    start = stop+1;
end

mat4AssociationLearning = array2table(binary_matrix,'VariableNames', ...
    binary_matrix_names);

end