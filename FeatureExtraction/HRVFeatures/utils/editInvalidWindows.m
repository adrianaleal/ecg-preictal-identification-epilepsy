function feature_dataset = editInvalidWindows(feat_names, feature_dataset)

% Here we can select the minimum window length to consider for each feature


%% put the invalid windows to NaN *************************************
% there are already some segment_size assigned to NaN...

[r,c] = size(feature_dataset);
if r~=c
    % structure is in a different format and it is necessary to concatenate
    % rows
    % structure_fields = fieldnames(feature_dataset);
    structure_fields = [feat_names; {'segment_size_seconds'; ...
        'percent_detected_noise'}];
    
    for ii = 1:numel(structure_fields)
        feature_dataset_restruct.(structure_fields{ii}) = vertcat(feature_dataset.(structure_fields{ii}));
    end
    feature_dataset = feature_dataset_restruct;
end

segment_size = feature_dataset.segment_size_seconds;

% these are the features which require a minimum of 180 s to be extracted:
features_180s = {'ApEn'};
ind_180s = ismember(feat_names, features_180s);
if any(ind_180s)
    mat_less_than_180s = segment_size<180;  
    ind_feat = find(ind_180s);
    for ss = ind_feat'
        feature_dataset.(feat_names{ss})(mat_less_than_180s) = NaN;
    end
    
end

% these are the features which require a minimum of 150 s to be extracted:
features_150s = {'NN50', 'pNN50', 'RMSSD', 'SDSD', 'RRMean', ...
    'RRMin', 'RRMax', 'RRVar', 'LF_POWER', 'HF_POWER', 'LF_NORM', ...
    'HF_NORM', 'LFtoHF', 'SD1', 'SD2', 'SD1toSD2', 'SampEn', ...
    'RQA_REC', 'RQA_DET', 'RQA_Lmax', ...
    'RQA_L', 'RQA_ENT', 'RQA_LAM', 'RQA_TT'};

ind_150s = ismember(feat_names, features_150s);
if any(ind_150s)
    mat_less_than_150s = segment_size<150;
    ind_feat = find(ind_150s);
    for ss = ind_feat'
        feature_dataset.(feat_names{ss})(mat_less_than_150s) = NaN;
    end
end



% these are the features which require a minimum of 300 s to be extracted:
percent_noise = feature_dataset.percent_detected_noise;
percent_noise_less_than_1 = percent_noise>1;

features_noisy = {'SDNN', 'TOTAL_POWER', 'VLF_POWER', ...
    'DFA_alpha1', 'DFA_alpha2', 'largLyapunov', 'corrDim'};

ind_noisy = ismember(feat_names,features_noisy);
if any(ind_noisy)
    ind_feat = find(ind_noisy);
    for ss = ind_feat'
        feature_dataset.(feat_names{ss})(percent_noise_less_than_1) = NaN;
    end
end



end
