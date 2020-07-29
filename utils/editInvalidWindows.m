function feature_dataset = editInvalidWindows(feat_names, feature_dataset)


 %% put the invalid windows to NaN *************************************
 % there are already some segment_size assigned to NaN...
 
[r,c] = size(feature_dataset); 
if r~=c
    % structure is in a different format and it is necessary to concatenate
    % rows
%     structure_fields = fieldnames(feature_dataset);
    structure_fields = [feat_names; {'segment_size_seconds'; ...
        'percent_detected_noise'}];
    
    for ii = 1:numel(structure_fields)
        feature_dataset_restruct.(structure_fields{ii}) = vertcat(feature_dataset.(structure_fields{ii}));
    end
    feature_dataset = feature_dataset_restruct;
end
     
%  abnormalRR1 = feature_dataset.abnormalRR;
 segment_size = feature_dataset.segment_size_seconds;

 ind_ApEn = strcmp(feat_names,'ApEn');
 if any(ind_ApEn)
     mat_less_than_180s = segment_size<180;
     feature_dataset.ApEn(mat_less_than_180s) = NaN;
 end
 
 
 features_150s = {'NN50', 'pNN50', 'RMSSD', 'SDSD', 'RRMean', ...
     'RRMin', 'RRMax', 'RRVar', 'LF_POWER', 'HF_POWER', 'LF_NORM', ...
     'HF_NORM', 'LFtoHF', 'SD1', 'SD2', 'SD1toSD2', 'SampEn', ...
     'corrDim', 'largLyapunov', 'RQA_REC', 'RQA_DET', 'RQA_Lmax', ...
     'RQA_L', 'RQA_ENT', 'RQA_LAM', 'RQA_TT'};
 
 ind_other = ismember(feat_names,features_150s);
 if any(ind_other)
     mat_less_than_150s = segment_size<150;
     ind_feat = find(ind_other);
     for ss = ind_feat'
%          feat_names{ss}
         feature_dataset.(feat_names{ss})(mat_less_than_150s) = NaN;
     end
 end
 
 
 
 
 percent_noise = feature_dataset.percent_detected_noise;
 percent_noise_less_than_1 = percent_noise>1;
 
 features_noisy = {'TOTAL_POWER', 'VLF_POWER', 'SDNN', ...
     'DFA_alpha1', 'DFA_alpha2'};
 
 ind_noisy = ismember(feat_names,features_noisy);
 if any(ind_noisy)
     ind_feat = find(ind_noisy);
     for ss = ind_feat'
%          feat_names{ss}
         feature_dataset.(feat_names{ss})(percent_noise_less_than_1) = NaN;
     end
 end
            


end
