function [featuresOut, comp_time] = getHRVFeatures(resultsRR, ...
    resultsOutput, plotFigure)


features_names = {'NN50', 'pNN50', 'SDNN', 'RMSSD', 'SDSD', 'RRMean', ...
    'RRMin', 'RRMax', 'RRVar', 'TOTAL_POWER', 'VLF_POWER', 'LF_POWER', ...
    'HF_POWER', 'VLF_NORM', 'LF_NORM', 'HF_NORM', 'LFtoHF', 'SD1', ...
    'SD2', 'SD1toSD2', 'ApEn', 'SampEn', 'DFA_alpha1', 'DFA_alpha2', ...
    'tau_rec_phase_space','embDim_rec_phase_space', 'corrDim', ...
    'largLyapunov', 'RQA_REC', 'RQA_DET', 'RQA_Lmax', 'RQA_L', ...
    'RQA_ENT', 'RQA_LAM', 'RQA_TT'}';

for ff = 1:length(features_names)
    [features(1:length(resultsRR)).(features_names{ff})] = deal(NaN);
end


% to analyse computational time
CP_features_names = {'time', 'frequency', 'poincare', 'ApEn', 'SampEn', ...
    'DFA', 'corrDim', 'largLyapunov', 'RQA'}';

for ff = 1:length(CP_features_names)
    [comp_time(1:length(resultsRR)).(CP_features_names{ff})] = deal(NaN);
end



for ii = 1:length(resultsRR) % [8, 10]%  % go through each 5 min window
    
    %     disp(['>> Analyzing window ' num2str(ii)])
    fprintf(2,['\n>>>>>> Analyzing window ' num2str(ii) ' <<<<<<<<<<<<<<<\n'])
    
    if isempty(resultsOutput(ii).clean_segment_indexes)% when the window has no R wave peaks
        continue
    end
    

    %%
    RR_interval_signal_segment = resultsRR(ii).RRI_240min_corrected;
    if ~isempty(RR_interval_signal_segment)

%         if (resultsRR(ii).nAbnormalRR/resultsRR(ii).nHB)>0.8
%             continue
%         end

        time_RR_intervals = resultsRR(ii).time_RRI_240min_corrected(2:end);
        
        concatenated_segments_size_seconds = resultsRR(ii).segment_size_seconds;
        ll = 1;
        [features(ii), comp_time(ii)] = HRV_features(RR_interval_signal_segment*1e3, ...
            time_RR_intervals, features(ii), comp_time(ii), ll, ...
            concatenated_segments_size_seconds, plotFigure);
        
        
%         if ii==1800
%             disp('')
%         end
    end
    
end


resultsRR_fieldnames = fieldnames(resultsRR);
resultsRR2Feat = rmfield(resultsRR, resultsRR_fieldnames(end-1:end));
featuresOut = table2struct([struct2table(resultsRR2Feat), struct2table(features)]);


end