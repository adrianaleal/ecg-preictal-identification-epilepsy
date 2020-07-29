function resultsOut = getRRInterval(resultsInput, time, sig_preprocess, ...
    fs, plotFigure)


% Information contained in resultsInput: 
% - resultsInput.window_indexes
% - resultsInput.percent_detected_noise
% - resultsInput.clean_segment_indexes
% - resultsInput.QRSindxs_240min



% Information contained in resultsOutput:
results_names = {'long_window', 'time', 'segment_size_seconds', ...
    'n_clean_segments', 'percent_detected_noise', 'nHB', ...
    'RRI_240min_corrected', 'time_RRI_240min_corrected'}';

n_wins = length(resultsInput);

for ff = 1:length(results_names)
    [resultsOut(1:n_wins).(results_names{ff})] = deal(NaN);
end


for ii = 1:n_wins % [8, 10]%  % go through each 5 min window
    
    %     disp(['>> Analyzing window ' num2str(ii)])
    fprintf(2,['\n>>>>>> Analyzing window ' num2str(ii) ' <<<<<<<<<<<<<<<\n'])
    
    ind_win1 = resultsInput(ii).window_indexes(1);
    ind_win2 = resultsInput(ii).window_indexes(2);
    time_win = time(ind_win1:ind_win2);
    
    resultsOut(ii).time = time_win(round(numel(time_win)/2));
    resultsOut(ii).percent_detected_noise = resultsInput(ii).percent_detected_noise;
    resultsOut(ii).long_window = 0;
    
    
    clean_segment_indexes = resultsInput(ii).clean_segment_indexes;
    % when the window has no R wave peaks:
    if isempty(clean_segment_indexes)
        resultsOut(ii).n_clean_segments = 0;
        continue
    end
    
    
%     if ii==2102
%         disp('')
%     end
                
    
    %%
    
    Rwave_indxs_240min_win = resultsInput(ii).QRSindxs_240min;    
    
    % results.QRSindxs_240min is obtained in define_clean_segments.m
    % this structure differs from the results.QRSindxs_entire_sig in the
    % sense that in the first, the location of the R peaks found in noisy
    % segments are assigned as NaN --> this will not affect the analysis in
    % feature extraction step as only the clean segments are analysed in
    % this step.
    
%     if plotFigure
%         segment = resultsInput(ii).sig_preprocess;
%         figure()
%         hax1 = axes;
%         plot(time_win, segment,'r'); hold on
%         axis tight
%         xlabel('Time (s)')
%         ylabel('mV')
%     end
    
    n_clean_segments = size(clean_segment_indexes,1);
    concatenated_Rpeak_segment = [];
    segment_size_seconds_vec = [];
    time_Rwave_segment_last = 0;
    

    for ll = 1:n_clean_segments
        ind_clean1 = clean_segment_indexes(ll,1);
        ind_clean2 = clean_segment_indexes(ll,2);
        ind1_Rwave1 = find(Rwave_indxs_240min_win>=ind_clean1);
        ind1_Rwave2 = find(Rwave_indxs_240min_win<=ind_clean2);
        
        if isempty(ind1_Rwave1) || isempty(ind1_Rwave2)
            continue
        end
        
        if ind1_Rwave1(1)>ind1_Rwave2(end)
            disp('the clean segment is located in between a RR interval!')
            continue
        end
        
     
        segment_size = length(sig_preprocess(ind_clean1:ind_clean2));
        
        segment_size_seconds = segment_size/fs;
        if isnan(segment_size_seconds)
            disp('see')
            pause
        end
        fprintf(2,['\n>>>>>> Analyzing segment ' num2str(ll) ' of length '...
            num2str(segment_size_seconds) ' secs \n'])
        
        if segment_size_seconds>0
            Rwave_indxs_win_segment = Rwave_indxs_240min_win(ind1_Rwave1(1):ind1_Rwave2(end));
            
            time_Rwave_segment = time(Rwave_indxs_win_segment);
            
            %% Plot segment
%             if plotFigure
%                 plot(hax1,time(ind_clean1:ind_clean2),sig_preprocess(ind_clean1:ind_clean2),'b')
%                 plot(hax1,time_Rwave_segment,sig_preprocess(Rwave_indxs_win_segment),'g*')
%             end

            % concatenate the time vector of the heart beats
            % by eliminating the interval between the end and start of
            % subsequent segments, resulting in overlapping two heart beats
            % and eliminating one of them
            if numel(time_Rwave_segment)>3
                
                % because sometimes at the beginning there are wrong peak
                % detection: **************************************************
                RR_interval_signal_segment = time_Rwave_segment(2:end)-time_Rwave_segment(1:end-1);
                ind_small_RR = RR_interval_signal_segment<0.4; 
                indexes = find(ind_small_RR);
                if ~isempty(indexes)
                    continuous_RR = isequal(indexes, indexes(1):indexes(end));
                    if continuous_RR
                        time_Rwave_segment = time_Rwave_segment(~ind_small_RR);
                    end
                end
                % typically peaks are too close to each other
                % *************************************************************
            
                RR_interval_signal_segment = time_Rwave_segment(2:end)-time_Rwave_segment(1:end-1);
                ind_large_RR = RR_interval_signal_segment>1.9;
                indexes = find(ind_large_RR);
                if ~isempty(indexes) && sum(ind_large_RR)>1
                    continuous_RR = isequal(indexes, indexes(1):indexes(end));
                    if continuous_RR
                        time_Rwave_segment = time_Rwave_segment(~ind_large_RR);
                    end
                end
                % typically peaks are too far from each other
                % *************************************************************
   
               
                if ll>1
                    RR_interval_signal_segment = diff(time_Rwave_segment);
                    time_Rwave_segment = time_Rwave_segment_last+cumsum(RR_interval_signal_segment);
                end

                concatenated_Rpeak_segment = [concatenated_Rpeak_segment, time_Rwave_segment];

                time_Rwave_segment_last = concatenated_Rpeak_segment(end);
                segment_size_seconds_vec = [segment_size_seconds_vec; segment_size_seconds];
                
            end
        end
    end
    
    resultsOut(ii).n_clean_segments = n_clean_segments;

    if numel(concatenated_Rpeak_segment)>3

        concatenated_segments_size_seconds = sum(segment_size_seconds_vec);
        if concatenated_segments_size_seconds>300 % 5 min
            disp('AVISO --> análise de sinais maiores que 300 secs')
            
            % concatenated_segments_size_seconds = 300;
            % time_win_300_segment = time_win(1:concatenated_segments_size_seconds*fs);
            % resultsOut(ii).time = time_win_300_segment(round(numel(time_win_300_segment)/2));
            % ind_end = find(concatenated_Rpeak_segment>time_win_300_segment(end), 1, 'first');
            % concatenated_Rpeak_segment = concatenated_Rpeak_segment(1:ind_end);
            
            % noise_long_window = resultsInput(ii).detected(1:concatenated_segments_size_seconds*fs);
            % percent_detected_noise_300secs = numel(find(noise_long_window))/numel(noise_long_window)*100;
%             resultsOut(ii).percent_detected_noise = percent_detected_noise_300secs;
            resultsOut(ii).long_window = 1;
        end
        
        
        resultsOut(ii).segment_size_seconds = concatenated_segments_size_seconds;
        resultsOut(ii).time_RRI_240min_corrected = concatenated_Rpeak_segment;
        resultsOut(ii).RRI_240min_corrected = diff(concatenated_Rpeak_segment);
        % resultsOut(ii).abnormalRR = nAbnormalHB;
        resultsOut(ii).nHB = numel(concatenated_Rpeak_segment);
        
%         if resultsOut(ii).abnormalRR/nHB>0.8
%             disp('')
%         end
        
%         figure()
%         plot(time_RR_intervals_original, RR_interval_signal_segment_original)
%         hold on
%         plot(time_RR_intervals_corrected, RR_interval_signal_segment_corrected, '--r')
%         hold off
%         axis tight
%         legend('RR intervals', 'RR intervals corrected')
%         ylabel('ms')
%         xlabel('Time (hour)')
%         
    end
    
end

% [features(:).time] = deal(results.time);


end