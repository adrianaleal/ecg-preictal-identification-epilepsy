function [resultsOutput, save2FigurePlotECG] = getSignalWithOverlapSlidingWindow( ...
    seiz_noise_det, time, gap_logical, gap_time, fs, plotFigure, file_name)

% Information contained in seiz_noise_det: 
% - seiz_noise_det.sig_preprocess
% - seiz_noise_det.QRSindxs_240min
% - seiz_noise_det.clean_segment_indexes_seiz
% - seiz_noise_det.window_indexes_seiz



% Information contained in resultsOutput:
% - resultsOutput.sig_preprocess
% - resultsOutput.percent_detected_noise
% - resultsOutput.window_indexes
% - resultsOutput.clean_segment_indexes
% - resultsOutput.QRSindxs_240min



% Note:
% the QRSindxs_entire_sig_seiz contains the information of 
% QRSindxs_entire_sig for the entire 240min
% the QRSindxs_entire_sig contains the information about the
% indexes of the R peaks for each recording since its beginning 
% this variable is obtained in GetAllFilesInDir4FeatureExtraction.m

window_indexes_seiz  = vertcat(seiz_noise_det.window_indexes_seiz);


close all

if any(gap_logical)
    % find the time indexes of the gaps existent in window:
    gap_samples_start_indexes = strfind(gap_logical',[0 1])';
    gap_samples_end_indexes = strfind(gap_logical',[1 0])';
    if gap_logical(1)==1
        gap_samples_start_indexes = [1; gap_samples_start_indexes];
    end
    if gap_logical(end)==1
        gap_samples_end_indexes = [gap_samples_end_indexes; length(gap_logical)];
    end
    
    % gap_samples = [gap_samples_start_indexes+1 gap_samples_end_indexes]
    
    indxs_no_gap_segments = [gap_samples_end_indexes(1:end-1)+1, ...
        gap_samples_start_indexes(2:end)];
    
    if gap_samples_start_indexes(1)~=1
        indxs_no_gap_segments = [1 gap_samples_start_indexes(1); ...
            indxs_no_gap_segments];
    end
    
    if gap_samples_end_indexes(end)~=length(gap_logical)
        indxs_no_gap_segments = [indxs_no_gap_segments;...
            gap_samples_end_indexes(end)+1 length(gap_logical)];
    end
    
    
    
    % get the no gap segments in terms of the time vector instead of the
    % gap_time:
    
    no_gap_segments_ingap_time = gap_time(indxs_no_gap_segments);
    ind_no_gap_segments_intime = zeros(size(no_gap_segments_ingap_time));
    window_indexes_seiz_time = time(window_indexes_seiz);
    
    for gg = 1:size(no_gap_segments_ingap_time)
        ind1 = gap_time(indxs_no_gap_segments(gg,1));
        ind2 = gap_time(indxs_no_gap_segments(gg,2));
        
        ind_clean1 = window_indexes_seiz_time(:,1)>=ind1 & ...
            window_indexes_seiz_time(:,1)<=ind2;
        ind_clean2 = window_indexes_seiz_time(:,2)>=ind1 & ...
            window_indexes_seiz_time(:,2)<=ind2;
        
        ver = [ind_clean1 ind_clean2];
        
        vec = find(sum(ver,2)>0);
        if ~isempty(vec)
            ind_no_gap_segments_intime(gg,:) = [window_indexes_seiz(vec(1),1) ...
                window_indexes_seiz(vec(end),2)];
        end
    end
    
else
    ind_no_gap_segments_intime = [window_indexes_seiz(1,1) ...
                window_indexes_seiz(end,2)];
end
% *************************************************************************
%% get the variables of interest:
sig_preprocess_seiz = [seiz_noise_det.sig_preprocess];
% sig_preprocess_seiz = horzcat(resultsInput.sig_preprocess);

QRSindxs_240min = [seiz_noise_det.QRSindxs_240min];
QRSindxs_240min = unique(QRSindxs_240min);



clean_segment_indexes_seiz = vertcat(seiz_noise_det.clean_segment_indexes_seiz);

% probably due to bad annotation, sometimes differences of one sample are
% found, and clean_segment_indexes_seiz(:,2) are smaller than
% clean_segment_indexes_seiz(:,1) which is corrected below:
clean_diff = clean_segment_indexes_seiz(:,2)-clean_segment_indexes_seiz(:,1);
indexes2remove = clean_diff==-1 | clean_diff==0;
clean_segment_indexes_seiz(indexes2remove,:) = [];

clean_sig = ones(size(sig_preprocess_seiz));
for cc = 1:size(clean_segment_indexes_seiz,1)
    clean_sig(clean_segment_indexes_seiz(cc,1):clean_segment_indexes_seiz(cc,2)) = 0;
end

% the following code is meant to concatenate the clean_segments that are
% close to each because of the 5min non overlapped previous analysis

if any(clean_sig)
    noisy_segments_start_indexes = strfind(clean_sig,[0 1])';
    noisy_segments_end_indexes = strfind(clean_sig,[1 0])';
    if clean_sig(1)==1
        noisy_segments_start_indexes = [1; noisy_segments_start_indexes];
    end
    if clean_sig(end)==1
        noisy_segments_end_indexes = [noisy_segments_end_indexes; length(clean_sig)];
    end
    % noisy_segments = [noisy_segments_start_indexes+1 noisy_segments_end_indexes]
    
    clean_segment_indexes_seiz = [noisy_segments_end_indexes(1:end-1)+1, ...
        noisy_segments_start_indexes(2:end)];
    
    if noisy_segments_start_indexes(1)~=1
        clean_segment_indexes_seiz = [1 noisy_segments_start_indexes(1); ...
            clean_segment_indexes_seiz];
    end
    
    if noisy_segments_end_indexes(end)~=numel(sig_preprocess_seiz)
        clean_segment_indexes_seiz = [clean_segment_indexes_seiz;...
            noisy_segments_end_indexes(end)+1 numel(sig_preprocess_seiz)];
    end
     
else
    clean_segment_indexes_seiz = [1 numel(sig_preprocess_seiz)];
end

%% plot Figure all seizure to save

% inmin = min(sig_preprocess_seiz(:));
% inmax = max(sig_preprocess_seiz(:));
% sig_preprocess_seiz = l + [(sig_preprocess_seiz-inmin)./(inmax-inmin)].*(u-l)
% sig_preprocess_seiz = (sig_preprocess_seiz - inmin) / (inmax - inmin);


save2FigurePlotECG.sig_preprocess_seiz = sig_preprocess_seiz;
save2FigurePlotECG.time_recording = time;
save2FigurePlotECG.clean_sig = clean_sig;
save2FigurePlotECG.QRSindxs_240min = QRSindxs_240min;
save2FigurePlotECG.gap_logical = gap_logical;
save2FigurePlotECG.window_indexes_seiz = window_indexes_seiz;


%% -----backSearch---------------------------------------------------------
QRSindxs_240min = RpeaksBackSearch(QRSindxs_240min, sig_preprocess_seiz, fs);

save2FigurePlotECG.QRSindxs_240min_original = QRSindxs_240min;

%% -----RR interval series editing-----------------------------------------
% edit only clean RRI_series (outside noisy and gap intervals)

% this was considered instead of performing correction in each overlapped
% 5-min window because the same correction was being made many times, and
% when the segment that was often corrected was in the end of a window a
% bad correction was performed





if plotFigure
   figure(124)
   plot(time, sig_preprocess_seiz), hold on
   plot(time, clean_sig*(2/3)*max(sig_preprocess_seiz), 'r')
   plot(gap_time, gap_logical*(2/3)*max(sig_preprocess_seiz), 'c')
   axis tight
end



RRI_original_final = {};
RRI_corrected_final = {};
time_Rpeaks_final = {};
time_Rpeaks_corrected_final = {};


% RRI_corrected_final = [];
% time_Rpeaks_corrected_final = [];
ind_Rpeaks_corrected_final = [];
% save_indAbnormalRR = [];

for gg = 1:size(ind_no_gap_segments_intime,1)
    
    ind_clean1 = clean_segment_indexes_seiz(:,1)>=ind_no_gap_segments_intime(gg,1) & ...
        clean_segment_indexes_seiz(:,1)<=ind_no_gap_segments_intime(gg,2);
    
    ind_clean2 = clean_segment_indexes_seiz(:,2)>=ind_no_gap_segments_intime(gg,1) & ...
        clean_segment_indexes_seiz(:,2)<=ind_no_gap_segments_intime(gg,2);
    
    ver = [ind_clean1 ind_clean2];
    ind_clean_segments = clean_segment_indexes_seiz(sum(ver,2)==2,:);

    
    if any(sum(ver,2)==1)
        ind_row = find(sum(ver,2)==1);
        ind_col = ver(ind_row,:);
        
        for pp = 1:numel(ind_row)
            ind_clean_segments = [ind_clean_segments; ...
                sort([clean_segment_indexes_seiz(ind_row(pp),ind_col(pp,:)) ...
                ind_no_gap_segments_intime(gg,~ind_col(pp,:))])];
        end
        ind_clean_segments = sort(ind_clean_segments);
    end

    
    % check if the no_gap_segment is inside a clean segment
    if isempty(ind_clean_segments)
        for cc = 1:size(clean_segment_indexes_seiz,1)
            if ind_no_gap_segments_intime(gg,1)>=clean_segment_indexes_seiz(cc,1) && ...
                    ind_no_gap_segments_intime(gg,2)<=clean_segment_indexes_seiz(cc,2)
                ind_clean_segments = [ind_no_gap_segments_intime(gg,1) ...
                    ind_no_gap_segments_intime(gg,2)];     
            end
        end
    end

    for cc = 1:size(ind_clean_segments,1)
        
        ind_start_clean = ind_clean_segments(cc,1);
        ind_end_clean = ind_clean_segments(cc,2);
        
        
        if plotFigure
            figure(124)
            plot(time(ind_start_clean:ind_end_clean), ...
                sig_preprocess_seiz(ind_start_clean:ind_end_clean)), hold on
        end
        
        ind_QRSindxs_clean = find(QRSindxs_240min>=ind_start_clean & ...
            QRSindxs_240min<=ind_end_clean);
        
        if numel(ind_QRSindxs_clean)>2
            RRIs_clean_segment = QRSindxs_240min(ind_QRSindxs_clean);

            [RRI_signal_segment_original, RRI_signal_segment_corrected, ...
                time_RRIs_original, time_RRIs_corrected, indAbnormalRR, ind_Rpeaks, ...
                ] = RRI_series_editing(time, RRIs_clean_segment, ...
                sig_preprocess_seiz, ind_start_clean, ind_end_clean, fs, ...
                plotFigure);
    %         RRI_corrected_final = [RRI_corrected_final RRI_signal_segment_corrected];
    %         time_Rpeaks_corrected_final = [time_Rpeaks_corrected_final time_RRIs_corrected];
    %         ind_Rpeaks_corrected_final = [ind_Rpeaks_corrected_final ind_Rpeaks];

            RRI_original_final = [RRI_original_final; RRI_signal_segment_original];
            RRI_corrected_final = [RRI_corrected_final; {RRI_signal_segment_corrected}];

            time_Rpeaks_final = [time_Rpeaks_final; {time_RRIs_original}];
            time_Rpeaks_corrected_final = [time_Rpeaks_corrected_final; {time_RRIs_corrected}];

            ind_Rpeaks_corrected_final = [ind_Rpeaks_corrected_final ind_Rpeaks];
            
            % save_indAbnormalRR = [save_indAbnormalRR; indAbnormalRR];
        end
    end
    
end

% save_indAbnormalRR = sort(save_indAbnormalRR);
QRSindxs_240min = ind_Rpeaks_corrected_final;


save2FigurePlotECG.QRSindxs_240min_corrected = ind_Rpeaks_corrected_final;

save2FigurePlotECG.time_Rpeaks_240min = time_Rpeaks_final;
save2FigurePlotECG.time_Rpeaks_240min_corrected = time_Rpeaks_corrected_final;

save2FigurePlotECG.RRI_240min = RRI_original_final;
save2FigurePlotECG.RRI_240min_corrected = RRI_corrected_final;

%% -----PLOT FIGURE--------------------------------------------------------

if plotFigure
    y = [min(sig_preprocess_seiz) max(sig_preprocess_seiz)];
    figure(101)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    h1 = plot(time, sig_preprocess_seiz); hold on
    h2 = plot(time, clean_sig*(1/3)*max(sig_preprocess_seiz));
    h3 = plot(time(QRSindxs_240min), ...
        sig_preprocess_seiz(QRSindxs_240min),'*');
    h4 = line([time(window_indexes_seiz(:)); time(window_indexes_seiz(:))], ...
        repmat(y,numel(window_indexes_seiz),1)','Color','m');
    h10 = plot(gap_time, gap_logical*(2/3)*max(sig_preprocess_seiz), 'c');
    
    legend_vec = [h1, h2, h3, h10, h4(1)];
    legend_vec_names = {'Signal preprocessed', 'Noise detection', 'R peaks', ...
        'Gap identification', '5 min windows'};
    axis tight
    legend(legend_vec,legend_vec_names)
    xlabel('Time (hour)')
    ylabel(regexprep(file_name,'_',' '))
    savefig(fullfile(cd, 'PrepareFeatureDataset2Clustering', ...
        'FiguresAnalysedECGSignals', [file_name '_signal']))
    close all
end

% % set the top axis to go from 0:240 min: *********************************
% ax = gca;
% ax1_pos = ax.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% set(ax2,'xlim',[0 240]);
% set(ax2,'XTick', 0:40:240)
% set(ax2,'XTickLabel', num2cell(0:40:240))
% set(ax2, 'FontSize', 11)
% c = lines(7);
% set(ax2, 'XColor', c(2,:))
% xlabel(ax2,'Time (min)')
% % *************************************************************************



%% Compute the windows data for feature selection

win_samp = 5*60*fs; % 5 min window
overlap_percent = 98.333;
overlap_samp = round(win_samp/100*overlap_percent);

% obtain a step of 5 secs in order to fuse the information with EEG:
% step_samp = 5*fs;
% step_samp2 = win_samp-overlap_samp;




n_points = numel(sig_preprocess_seiz);
ind = 1;
count_wins = 0;

for ii = 1:win_samp-overlap_samp:n_points
    
    ind_end = ii+win_samp-1;
    
    if ind_end>n_points
        break;
    end
    count_wins = count_wins+1;
    
%     if count_wins==2102
%         pause
%     end

    resultsOutput(ind).sig_preprocess = sig_preprocess_seiz(ii:ind_end);
    percent_detected_noise = sum(clean_sig(ii:ind_end)==1)/length(clean_sig(ii:ind_end))*100;
    resultsOutput(ind).percent_detected_noise = percent_detected_noise;
    resultsOutput(ind).window_indexes = [ii, ind_end];
    % resultsOutput(ind).detected = clean_sig(ii:ind_end);
    
    
    %% COMPUTE THE clean_segment_indexes in each 5-min overlapped window
    if percent_detected_noise<100
        ind_clean_start = ii;
        ind_clean_end = ind_end;
        
        ind_check1 = ind_clean_end>clean_segment_indexes_seiz(:,2);
        ind_check2 = ind_clean_end<clean_segment_indexes_seiz(:,1);
        if all(ind_check1+ind_check2)
            ind_clean_end = clean_segment_indexes_seiz(find(ind_check2,1,'first')-1,2);
            % when there's noise at the end of the signal
            if ind_end>clean_segment_indexes_seiz(end,2)
                % the window will have lower samples that 7800, and should
                % therefore be the last window
                ind_clean_end = clean_segment_indexes_seiz(end,2);
            end
        end
        
        ind_check3 = ind_clean_start>clean_segment_indexes_seiz(:,2);
        ind_check4 = ind_clean_start<clean_segment_indexes_seiz(:,1);
        if all(ind_check3+ind_check4)
            ind_clean_start = clean_segment_indexes_seiz(find(ind_check4,1,'first'),1);
        end
        
        %     ver = [ind_check1, ind_check2, ind_check3, ind_check4];
        
        s = ind_check1+ind_check2+ind_check3+ind_check4;
        ind_s = find(s==1);
        if numel(ind_s)>1
            if ind_s(2)~=ind_s(1)+1
                
                ind_clean_start = [ind_clean_start; clean_segment_indexes_seiz(ind_s(1)+2:ind_s(2),1)];
                ind_clean_end = [clean_segment_indexes_seiz(ind_s(1)+1:ind_s(2)-1,2); ind_clean_end];
                
            end
        end
        
        check_equal_wins = ismember(clean_segment_indexes_seiz(:,1),ii)+ ...
            ismember(clean_segment_indexes_seiz(:,2),ind_end);
        if numel(find(check_equal_wins))>1
            % means that the window defined by [ind_clean1, ind_clean2] is
            % exactly the same as one of the 5 min windows obtained without
            % overlap and also that it contains more than one segment
            ind_clean_start = clean_segment_indexes_seiz(find(check_equal_wins),1);
            ind_clean_end = clean_segment_indexes_seiz(find(check_equal_wins),2);
        end
        
        checkClean = numel([ind_clean_start, ind_clean_end]);
        if checkClean<2
            pause
        end
        
        % check for intervals in the middle:
        % encontrar valores de clean_segments no meio dos índices que encontrei
        % porque até aqui só está a introduzir SEMPRE dois clean_segments mesmo
        % que exista mais do que um
        
        
        middle_clean_segments_ind1 = min([find(clean_segment_indexes_seiz(:,1)>=ind_clean_start(1),1,'first'), ...
            find(clean_segment_indexes_seiz(:,2)>=ind_clean_start(1),1,'first')]);
        
        middle_clean_segments_ind2 = find(clean_segment_indexes_seiz(:,2)<=ind_clean_end(end),1,'last');
        
        get_clean_segments = clean_segment_indexes_seiz(middle_clean_segments_ind1:middle_clean_segments_ind2,:);
        
        if ~isempty(get_clean_segments)
            if ind_clean_start(1)<get_clean_segments(1,2)
                get_clean_segments(1,1) = ind_clean_start(1);
            end
            
            if ind_clean_end(end)>get_clean_segments(end,2)
                get_clean_segments = [get_clean_segments; ...
                    clean_segment_indexes_seiz(middle_clean_segments_ind2+1,1) ind_clean_end(end)];
            end
            
            ind_clean_start = get_clean_segments(:,1);
            ind_clean_end = get_clean_segments(:,2);
        end
        
        % até aqui estava a seleccionar intervalos contíguos, mas com o código
        % seguinte junto esses intervalos
        diff_between_clean = ind_clean_start(2:end)-ind_clean_end(1:end-1);
        ind_single_sample_diff = find(diff_between_clean==1);
        if ~isempty(ind_single_sample_diff)
            ind_clean_start(ind_single_sample_diff+1) = [];
            ind_clean_end(ind_single_sample_diff) = [];
        end
        
        
    else
        ind_clean_start = [];
        ind_clean_end = [];
    end
    
   
    resultsOutput(ind).clean_segment_indexes = [ind_clean_start, ind_clean_end];
    
    %% COMPUTE THE QRSindxs_entire_sig_nan
    
    
    %     if size([ind_clean1, ind_clean2],1)>1
    %         pause
    %     end
    %
    ind_QRS_start = find(QRSindxs_240min>=ii, 1, 'first');
    ind_QRS_end = find(QRSindxs_240min<=ind_end, 1, 'last');
    
    QRSindxs_240min_window = QRSindxs_240min(ind_QRS_start:ind_QRS_end);
    
    % ind_abnormal_start = find(save_indAbnormalRR>=ii, 1, 'first');
    % ind_abnormal_end = find(save_indAbnormalRR<=ind_end, 1, 'last');
    % timeAbnormalRR_win = save_indAbnormalRR(ind_abnormal_start:ind_abnormal_end);
    
    % in order to avoid error when computing the features, when
    % QRSindxs_240min_window results in only one peak that vector is set 
    % to empty
    if numel(QRSindxs_240min_window)==1
        QRSindxs_240min_window = [];
    end
    
    resultsOutput(ind).QRSindxs_240min = QRSindxs_240min_window;
    % resultsOutput(ind).nAbnormalHB = numel(timeAbnormalRR_win);
    
    %% UPDATE clean_segment_indexes if gap exists
    
    time_window = time(ii:ind_end);
    
    ind_win_gap1 = find(gap_time>=time_window(1), 1, 'first');
    ind_win_gap2 = find(gap_time<=time_window(end), 1, 'last');
    gap_win_logical = gap_logical(ind_win_gap1:ind_win_gap2);
    
    
    
    % sometimes only gap_win_logical(1) = 1 --> no need to compute the
    % remaining:
    if any(gap_win_logical(2:end-1))
        
        % there is gap in the window
        
        gap_win_time = gap_time(ind_win_gap1:ind_win_gap2);
        
         % find the time indexes of the gaps existent in window:
        gap_samples_start_indexes = strfind(gap_win_logical',[0 1])';% add the zero for the cases in that the vector starts with zero
        gap_samples_end_indexes = strfind(gap_win_logical',[1 0])';
        if gap_win_logical(1)==1
            gap_samples_start_indexes = [1; gap_samples_start_indexes];
        end
        if gap_win_logical(end)==1
            gap_samples_end_indexes = [gap_samples_end_indexes; length(gap_win_logical)];
        end
        
        
        
        
        % find the indexes of start and end of gap:
        Rpeaks_time = time(QRSindxs_240min_window);
        RR_intervals_segment = Rpeaks_time(2:end)-Rpeaks_time(1:end-1);
        
        gap_samples_start_indexes(gap_samples_start_indexes==1) = [];
        indx_Rpeaks_left = [];
        
        
        % find the indexes of the R peaks near the gap
        for rr = 1:numel(gap_samples_start_indexes)
            ind_left = find(Rpeaks_time<gap_win_time(gap_samples_start_indexes(rr)), ...
                1, 'last');
            if ~isempty(ind_left)
                indx_Rpeaks_left = [indx_Rpeaks_left; ind_left];
            end
        end
        
        if any(indx_Rpeaks_left==numel(Rpeaks_time))
            indx_Rpeaks_left(indx_Rpeaks_left==numel(Rpeaks_time)) = [];
        end
        

        if ~isempty(indx_Rpeaks_left)
            other_intervals = resultsOutput(ind).clean_segment_indexes;
            for aa = 1:numel(indx_Rpeaks_left)
                
                ind_time_gap = Rpeaks_time(indx_Rpeaks_left(aa):indx_Rpeaks_left(aa)+1);

                % ind_time_gap refers to the Rpeaks that fall in the in the
                % time gap existing between files
                
                ind_gap_start = find(time>=(ind_time_gap(1)-0.35), 1,'first');
                
                ind_gap_end = find(time<=(ind_time_gap(2)+0.35), 1,'last');
                
                if ind_gap_start>ind_gap_end
                    disp('PROBLEM: not possible')
                    pause
                end
                
                interval2add = [ind_clean_start(1) ind_gap_start-1; ...
                    ind_gap_end+1 ind_clean_end(end)];
                plotFigureIntervals = 0;
                new_intervals = insertIntervalsInOtherIntervals(other_intervals, ...
                    interval2add, plotFigureIntervals);
                resultsOutput(ind).clean_segment_indexes = new_intervals;
            end
        end
        
    end
    
    
    %% PLOT FIGURE
    
    
    %     plotFigure = 1;
    if plotFigure
        % rescale sig_preprocess_seiz since different 5-min windows have
        % different amplitude:
        %         sig_preprocess_seiz = sig_preprocess_seiz/std(sig_preprocess_seiz);
        
        
        window_test = sig_preprocess_seiz(ii:ind_end);
        % time_window = time(ii:ind_end); %time(hh+round(win/2));
        figure(100)
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        h1 = plot(time, sig_preprocess_seiz); hold on
        h2 = plot(time, clean_sig*(1/3)*max(sig_preprocess_seiz));
        h3 = plot(time(QRSindxs_240min), ...
            sig_preprocess_seiz(QRSindxs_240min),'*');
        h4 = line([time(window_indexes_seiz(:)); time(window_indexes_seiz(:))], ...
            repmat(y,numel(window_indexes_seiz),1)','Color','m');
        h5 = line([time(ii) time(ind_end); time(ii) time(ind_end)], repmat(y',1,2),'Color','k');
        
        legend_vec = [h1, h2, h3, h4(1), h5(1)];
        legend_vec_names = {'Signal preprocessed', 'Noise detection', 'R peaks', ...
            '5 min windows', 'window under inspection indexes'};
        if ~isempty(ind_clean_start)
            for ll = 1:numel(ind_clean_start)
                h6 = line([time(ind_clean_start(ll)) time(ind_clean_start(ll))], y,'Color','r','LineStyle','--');
            end
            legend_vec = [legend_vec, h6];
            legend_vec_names = [legend_vec_names, {'left clean segments under inspection'}];
        end
        if ~isempty(ind_clean_end)
            for ll = 1:numel(ind_clean_end)
                h7 = line([time(ind_clean_end(ll)) time(ind_clean_end(ll))], y,'Color','g','LineStyle','--');
            end
            legend_vec = [legend_vec, h7];
            legend_vec_names = [legend_vec_names, {'right clean segments under inspection'}];
        end
        h8 = plot(time_window,window_test,'--');
        h9 = plot(time(QRSindxs_240min(ind_QRS_start:ind_QRS_end)), ...
            sig_preprocess_seiz(QRSindxs_240min(ind_QRS_start:ind_QRS_end)),'o');
        h10 = plot(linspace(time(1),time(end), numel(gap_logical)), gap_logical*(2/3)*max(sig_preprocess_seiz));
        hold off
        xlabel('Time (hour)')
        axis tight
        % xlim([time(2073600) time(end)])
        legend_vec = [legend_vec, h8, h9, h10];
        legend(legend_vec, [legend_vec_names, ...
            {'window under inspection','R peaks under inspection', ...
            'Gap identification'}])
        title(regexprep(file_name,'_',' '))
        
    end
    
%     if ind_end>=1007360%>clean_segment_indexes_seiz(end,2)
%         % it means that there are no more clean segments
%         % disp('')
%         % the window will have lower samples that 7800, and should
%         % therefore be the last window
%         % break
%         % pause
%     end
%     
    ind = ind+1;
    
end






end