function [] = plotSignal2Features(feat_names, patients_name, ...
    seizure_names, features_folder_path, patient_info)



% Plot the original ECG time series with the R-peak identification, along
% with the RR interval series and the features obtained

% Inputs:
% - feat_names (cell): names of the features to analyse
% - patients_name (cell): patients index name
% - seizure_names (cell): seizures index name
% - features_folder_path (char):
% - patient_info (struct): structure containing information for each
%   patient

close all
figureFolder2save = fullfile(cd, 'DataFiguresECG2Features');


% CHECK THE FOLLOWING: ****************************************************

plotNonOverlappedWindows = 0;
plotSpecificWindows = 0;
plotNormalizedFeatures = 1;

plotTimeMin = 1;
plotTimeMinReverse = 1;
% plotTimeMin = 0 --> plot original time vector in hours
% plotTimeMin = 1 && plotTimeMinReverse = 0 --> plot 0-240 time vector in min
% plotTimeMinReverse = 1 && plotTimeMinReverse = 1 --> plot 240-0 time vector in min

% *************************************************************************

if plotTimeMinReverse
    vec_time_min = -240:40:0;
    vec_time_min_label = 240:-40:0;
else
    vec_time_min = 0:40:240;
    vec_time_min_label = 0:40:240;
end


for pp = 37%:numel(patients_name)
    
    
    file_name = patients_name{pp}
    disp(['Patient ' file_name])
    patFolderName = ['pat_' file_name];
    seizFoldersFeatures = dir(fullfile(features_folder_path, patFolderName));
    seizFoldersFeatures = seizFoldersFeatures(~ismember({seizFoldersFeatures.name},{'.','..'}));
    
    n_seizures_pat = sum(strcmp(seizure_names, file_name));
    
    if plotSpecificWindows
        fs = patient_info(pp).samp_rate;
    end
    
    for ss = 6%:n_seizures_pat
        
        disp(['Seizure ' num2str(ss)])
        
        
        hrv_features = parload(fullfile(seizFoldersFeatures(ss).folder, ...
            seizFoldersFeatures(ss).name), ...
            [patFolderName '_seizure' num2str(ss)], 'hrv_features');
        hrv_features_original = hrv_features;
        
        clear hrv_features
        
        
        hrv_features = parload(fullfile(seizFoldersFeatures(ss).folder, ...
            seizFoldersFeatures(ss).name), ...
            [patFolderName '_seizure' num2str(ss)], 'hrv_features');
        
        time_vec_original = vertcat(hrv_features.time);
        
        % load the ECG signal, the R-peaks and the RR interval series
        load(fullfile(figureFolder2save, ['pat_' file_name '_seizure' ...
            num2str(ss), '_save2FigurePlotECG.mat']), ...
            'save2FigurePlotECG')
        
        if plotSpecificWindows
            % when I want to observe the location of certain windows in plot:
            % example: 5-min window after 80 min before seizure and 5-min
            % window after 25 min before seizure onset
            window_indexes= save2FigurePlotECG.overlapped_window_indexes;
            window_secs = window_indexes/fs;
            window_min = window_secs/60;
            win1 = window_min(1940,:)
            win2 = window_min(2569,:)
        end
        
        % get the ECG signal:
        sig_preprocess_seiz = save2FigurePlotECG.sig_preprocess_seiz;
        % get the time vector:
        time_recording = save2FigurePlotECG.time_recording;
        % get the binary vector indicating the noisy segments:
        clean_sig = save2FigurePlotECG.clean_sig;
        % get the R-peaks location:
        QRSindxs_entire_sig_seiz = save2FigurePlotECG.QRSindxs_240min;
        % get the binary vector indicating the existence of time gaps
        % between the recordings:
        gap_logical = save2FigurePlotECG.gap_logical;
        
        duration_seiz_min = (time_recording(end)-time_recording(1))/60;
        
        if plotNonOverlappedWindows
            % get the nonoverlapped windows used in noise detection:
            window_indexes_seiz = save2FigurePlotECG.window_indexes_seiz;
        end
        
        if plotTimeMin % plot 0-240 time vector in min
            time = ([0 cumsum(diff(time_recording))])/60;
            time_features = linspace(0, time(end), numel(time_vec_original));
            
            xlabel_name = 'Time (min)';
            
            if plotTimeMinReverse % plot 240-0 time vector in min
                time = time-duration_seiz_min;
                time_features = time_features-duration_seiz_min;
            end
        else
            time = time_recording;
            time_features = time_vec_original';
            xlabel_name = 'Time (hour)';
        end
        
        time_gap = linspace(time(1),time(end), numel(gap_logical));
        
        figure(100)
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        %% plot the ECG signal and the R-peaks *****************************
        sp1 = subplot_tight(3,1,1, [0.06, 0.04]);
        
        h1 = plot(time, sig_preprocess_seiz); hold on
        h2 = plot(time, clean_sig*(1/3)*max(sig_preprocess_seiz));
        h3 = plot(time(QRSindxs_entire_sig_seiz), sig_preprocess_seiz(QRSindxs_entire_sig_seiz),'*');
        h10 = plot(time_gap, gap_logical*(2/3)*max(sig_preprocess_seiz));
        legend_vec_names = {'ECG signal preprocessed', 'Noise detection', ...
            'R peak detection', 'File gap detection'};
        legend_vec = [h1, h2, h3, h10];
        
        
        if plotNonOverlappedWindows
            % plot the nonoverlapped windows used in noise detection:
            y = [min(sig_preprocess_seiz) max(sig_preprocess_seiz)];
            h4 = line([time(window_indexes_seiz(:)); time(window_indexes_seiz(:))], ...
                repmat(y,numel(window_indexes_seiz),1)','Color','m');
            legend_vec = [legend_vec, h4(1)];
            legend_vec_names = [legend_vec_names, {'5 min windows'}];
        end
        
        axis tight, grid on
        legend(legend_vec, legend_vec_names)
        ylabel('mV')
        ylim([-21 max(sig_preprocess_seiz)])
        if plotTimeMin
            set(gca,'XTick', vec_time_min)
            set(gca,'XTickLabel', num2cell(vec_time_min_label))
        end
        set(gca,'FontSize',14)
        
        title(['Patient ' file_name ', Seizure ' num2str(ss)])
        
        %% plot the RR interval series ************************************
        sp2 = subplot_tight(3,1,2, [0.06, 0.04]);
        for ll = 1:size(save2FigurePlotECG.RRI_240min)
            if plotTimeMin
                
                time_RR_hour = save2FigurePlotECG.time_Rpeaks_240min{ll};
                init_hour = time_RR_hour(1)-save2FigurePlotECG.time_Rpeaks_240min{1}(1);
                
                time_RR_secs = [init_hour init_hour+cumsum(diff(time_RR_hour))];
                time_RR = time_RR_secs/60; % in min
                time_RR_corrected_hour = save2FigurePlotECG.time_Rpeaks_240min_corrected{ll};
                time_RR_corrected_secs = [init_hour init_hour+cumsum(diff(time_RR_corrected_hour))];
                time_RR_corrected = time_RR_corrected_secs/60;% in min
                
                if plotTimeMinReverse % plot 240-0 time vector in min
                    time_RR = time_RR-duration_seiz_min;
                    time_RR_corrected = time_RR_corrected-duration_seiz_min;
                end
                
                h1 = plot(time_RR, save2FigurePlotECG.RRI_240min{ll}*1e3, 'b');
                hold on
                h2 = plot(time_RR_corrected, ...
                    save2FigurePlotECG.RRI_240min_corrected{ll}*1e3, '--r');
            else
                h1 = plot(save2FigurePlotECG.time_Rpeaks_240min{ll}, ...
                    save2FigurePlotECG.RRI_240min{ll}*1e3, 'b');
                hold on
                h2 = plot(save2FigurePlotECG.time_Rpeaks_240min_corrected{ll}, ...
                    save2FigurePlotECG.RRI_240min_corrected{ll}*1e3, '--r');
            end
        end
        
        
        axis tight, grid minor
        ylabel('ms')
        if plotTimeMin
            % if plot the lines identifying the 5-min windows from which
            % the features images were taken
            if plotSpecificWindows
                xline(win1(1), 'k--');
                xline(win1(2), 'k--');
                xline(win2(1), 'k--');
                xline(win2(2), 'k--');
            end
            
            set(gca,'XTick', vec_time_min)
            set(gca,'XTickLabel', num2cell(vec_time_min_label))
        end
        set(gca,'FontSize',14)
        legend([h1 h2], {'RR intervals', 'RR intervals corrected'})
        
        hold off
        
        
        
        %% plot the features **********************************************
        
        disp(['Features to plot: ' char(join(feat_names, ', '))])
        
        
        % load of the original data:
        feature_dataset = editInvalidWindows(feat_names', hrv_features);
        feat_data = [];
        for ff = 1:numel(feat_names)
            feat_data = [feat_data, ...
                vertcat(feature_dataset.(feat_names{ff}))];
        end
        feat_data2plot = feat_data;
        
        if plotNormalizedFeatures
            feat_data_norm = (feat_data-min(feat_data))./(max(feat_data)-min(feat_data));
            feat_data2plot = feat_data_norm;
        end
        
        time_data = repmat(time_features',1,numel(feat_names));
        
        
        
        sp3 = subplot_tight(3,1,3, [0.06, 0.04]);
        plot(time_data, feat_data2plot)
        axis tight, grid minor
        legend(replaceFeatureNames(regexprep(feat_names,'_',' ')))
        if plotTimeMin
            set(gca,'XTick', vec_time_min)
            set(gca,'XTickLabel', num2cell(vec_time_min_label))
            linkaxes([sp1 sp2 sp3],'x')
        end
        xlabel(xlabel_name)
        ylabel('n. u.')
        set(gca,'FontSize',14)
        
        
        
        % figPathName = fullfile(cd,'FiguresAnalysedECGSignals', ...
        %     [file_name '_seizure' num2str(ss) '_feats_' ...
        %     char(join(feat_names, '_')) '_ECG2Feat']);
        % savefig(figPathName)
        % tightfig
        % export_fig(gcf, [figPathName '.pdf'], '-painters', '-nocrop', ...
        %     '-transparent')
        
        
    end
    
end



end
