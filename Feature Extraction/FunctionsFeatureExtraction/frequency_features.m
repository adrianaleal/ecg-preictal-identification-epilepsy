function [hrv_fd, interp_RRI_signal_segment] = frequency_features(RRI_signal_segment,...
    time_RR_intervals, segment_size_seconds, plotFigure, varargin)

% (1) Compute the power spectrum of HRV.
% (2) Extract HRV features from that spectrum.


% READ:
% (1) TaskForce1996
% (2) Tripathi2013 ("Significance of Heart Rate Variability in Patients
% with Epilepsy"), from page 14
% (3) Clifford2006 ("Advanced Methods and Tools for ECG Data Analysis"),
% from page 74




% The spectrum must be calculated from short-term recordings of 2 (120 secs) 
% to 5 min (TaskForce1996)

hrv_fd = [];
interp_RRI_signal_segment = [];

% number of progressive beats
beats_vec = 1:length(RRI_signal_segment);

% irregularly time-sampled signal
RRI_signal_segment = RRI_signal_segment/1e3; % from ms to s

% if plotFigure
%     figure()
%     hist(RRI_signal_segment)
%     grid
%     xlabel('Sampling interval (s)')
%     ylabel('RR distribution')
% end

if length(RRI_signal_segment) < 9
    disp('The HR segment does not have a minimum of 10 beats!')
    return;
end


%% Get instantaneous Heart Rate
HR_signal_segment = 60./(RRI_signal_segment); % units: bpm

HR_signal_segment2 = 1./(RRI_signal_segment); % units: s-1 (see Sornmo2005a)

%% compute heart rate
HR = mean(HR_signal_segment); % bpm

% Zero mean to remove DC component (REF?)
% HR_signal_segment = HR_signal_segment - HR;



%% Initialize outputs

methods = {'lomb'}; % 'ar', 'welch', 'fft'


pxx_lomb  = []; calc_lomb  = false;
pxx_ar    = []; calc_ar    = false;
pxx_welch = []; calc_welch = false;
pxx_fft   = []; calc_fft   = false;

if (any(strcmp(methods, 'lomb')))
    calc_lomb = true;
end
if (any(strcmp(methods, 'ar')))
    calc_ar = true;
end
if (any(strcmp(methods, 'welch')))
    calc_welch = true;
end
if (any(strcmp(methods, 'fft')))
    calc_fft = true;
end


%% Compute Lomb-Scargle Periodogram
% estimate and plot the PSD of the nonuniformly sampled signal
% read The Lomb-Scargle Periodogram in Clifford2006, page 79

% In  reality, the highest frequency has to be somewhat lower than half the
% mean heart rate since the heart rate may fluctuate considerably so that
% the length of the longest RR interval bounds the highest frequency
% Sornmo2005a

if calc_lomb
    % 0.4 Hz is the minimum frequency resolution required to compute HF band
%     min_fn = 0.4;
%     if HR/60/2 < min_fn
%         disp('There is no frequency resolution to conduct Lomb-Scargle spectrum analysis')
%         return;
%     end

    %     [pxx_lomb,f_lomb] = plomb(HR_signal_segment,time_RR_intervals);
    [pxx_lomb,f_lomb] = plomb(RRI_signal_segment,beats_vec);
end

%% Interpolate instantaneous HR to obtain an uniformly spaced time-series
% Interpolation of the instaneous heart rate
% interpolation adds extra (erroneous) information into the time series
% and pads the FFT (in the time domain), tricking the user into assuming
% that there is a signal there, when really, there are simply not enough
% samples within a given range to allow the detection of a signal (in a
% statistically significant sense). Clifford2006a


    
new_fs = 4; % in Hz
% fn = new_fs/2;
% resampling at 4 Hz to capture oscillations up to 2 Hz according to
% the Nyquist theorem (Ernst2017,Clifford2006a,Kuusela2012,)

% highest relevant frequency (see Kuusela2012, page 14)
% fc = HR/60/2;

% instead of fc, the average Nyquist frequency can also be given by the
% following ratio: (see Clifford2006, page 79)
% fc2 = (numel(beats_vec)+1)/(2*segment_size_seconds); % =fn

%     if length(RRI_signal_segment) < (new_fs * 2) % ????
%         disp('There is no frequency resolution to conduct spectrum analysis')
%         return;
%     end

% Cubic spline interpolation
% new_fs = fc*2;
new_ts = 1 / new_fs;
interp_time_Rwave = ceil(time_RR_intervals(1)):new_ts:floor(time_RR_intervals(end));
interpMethod = 'spline';
if plotFigure
    interp_HR_signal_segment = interp1(time_RR_intervals, HR_signal_segment, ...
        interp_time_Rwave, interpMethod);
end

interp_RRI_signal_segment = interp1(time_RR_intervals, RRI_signal_segment, ...
    interp_time_Rwave, interpMethod);
    


% Build a frequency axis
% f_max = hf_band(2);
% f_axis = (new_samp_freq : new_samp_freq : f_max)';


n_samples = length(interp_RRI_signal_segment);
nfft = [];
%% Parametric AR Method
% Autoregressive power spectral density estimate — Yule-Walker method
% read acharya2007, page 136
% Advantages: it doesn't require to have analyzed data series to be in
% steady state. Thus any HRV data can be analyzed and fair HRV information
% still derived. Such analysis can be also performed at relatively shorter
% time intervals (less than 5 minutes) without missing meaningful HRV
% information. Finally this method is sensitive to rapid changes in HR
% properly showing tiny changes in autonomic balance.
% Drawback: necessity to perform massive calculations to find best order of
% autoregression model.(acharya2007, page 135)



if calc_ar
    % AR periodogram
    ar_order = 16;% see Acharya2006
    
    if length(interp_RRI_signal_segment)<ar_order
        ind_ar = cellfun('isempty', strfind(methods,'ar'));
        methods(~ind_ar) = [];
    else
        %         [pxx_ar, f_ar] = pyulear(interp_HR_signal_segment, ar_order, nfft, new_samp_freq);
        %         [pxx_ar, f_ar] = pburg(interp_HR_signal_segment, ar_order, nfft, new_samp_freq);
        %         [pxx_ar, f_ar] = pburg(interp_HR_signal_segment, ar_order);
        [pxx_ar, f_ar] = pburg(interp_RRI_signal_segment, ar_order, nfft, new_fs);
    end
    
    %     [pxx,f] = pyulear(x,order,nfft,fs)
end

%% Non-parametric methods (based upon the FFT algorithm)

%% Welch Method
if calc_welch
    welch_overlap = 50;
    %  Hamming window
%     win_length = 2^6;
    welch_overlap_samples = floor(n_samples * welch_overlap / 100);
    [pxx_welch, f_welch] = pwelch(interp_RRI_signal_segment, [],...
        [], nfft, new_fs);
    
    %     [pxx,f] = pwelch(x,window,noverlap,nfft,fs)
end


%% Non-parametric FFT periodogram
if calc_fft
    win_func = hamming(n_samples);
    [pxx_fft, f_fft] = periodogram(interp_RRI_signal_segment, win_func,...
        nfft, new_fs);
    % [pxx,f] = periodogram(x,window,nfft,fs)
end

%% FFT
% - Analyse a minimum of 5 min ECG data
% FFT assumes that time series represents a steady-state process. Because
% of that all data recordings should be conducted at highly stable
% standardized conditions, when no other factors other than current auto-
% nomic tone contributes in HRV. One of the most serious disadvantages is,
% its insensitivity to rapid transitory processes, which often possess very
% valuable information about how physiology or certain pathological
% processes behave.(acharya2007, page 135)


%%

hrv_fd = array2table(zeros(length(methods),11),'VariableNames',{'TOTAL_POWER',...
    'VLF_POWER', 'LF_POWER', 'HF_POWER', 'VLF_NORM', 'LF_NORM', 'HF_NORM',...
    'LF_TO_HF', 'LF_PEAK', 'HF_PEAK', 'BETA'},'RowNames',lower(methods));

hrv_fd.Properties.Description = 'Frequency Domain HRV Metrics';
% Method to compute power according to epilab
% TP = sum(pxx_lomb)
% VLF = sum(pxx_lomb(f_lomb<=vlf_band(2)))/TP
% LF = sum(pxx_lomb(f_lomb>lf_band(1)&f_lomb<=lf_band(2)))/TP
% HF = sum(pxx_lomb(f_lomb>hf_band(1)&f_lomb<=hf_band(2)))/TP



% Loop over power methods and calculate metrics based on each one. Loop in reverse order so that
% the first power method is the last and it's variables retain their values (pxx, column names).

vlf_band = [0.003, 0.04];
lf_band = [0.04, 0.15];
hf_band = [0.15, 0.4];
beta_band = [0, 0.1]; % see Clifford2006 plot in page 75

% Read about ULF(requires a recording period of at least 24h), VLF, LF, HF
% in Shaffer2017, Clifford2006(page 75), Acharya2007(page 134)

for ii = 1:length(methods)
    % Current PSD for metrics calculations
    pxx = eval(['pxx_' lower(methods{ii})]);
    f_axis = eval(['f_' lower(methods{ii})]);
    
    % Absolute power in each band
    % Get entire frequency range
    total_band = [f_axis(1), f_axis(end)];
    
    col_total_power = 'TOTAL_POWER';
    hrv_fd{ii,1} = freqband_power(pxx, f_axis, total_band) * 1e6;
    hrv_fd.Properties.VariableUnits{col_total_power} = 'ms^2';
    hrv_fd.Properties.VariableDescriptions{col_total_power} = sprintf('Total power (%s)', methods{ii});
    
    col_vlf_power = 'VLF_POWER';
    hrv_fd{ii,2} = freqband_power(pxx, f_axis, vlf_band) * 1e6;
    hrv_fd.Properties.VariableUnits{col_vlf_power} = 'ms^2';
    hrv_fd.Properties.VariableDescriptions{col_vlf_power} = sprintf('Power in VLF band (%s)', methods{ii});
    % This feature comprises a non-harmonic component which is affected by
    % algorithms of baseline or trend removal. This is a reliable measure
    % when extracted from a minimum of 5 min recordings. TaskForce1996
    
    col_lf_power = 'LF_POWER';
    hrv_fd{ii,3} = freqband_power(pxx, f_axis, lf_band) * 1e6;
    hrv_fd.Properties.VariableUnits{col_lf_power} = 'ms^2';
    hrv_fd.Properties.VariableDescriptions{col_lf_power} = sprintf('Power in LF band (%s)', methods{ii});
    
    col_hf_power = 'HF_POWER';
    hrv_fd{ii,4} = freqband_power(pxx, f_axis, [hf_band(1) f_axis(end)]) * 1e6;
    hrv_fd.Properties.VariableUnits{col_hf_power} = 'ms^2';
    hrv_fd.Properties.VariableDescriptions{col_hf_power} = sprintf('Power in HF band (%s)', methods{ii});
    
    % Calculate normalized power in each band (normalize by TOTAL_POWER)
    total_power = hrv_fd{ii,col_total_power};
    
    col_vlf_norm = 'VLF_NORM';
    hrv_fd{ii,5} = 100 * hrv_fd{ii,col_vlf_power} / total_power;
    hrv_fd.Properties.VariableUnits{col_vlf_norm} = 'n.u.';
    hrv_fd.Properties.VariableDescriptions{col_vlf_norm} = sprintf('VLF to total power ratio (%s)', methods{ii});
    
    col_lf_norm = 'LF_NORM';
    %     hrv_fd{ii,6} = 100 * hrv_fd{ii,col_lf_power} / total_power;
    hrv_fd{ii,6} = 100 * (hrv_fd{ii,col_lf_power}-hrv_fd{ii,col_vlf_power}) / total_power;% TaskForce1996
    hrv_fd.Properties.VariableUnits{col_lf_norm} = 'n.u.';
    hrv_fd.Properties.VariableDescriptions{col_lf_norm} = sprintf('LF to total power ratio (%s)', methods{ii});
    
    col_hf_norm = 'HF_NORM';
    %     hrv_fd{ii,7} = 100 * hrv_fd{ii,col_hf_power} / total_power;
    hrv_fd{ii,7} = 100 * (hrv_fd{ii,col_hf_power}-hrv_fd{ii,col_vlf_power}) / total_power;% TaskForce1996
    hrv_fd.Properties.VariableUnits{col_hf_norm} = 'n.u.';
    hrv_fd.Properties.VariableDescriptions{col_hf_norm} = sprintf('HF to total power ratio (%s)', methods{ii});
    
    % Calculate LF/HF ratio
    col_lf_to_hf = 'LF_TO_HF';
    hrv_fd{ii,8} = hrv_fd{ii,col_lf_power}  / hrv_fd{ii,col_hf_power};
    hrv_fd.Properties.VariableUnits{col_lf_to_hf} = 'n.u.';
    hrv_fd.Properties.VariableDescriptions{col_lf_to_hf} = sprintf('LF to HF power ratio (%s)', methods{ii});
    
    % Find peaks in the spectrum
    num_peaks = 5;
    lf_band_idx = f_axis >= lf_band(1) & f_axis <= lf_band(2);
    hf_band_idx = f_axis >= hf_band(1) & f_axis <= hf_band(2);
    [~, f_peaks_lf, ~, p_peaks_lf] = findpeaks(pxx(lf_band_idx), f_axis(lf_band_idx));
    [~, f_peaks_hf, ~, p_peaks_hf] = findpeaks(pxx(hf_band_idx), f_axis(hf_band_idx));
    if isempty(f_peaks_lf) % Assign the peak to one of the extreme values in frequency band
        f_peaks_lf(1) = max(f_axis(lf_band_idx));
    else % Sort by peak prominence & pad to length num_peaks
        [~,sort_idx] = sort(p_peaks_lf, 'descend');
        f_peaks_lf = f_peaks_lf(sort_idx)';
        f_peaks_lf = padarray(f_peaks_lf, [0,max([0,num_peaks-length(f_peaks_lf)])],NaN,'post');
        f_peaks_lf = f_peaks_lf(1:num_peaks);
    end
    if isempty(f_peaks_hf) % Assign the peak to one of the extreme values in frequency band
        f_peaks_hf(1) = max(f_axis(hf_band_idx));
    else % Sort by peak prominence & & pad to length num_peaks
        [~,sort_idx] = sort(p_peaks_hf, 'descend');
        f_peaks_hf = f_peaks_hf(sort_idx)';
        f_peaks_hf = padarray(f_peaks_hf, [0,max([0, num_peaks-length(f_peaks_hf)])],NaN,'post');
        f_peaks_hf = f_peaks_hf(1:num_peaks);
    end
    
    col_lf_peak = 'LF_PEAK';
    hrv_fd{ii,9} = f_peaks_lf(1);
    hrv_fd.Properties.VariableUnits{col_lf_peak} = 'Hz';
    hrv_fd.Properties.VariableDescriptions{col_lf_peak} = sprintf('LF peak frequency (%s)', methods{ii});
    
    col_hf_peak = 'HF_PEAK';
    hrv_fd{ii,10} = f_peaks_hf(1);
    hrv_fd.Properties.VariableUnits{col_hf_peak} = 'Hz';
    hrv_fd.Properties.VariableDescriptions{col_hf_peak} = sprintf('HF peak frequency (%s)', methods{ii});
    
    % Calculate beta (VLF slope)
    % Take the log of the spectrum in the VLF frequency band
    
    % Default beta_band to vlf_band if unspecified
    
    beta_idx = f_axis >= beta_band(1) & f_axis <= beta_band(2);
    pxx_beta_log = log10(pxx(beta_idx));
    f_axis_beta_log = log10(f_axis(beta_idx));
    
    % Fit a line and get the slope
    degree = 1;
    if numel(f_axis_beta_log)<2
        pxx_fit_beta = NaN;
    else
        pxx_fit_beta = polyfit(f_axis_beta_log, pxx_beta_log, degree);
    end
    
    col_beta = 'BETA';
    hrv_fd{ii,11} = pxx_fit_beta(1);
    hrv_fd.Properties.VariableUnits{col_beta} = 'n.u.';
    hrv_fd.Properties.VariableDescriptions{col_beta} = 'Log-log slope of frequency spectrum in the VLF band after linear regression';
end



%%
if plotFigure
    
    normalize = 0;
    detailed_legend = 1;
    plot_peaks = 1;
    yscale = 'sedhgsrh';
    
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
    subplot(2,1,1)
%     subplot(2,2,1:2)
%     string2eval = '[ax,hlines] = plotyyy(';
    legend_entries = {};
    max_vec = [];
    for ii = length(methods):-1:1
        
        pxx = eval(['pxx_' lower(methods{ii})]);
        f_axis = eval(['f_' lower(methods{ii})]);
        
        % Get PSD and normalize if requested
        if normalize
            total_power = freqband_power(pxx, f_axis, [min(f_axis), max(f_axis)]);
            pxx = pxx ./ total_power;
            pxx = pxx ./ max(pxx);
        end
        
%         if strcmpi(yscale, 'log')
%             pxx = pow2db(pxx);
%         end
%         figure()
        plot(f_axis, pxx, 'LineWidth', 1.5, 'Tag', 'freq'); hold on
        
        max_vec = [max_vec, max(pxx)];
%         string2eval = [string2eval, 'f_' lower(methods{ii}),',pxx_' lower(methods{ii}) ','];
        
        % Create legend label
        if detailed_legend
            switch lower(methods{ii})
                case 'lomb'
                    legend_entries{end+1} = 'Lomb Scargle';
                case 'fft'
                    legend_entries{end+1} = sprintf('%s ', upper(methods{ii}));
                case 'welch'
                    legend_entries{end+1} = sprintf('Welch (%d\\%% ovl.)', welch_overlap);
                case 'ar'
                    legend_entries{end+1} = sprintf('%s (order $=$ %d)', upper(methods{ii}), ar_order);
            end
        else
            legend_entries{end+1} = methods{ii};
        end
        
        
    end
    
%     string2eval = [string2eval, 'legend_entries'')'];
%     eval(string2eval)
    %% Peaks
    
%     % f_axis = eval(['f_' lower(methods{1})]);
%     % pxx = eval(['pxx_' lower(methods{1})]);
%     
%     lf_peak = hrv_fd{1,col_lf_peak};
%     
%     if plot_peaks && ~isnan(lf_peak)
%         plot(lf_peak, pxx(f_axis==lf_peak), 'bv', 'MarkerSize', 8, 'MarkerFaceColor', 'blue');
%         legend_entries{end+1} = sprintf('%.3f Hz', lf_peak);
%     end
%     hf_peak = hrv_fd{1,col_hf_peak};
%     if plot_peaks && ~isnan(hf_peak)
%         plot(hf_peak, pxx(f_axis==hf_peak), 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
%         legend_entries{end+1} = sprintf('%.3f Hz', hf_peak);
%     end
    
   
    
    %% Vertical area of frequency ranges 
    
    yrange = ylim;
    color = [0, 0, 1];
    vlf_h = fill([vlf_band(1) vlf_band(1) vlf_band(2) vlf_band(2)],...
        [yrange(1) yrange(2) yrange(2) yrange(1)],color);
    vlf_h.FaceAlpha = 0.3;
    vlf_h.EdgeAlpha = 0.3;
    legend_entries = [legend_entries,{'VLF'}];
    
    color = [1, 0, 0];
    lf_h = fill([lf_band(1) lf_band(1) lf_band(2) lf_band(2)],...
        [yrange(1) yrange(2) yrange(2) yrange(1)], color);
    lf_h.FaceAlpha = 0.3;
    lf_h.EdgeAlpha = 0.3;
    legend_entries = [legend_entries,{'LF'}];
    
    color = [1, 1, 0];
    hf_h = fill([hf_band(1) hf_band(1) hf_band(2) hf_band(2)],...
        [yrange(1) yrange(2) yrange(2) yrange(1)], color);
    hf_h.FaceAlpha = 0.3;
    hf_h.EdgeAlpha = 0.3;
    legend_entries = [legend_entries,{'HF'}];
 
    axis tight
    xlabel('Frequency (Hz)')
    hold off
    
    %% Legend
    legalpha(legend_entries)

    %% Labels
    
    % Y
    if normalize
        ylabel_prefix = 'Normalized';
        ylabel_units = 'n.u.';
    else
        ylabel_prefix = '';
        ylabel_units = 's$^2$/Hz';
        % When a signal is defined in terms only of a voltage, for instance,
        % there is no unique power associated with the stated amplitude. In
        % this case "power" is simply reckoned in terms of the square of the
        % signal, as this would always be proportional to the actual power
        % delivered by that signal into a given impedance. So one might use
        % units of V2 Hz?1 for the PSD and even though no actual "power" or
        % "energy" is specified. Source: wikipedia
    end
    if strcmpi(yscale, 'log')
        ylabel_prefix = ['Log ' ylabel_prefix];
        ylabel_units = 'dB/Hz';
    end
    
    % ylabel('Power/frequency (dBW/Hz)')
    
    ylabel(sprintf('%s Power Spectral Density (%s)', ylabel_prefix, ylabel_units));
    
    %%
    %         plomb(HR_signal_segment,time_Rwave_segment);
    %         title('Lomb-Scargle Nonuniform Sampling')
    %         ylabel('Power (dBW)')
    %         xlabel('Frequency (Hz)')
    
    subplot(212)
%     subplot(223)
    % plot(time_RR_intervals,RR_interval_signal_segment,'*-')
    plot(beats_vec,RRI_signal_segment*1e3,'*-')%, hold on
%     plot(interp_beats_vec,interp_RR_signal_segment,'*-r'), hold off
%     legend('Interpolated RRIs','Original RRIs')
    ylabel('$RRI_k = R_{i}-R_{i-1} \, (ms)$')
    xlabel('Beat \#')
    % xlabel('Time (s)')
    % legend(['Number of beats = ' num2str(length(RR_interval_signal_segment)+1)])
    axis tight
    
%     subplot(224)
%     yyaxis left
%     plot(interp_time_Rwave,interp_HR_signal_segment,'*-'), hold on
%     plot(time_RR_intervals,HR_signal_segment,'o-'), hold off
%     ylabel('Instantaneous HR (bpm)')
%     
%     yyaxis right
%     plot(interp_time_Rwave,interp_RRI_signal_segment,'*-'), hold on
%     plot(time_RR_intervals,RRI_signal_segment,'o-'), hold off
%     xlabel('Time (s)')
%     ylabel('$RRI_k$ (s)')
%     axis tight
%     
%     legend('Interpolated (cubic spline) HR signal ($f_s = 4 Hz$)',...
%            'Original HR signal',...
%            'Interpolated (cubic spline) RRI signal ($f_s = 4 Hz$)',...
%            'Original RRI signal')
% %     ylim([0 min(max_vec)])
    
end



end