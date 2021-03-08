function [hrv_out, interp_RRI_signal_segment] = frequency_features(RRI_signal_segment,...
    time_RR_intervals, plotFigure)

% (1) Compute the power spectrum of HRV.
% (2) Extract HRV features from that spectrum.

% this code was adapted from function mhrv.hrv.hrv_freq.m in 
% https://github.com/physiozoo/mhrv and uses function freqband_power.m from
% this repository


% Inputs:
% - RRI_signal_segment (double): RR interval series
% - time_RR_intervals (double): irregularly time-sampled vector for each RR
%                               interval
% - plotFigure (double): flag to plot figure

% Outputs:
% - hrv_out: Table containing the following HRV metrics computed using the 
%            Lomb-Scargle Periodogram:
%           - TOTAL_POWER: Total power in all three bands combined.
%           - VLF_POWER: Power in the VLF band.
%           - LF_POWER: Power in the LF band.
%           - HF_POWER: Power in the HF band.
%           - VLF_NORM: Ratio between VLF power and total power.
%           - LF_NORM: Ratio between LF power and total power.
%           - HF_NORM: Ratio between HF power and total power.
%           - LF_TO_HF: Ratio between LF and HF power.
% - interp_RRI_signal_segment: interpolated RR interval series to compute
%   the nonlinear features


% READ:
% (1) TaskForce1996: Electrophysiology, T. F. o. t. E. S. o. C., American,
% T. N., & Electrophysiology, S. o. P. Heart Rate Variability. Circulation
% 93, 1043–1065, DOI: 10.1161/01.CIR.93.5.1043 (1996).
% (2) Kamath2013: Kamath, M.V., Watanabe, M., & Upton, A. (Eds.). (2013). 
% Heart Rate Variability (HRV) Signal Analysis: Clinical Applications (1st 
% ed.). CRC Press. DOI: doi.org/10.1201/b12756
% (3) Clifford2006: Clifford, G. D., Azuaje, F. & Mcsharry, P. E. Advance
% Methods and Tools for ECG Data Analysis (Artech House, 2006)
% (4) Sornmo2005: Sornmo, L. & Laguna, P. Bioelectrical signal processing
% in cardiac and neurological applications (Elsevier Academic Press, 2005).
% (5) Rajendra Acharya, U., Paul Joseph, K., Kannathal, N., Lim, C. M. &
% Suri, J. S. Heart rate variability: a review. Med. & Biol. Eng. & Comput.
% 44, 1031–1051, DOI: 10.1007/s11517-006-0119-0 (2006). s11517-006-0119-0.
% (6) Acharya2007: Rajendra Acharya, U., Suri, J. S., Spaan, J. A. &
% Krishnan, S. M. Advances in Cardiac Signal Processing (Springer Berlin
% Heidelberg, Berlin, Heidelberg, 2007).
% (7) Shaffer2017: Shaffer, F. & Ginsberg, J. P. An Overview of Heart Rate
% Variability Metrics and Norms. Front. Public Heal. 5, 1–17, DOI:
% 10.3389/fpubh.2017.00258 (2017).
% (8) Ernst2017: Ernst, G. Hidden Signals—The History and Methods of Heart
% Rate Variability. Front. Public Heal. 5, 265, DOI:
% 10.3389/fpubh.2017.00265 (2017).


% The spectrum must be calculated from short-term recordings of 2 to 5 min
% (TaskForce1996)

hrv_out = [];
interp_RRI_signal_segment = [];

% number of progressive beats
beats_vec = 1:length(RRI_signal_segment);

% irregularly time-sampled RR interval series
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

% HR_signal_segment = 1./(RRI_signal_segment); % units: s-1 (see
% Sornmo2005ch8)

%% Initialize methods to compute the RR interval series frequency spectrum

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
% read The Lomb-Scargle Periodogram in Clifford2006, chapter 3, page 79

% "In reality, the highest frequency has to be somewhat lower than half the
% mean heart rate since the heart rate may fluctuate considerably so that
% the length of the longest RR interval bounds the highest frequency"
% (Sornmo2005, chapter 8, page 590)

if calc_lomb
    % 0.4 Hz is the minimum frequency resolution required to compute HF band
    % min_fn = 0.4;
    % if HR/60/2 < min_fn
    %     disp('There is no frequency resolution to conduct Lomb-Scargle spectrum analysis')
    %     return;
    % end
    
    % [pxx_lomb,f_lomb] = plomb(HR_signal_segment,time_RR_intervals);
    [pxx_lomb,f_lomb] = plomb(RRI_signal_segment, beats_vec);
end

%% Interpolate instantaneous HR to obtain an uniformly spaced time-series
% Interpolation of the instantaneous heart rate
% "interpolation adds extra (erroneous) information into the time series
% and pads the FFT (in the time domain), tricking the user into assuming
% that there is a signal there, when really, there are simply not enough
% samples within a given range to allow the detection of a signal (in a
% statistically significant sense)." Clifford2006, chapter 3, page 79

% NOTE: THIS INTERPOLATION WILL BE NEEDED  TO COMPUTE THE NONLINEAR
% FEATURES
    
new_fs = 4; % in Hz
% fn = new_fs/2;
% resampling at 4 Hz to capture oscillations up to 2 Hz according to
% the Nyquist theorem (Ernst2017; Clifford2006, chapter 3; Kamath2013,
% chapter 2)

% highest relevant frequency (see Kamath2013, chapter 2, page 14)
% fc = HR/60/2;


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
    

n_samples = length(interp_RRI_signal_segment);
nfft = [];

%% Parametric AR Method
% Autoregressive power spectral density estimate — Yule-Walker method
% read Acharya2007, chapter 5, page 136
% Advantages: "it doesn't require to have analyzed data series to be in
% steady state. Thus any HRV data can be analyzed and fair HRV information
% still derived. Such analysis can be also performed at relatively shorter
% time intervals (less than 5 minutes) without missing meaningful HRV
% information. Finally this method is sensitive to rapid changes in HR
% properly showing tiny changes in autonomic balance."
% Drawback: "necessity to perform massive calculations to find best order
% of autoregression model." (Acharya2007, chapter 5, page 135)



if calc_ar
    % AR periodogram
    ar_order = 16; % see Acharya2006, chapter 5, page 136
    
    if length(interp_RRI_signal_segment)<ar_order
        ind_ar = cellfun('isempty', strfind(methods,'ar'));
        methods(~ind_ar) = [];
    else
        % [pxx_ar, f_ar] = pyulear(interp_HR_signal_segment, ar_order, nfft, new_fs);
        % [pxx_ar, f_ar] = pburg(interp_HR_signal_segment, ar_order, nfft, new_fs);
        % [pxx_ar, f_ar] = pburg(interp_HR_signal_segment, ar_order);
        [pxx_ar, f_ar] = pburg(interp_RRI_signal_segment, ar_order, nfft, new_fs);
    end
end

%% Non-parametric methods (based upon the FFT algorithm)

%% Welch Method
if calc_welch
    welch_overlap = 50;
    % Hamming window
    % win_length = 2^6;
    welch_overlap_samples = floor(n_samples * welch_overlap / 100);
    [pxx_welch, f_welch] = pwelch(interp_RRI_signal_segment, [],...
        [], nfft, new_fs);
end

%% FFT periodogram
if calc_fft
    win_func = hamming(n_samples);
    [pxx_fft, f_fft] = periodogram(interp_RRI_signal_segment, win_func,...
        nfft, new_fs);
end

%% FFT
% - Analyse a minimum of 5 min ECG data
% "FFT assumes that time series represents a steady-state process. Because
% of that all data recordings should be conducted at highly stable
% standardized conditions, when no other factors other than current auto-
% nomic tone contributes in HRV. One of the most serious disadvantages is,
% its insensitivity to rapid transitory processes, which often possess very
% valuable information about how physiology or certain pathological
% processes behave dynamically." (Acharya2007, chapter 5, page 135)

%%

var_names = {'TOTAL_POWER', 'VLF_POWER', 'LF_POWER', 'HF_POWER', ...
    'VLF_NORM', 'LF_NORM', 'HF_NORM', 'LF_TO_HF'};

hrv_out = array2table(zeros(length(methods),numel(var_names)), ...
    'VariableNames', var_names, 'RowNames', lower(methods));

hrv_out.Properties.Description = 'Frequency Domain HRV Metrics';



vlf_band = [0.003, 0.04];
lf_band = [0.04, 0.15];
hf_band = [0.15, 0.4];

% Read about ULF(requires a recording period of at least 24h), VLF, LF, HF
% in Shaffer2017; Clifford2006, chapter 3, page 75; Acharya2007, chapter 5,
% page 134

for ii = 1:length(methods)
    % Current PSD for metrics calculations
    pxx = eval(['pxx_' lower(methods{ii})]);
    f_axis = eval(['f_' lower(methods{ii})]);
    
    % Absolute power in each band
    % Get entire frequency range
    total_band = [f_axis(1), f_axis(end)];
    
    col_total_power = 'TOTAL_POWER';
    hrv_out{ii,1} = freqband_power(pxx, f_axis, total_band) * 1e6;
    hrv_out.Properties.VariableUnits{col_total_power} = 'ms^2';
    hrv_out.Properties.VariableDescriptions{col_total_power} = sprintf('Total power (%s)', methods{ii});
    
    col_vlf_power = 'VLF_POWER';
    hrv_out{ii,2} = freqband_power(pxx, f_axis, vlf_band) * 1e6;
    hrv_out.Properties.VariableUnits{col_vlf_power} = 'ms^2';
    hrv_out.Properties.VariableDescriptions{col_vlf_power} = sprintf('Power in VLF band (%s)', methods{ii});
    % This feature comprises a non-harmonic component which is affected by
    % algorithms of baseline or trend removal. This is a reliable measure
    % when extracted from a minimum of 5 min recordings. TaskForce1996
    
    col_lf_power = 'LF_POWER';
    hrv_out{ii,3} = freqband_power(pxx, f_axis, lf_band) * 1e6;
    hrv_out.Properties.VariableUnits{col_lf_power} = 'ms^2';
    hrv_out.Properties.VariableDescriptions{col_lf_power} = sprintf('Power in LF band (%s)', methods{ii});
    
    col_hf_power = 'HF_POWER';
    hrv_out{ii,4} = freqband_power(pxx, f_axis, [hf_band(1) f_axis(end)]) * 1e6;
    hrv_out.Properties.VariableUnits{col_hf_power} = 'ms^2';
    hrv_out.Properties.VariableDescriptions{col_hf_power} = sprintf('Power in HF band (%s)', methods{ii});
    
    % Calculate normalized power in each band (normalize by TOTAL_POWER)
    total_power = hrv_out{ii,col_total_power};
    
    col_vlf_norm = 'VLF_NORM';
    hrv_out{ii,5} = 100 * hrv_out{ii,col_vlf_power} / total_power;
    hrv_out.Properties.VariableUnits{col_vlf_norm} = 'n.u.';
    hrv_out.Properties.VariableDescriptions{col_vlf_norm} = sprintf('VLF to total power ratio (%s)', methods{ii});
    
    col_lf_norm = 'LF_NORM';
    % hrv_out{ii,6} = 100 * hrv_out{ii,col_lf_power} / total_power;
    hrv_out{ii,6} = 100 * (hrv_out{ii,col_lf_power}-hrv_out{ii,col_vlf_power}) / total_power;% TaskForce1996
    hrv_out.Properties.VariableUnits{col_lf_norm} = 'n.u.';
    hrv_out.Properties.VariableDescriptions{col_lf_norm} = sprintf('LF to total power ratio (%s)', methods{ii});
    
    col_hf_norm = 'HF_NORM';
    % hrv_out{ii,7} = 100 * hrv_out{ii,col_hf_power} / total_power;
    hrv_out{ii,7} = 100 * (hrv_out{ii,col_hf_power}-hrv_out{ii,col_vlf_power}) / total_power;% TaskForce1996
    hrv_out.Properties.VariableUnits{col_hf_norm} = 'n.u.';
    hrv_out.Properties.VariableDescriptions{col_hf_norm} = sprintf('HF to total power ratio (%s)', methods{ii});
    
    % Calculate LF/HF ratio
    col_lf_to_hf = 'LF_TO_HF';
    hrv_out{ii,8} = hrv_out{ii,col_lf_power}  / hrv_out{ii,col_hf_power};
    hrv_out.Properties.VariableUnits{col_lf_to_hf} = 'n.u.';
    hrv_out.Properties.VariableDescriptions{col_lf_to_hf} = sprintf('LF to HF power ratio (%s)', methods{ii});
end



%%
if plotFigure
    
    normalize = 0;
    detailed_legend = 1;
    plot_peaks = 1;
    yscale = ' ';
    
    figure()
    set(gcf,'units','normalized','outerposition',[0 0 0.5 1])
    subplot(211)
    legend_entries = {};
    for ii = length(methods):-1:1
        
        pxx = eval(['pxx_' lower(methods{ii})]);
        f_axis = eval(['f_' lower(methods{ii})]);
        
        % Get PSD and normalize if requested
        if normalize
            total_power = freqband_power(pxx, f_axis, [min(f_axis), max(f_axis)]);
            pxx = pxx ./ total_power;
            pxx = pxx ./ max(pxx);
        end
        
        if strcmpi(yscale, 'log')
            pxx = pow2db(pxx);
        end
        plot(f_axis, pxx, 'LineWidth', 1.5, 'Tag', 'freq'); hold on
        
        
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
    
    if normalize
        ylabel_prefix = 'Normalized';
        ylabel_units = 'n.u.';
    else
        ylabel_prefix = '';
        % ylabel_units = 's$^2$/Hz'; % if latex
        ylabel_units = 's^2/Hz';
    end
    if strcmpi(yscale, 'log')
        ylabel_prefix = ['Log ' ylabel_prefix];
        ylabel_units = 'dB/Hz';
    end
    
    ylabel(sprintf('%s Power Spectral Density (%s)', ylabel_prefix, ylabel_units));
    
    %%
    
    subplot(212)
    plot(beats_vec,RRI_signal_segment*1e3,'*-')
    
    % ylabel('$RRI_k = R_{i}-R_{i-1} \, (ms)$') % if latex
    ylabel('RRI_k = R_i-R_{i-1} (ms)')
    
    % xlabel('Beat \#')% if latex
    xlabel('Beat #')
    
    axis tight
    
end

end