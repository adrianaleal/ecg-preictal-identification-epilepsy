function hrv_out = time_features(RRI_series_segment)

% Compute time domain HRV features


% Inputs:
% - RRI_signal_segment (double): RR interval series
% - time_RR_intervals (double): irregularly time-sampled vector for each RR
%                               interval
% - plotFigure (double): flag to plot figure

% Outputs:
% - hrv_out: Table containing the following HRV metrics:
%           - NN50: number of RRIs that last more than 50 ms.
%           - pNN50: percentage of RRIs that last more than 50 ms.
%           - SDNN: standard deviation of RRIs.
%           - RMSSD: square root of the mean squared di?erences of 
%             successive RRIs.
%           - SDSD: standard deviation of the di?erences between successive 
%             RRIs.
%           - RRMean: mean of RRI series.
%           - RRMin: minimum of RRI series.
%           - RRMax: maximum of RRI series.
%           - RRVar: variance of RRI series.


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
% (5) Acharya2007: Rajendra Acharya, U., Suri, J. S., Spaan, J. A. &
% Krishnan, S. M. Advances in Cardiac Signal Processing (Springer Berlin
% Heidelberg, Berlin, Heidelberg, 2007).
% (6) Shaffer2017: Shaffer, F. & Ginsberg, J. P. An Overview of Heart Rate
% Variability Metrics and Norms. Front. Public Heal. 5, 1–17, DOI:
% 10.3389/fpubh.2017.00258 (2017).


% Warning: it is inappropriate to compare time-domain measures, especially
% those expressing overall HRV, obtained from recordings of different
% durations. Kamath2013, chapter 2, page 13


hrv_out = [];

var_names = {'NN50', 'pNN50',  'RMSSD', 'SDNN', 'SDSD', 'RRMean', ...
    'RRMin', 'RRMax', 'RRVar'};

hrv_out = array2table(zeros(1,numel(var_names)), 'VariableNames', var_names);

hrv_out.Properties.Description = 'Time Domain HRV Metrics';


%% the difference operation is meant to accentuate the high-frequency
% content of the NN interval series (Sornmo2005, chapter 8, page 570)
% it is used to compute dispersion measures that reflect short-term
% variability (Sornmo2005, chapter 8, page 570)
RR_intervals_diff_segment = diff(RRI_series_segment);


%% NN: interbeat interval of normal sinus beats (excluding ectopic
% beats) Shaffer2017

over_50ms_diff = abs(RR_intervals_diff_segment)>50;
NN50 = sum(over_50ms_diff);
col_total_power = 'NN50';
hrv_out{1,1} = NN50;
hrv_out.Properties.VariableUnits{col_total_power} = 'count';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'NN50';

%% pNN50 (ms): Percentage of successive RR intervals that differ by more
% than 50 ms

col_total_power = 'pNN50';
hrv_out{1,2} = NN50/numel(RR_intervals_diff_segment)*100;
hrv_out.Properties.VariableUnits{col_total_power} = '%';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'pNN50';

%% square root of the mean squared differences of successive NN intervals (RMSSD) 
% estimate of short-term components of HRV
% reflects the beat-to-beat variance in heart rate (Shaffer2017)
% it is used to estimate the vagally mediated changes reflected in HRV 
% (Shaffer2017)
% RMSSD method is preferred to pNN50 and NN50 because it has better
% statistical properties

col_total_power = 'RMSSD';
hrv_out{1,3} = sqrt(mean(RR_intervals_diff_segment.^2));
hrv_out.Properties.VariableUnits{col_total_power} = 'ms';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'RMSSD';


% pNN50 and RMSSD are measurements predominantly reflecting
% parasympathetic modulation of the heart
% RMSSD provides a more detailed description of short-term variability 
% whereas pNN50 is much less vulnerable than rMSSD to artifacts that may be
% present in the RR interval series (Sornmo2005, chapter 8, page 572)

%% standard deviation of differences between adjacent RR intervals (SDSD)

col_total_power = 'SDNN';
hrv_out{1,4} = std(RRI_series_segment);
hrv_out.Properties.VariableUnits{col_total_power} = 'ms';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'SDNN';

col_total_power = 'SDSD';
hrv_out{1,5} = std(RR_intervals_diff_segment);
hrv_out.Properties.VariableUnits{col_total_power} = 'ms';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'SDSD';

% it is the same as SDNN: standard deviation of the NN interval
% characterized by the following:
% depends on the length of the recording period
% estimate of overall HRV
% it is more accurate when calculated over 24h (Shaffer2017)

%% Standard deviation of the average normal-to-normal intervals (SDANN)
% estimate of the changes in the heart rate due to cycles longer than 5 min
% estimate of long-term components of HRV
% primarily reflects very slow, circadian variations in HR (Sornmo2005, 
% chapter 8, page 570)
% calculate from segments of the total monitoring period (TaskForce1996)
% SDANN = std(mean(RR_interval_signal_segment));




%% triangular index: TINN 
% too sensitive to outliers and artifacts (Shaffer2017; Clifford2006, 
% chapter 3, page 74; Acharya2007, chapter 5, page 132)
% has high correlation with the standard deviation of all RR intervals 
% (Acharya2007, chapter 5, page 132)
% the use of triangular methods is no longer suitable since they tend to 
% overestimate the variability in heart rate (Sornmo2005, chapter 8, page 
% 573) 

%%

col_total_power = 'RRMean';
hrv_out{1,6} =  mean(RRI_series_segment);
hrv_out.Properties.VariableUnits{col_total_power} = 'ms';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'RRMean';

col_total_power = 'RRMin';
hrv_out{1,7} = min(RRI_series_segment);
hrv_out.Properties.VariableUnits{col_total_power} = 'ms';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'RRMin';

col_total_power = 'RRMax';
hrv_out{1,8} = max(RRI_series_segment);
hrv_out.Properties.VariableUnits{col_total_power} = 'ms';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'RRMax';

col_total_power = 'RRVar';
hrv_out{1,9} = var(RRI_series_segment);
hrv_out.Properties.VariableUnits{col_total_power} = 'ms';
hrv_out.Properties.VariableDescriptions{col_total_power} = 'RRVar';

end