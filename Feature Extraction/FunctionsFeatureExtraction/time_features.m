function features_segment = time_features(RR_interval_signal_segment,...
    features_segment,ll)

% READ: 
% (1) TaskForce1996
% (2) Tripathi2013 ("Significance of Heart Rate Variability in Patients
% with Epilepsy"), from page 14
% (3) Shaffer2017 ("An Overview of Heart Rate Variability Metrics and
% Norms")
% (4) Sornmo2005 ("Bioelectrical signal processing in cardiac and 
% neurological applications")

% Warning: it is inappropriate to compare time-domain measures, especially
% those expressing overall HRV, obtained from recordings of different
% durations.

%% the difference operation is meant to accentuate the high-frequency
% content of the NN interval series (Sornmo2005, page 570)
% it is used to compute dispersion measures that reflect short-term
% variability (Sornmo2005, page 570)
RR_intervals_diff_segment = diff(RR_interval_signal_segment);


%% NN: interbeat interval of normal sinus beats (excluding ectopic
% beats) Shaffer2017
over_50ms_diff = abs(RR_intervals_diff_segment)>50;
features_segment.NN50(ll) = sum(over_50ms_diff);

%% pNN50 (ms): Percentage of successive RR intervals that differ by more
% than 50 ms
features_segment.pNN50(ll) = features_segment.NN50(ll)/numel(RR_intervals_diff_segment)*100;


%% square root of the mean squared differences of successive NN intervals (RMSSD) 
% estimate of short-term components of HRV
% reflects the beat-to-beat variance in heart rate (Shaffer2017)
% it is used to estimate the vagally mediated changes reflected in HRV (Shaffer2017)
% RMSSD method is preferred to pNN50 and NN50 because it has better
% statistical properties
RMSSD = sqrt(mean(RR_intervals_diff_segment.^2));
features_segment.RMSSD(ll) = RMSSD;

% pNN50 and RMSSD are measurements predominantly reflecting
% parasympathetic modulation of the heart
% RMSSD provides a more detailed description of short-term variability 
% whereas pNN50 is much less vulnerable than rMSSD to artifacts that may be
% present in the RR interval series (Sornmo2005)

%% standard deviation of differences between adjacent RR intervals (SDSD)


SDNN = std(RR_interval_signal_segment);
SDSD = std(RR_intervals_diff_segment);
features_segment.SDSD(ll) = SDSD;
features_segment.SDNN(ll) = SDNN;

% it is the same as SDNN: standard deviation of the NN interval
% characterized by the following:
% depends on the length of the recording period
% estimate of overall HRV
% it is more accurate when calculated over 24h (Shaffer2017)

%% Standard deviation of the average normal-to-normal intervals (SDANN)
% estimate of the changes in the heart rate due to cycles longer than 5 min
% estimate of long-term components of HRV
% primarily reflects very slow, circadian variations in HR (Sornmo2005)
% calculate from segments of the total monitoring period (TaskForce1996)
% features_segment.SDANN(ll) = std(mean(RR_interval_signal_segment));




%% triangular index: TINN 
% too sensitive to outliers and artifacts (Shaffer2017, Clifford2006,acharya2007)
% has high correlation with the standard deviation of all RR intervals (acharya2007)
% the use of triangular methods is no longer suitable since they tend to 
% overestimate the variability in heart rate (Sornmo2005) 

%%

features_segment.RRMean(ll) = mean(RR_interval_signal_segment);
features_segment.RRMin(ll) = min(RR_interval_signal_segment);
features_segment.RRMax(ll) = max(RR_interval_signal_segment);
features_segment.RRVar(ll) = var(RR_interval_signal_segment);

end