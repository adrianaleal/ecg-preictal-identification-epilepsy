# Heart-Rate-Variability-Preictal-Identification-Epilepsy
Heart Rate Variability Analysis for the Identiﬁcation of the Preictal Interval in Patients with Drug-resistant Epilepsy

Functions to include in folder to run main_extract_features_overlap.m:

	- HRVFeatures folder containing:
		- CRPtool folder from epi2pub_with_changes\Toolbox
		- approximate_entropy.m (if approximateEntropy.m from Matlab is not used)
		- cerecurr_y.m, recurrence_quantification_analysis.m and tdrecurr_y.m from HRVFeatures\RecurrenceQuantificationAnalysis\RecurrencePlot_ToolBox
		- freqband_power.m from HRVFeatures\rhrv-master\src\util
		- detrended_fluctuation_analysis.m
		- frequency_features.m
		- getAverageMutualInformation.m
		- HRV_features.m
		- lyap_compu.m
		- poincare_plot.m
		- recurrenceAnalysis.m
		- remove_ect.m
		- sample_entropy.m
		- time_features.m
		- HRV_features.m
		- RR_intervals.m
	- ResultsEPILEPSIAE folder containing NoiseDetection folder (with all patients that appear in patients_info_final_preictal_study_240min_before_seiz.mat. otherwise it is not computing correctly)
	- main_extract_features_overlap.m	
	- GetAllFilesInDir4FeatureExtraction.m 
	- getSignalWithOverlapSlidingWindow.m	
	- getRRInterval.m
	- insertIntervalsInOtherIntervals.m
	- getHRVfeatures.m
	- patients_info_final_preictal_study_240min_before_seiz.mat
	- checkFeatureExtractionFiles.m
	- parsave.m
	- parload.m

Toolboxes required to be installed:

- Signal Processing
- Image Processing
- Predictive Maintenance
