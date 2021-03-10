# Heart Rate Variability Analysis for the Identification of the Preictal Interval in Patients with Drug-resistant Epilepsy

__FEATURE EXTRACTION:__

Code to extract HRV features from the irregularly sampled RR interval series.


__How to run the code:__

Run main_extract_features.m script. This is an example code to run feature extraction using data from two 5-min RR interval series (corresponding to the ones depicted in Figures S5, S6 and S7 in Supplementary Material, extracted from the first
seizure of patient 8902). 


Functions included in folder HRVFeatures to run main_extract_features.m:
		
	- freqband_power.m
	- frequency_features.m
	- monofractal_detrended_fluctuation_analysis.m
	- poincare_plot.m
	- recurrence_quantification_analysis.m
	- recurrenceAnalysis.m
	- sample_entropy.m
	- time_features.m

The script subplot_tight.m is used to automatically plot the last figure in code, with minimum margins. The code was developed by Nikolay S. (2021). Controllable tight subplot (https://www.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot), MATLAB Central File Exchange. Retrieved March 9, 2021.

The script legalpha.m was used to produce the frequency figure. The code was developed by Chad Greene (2021). legalpha (https://www.mathworks.com/matlabcentral/fileexchange/47550-legalpha), MATLAB Central File Exchange. Retrieved March 9, 2021.


Matlab toolboxes required:

- Signal Processing
- Statistics and Machine Learning
- Predictive Maintenance


