# Heart-Rate-Variability-Preictal-Identification-Epilepsy
Heart Rate Variability Analysis for the Identiﬁcation of the Preictal Interval in Patients with Drug-resistant Epilepsy

1 - Open a folder e.g., mainFolderRunFeatureRedundancy that will contain the following folders and files:
	- PrepareFeatureDataset2Clustering
	- patients_info_final_preictal_study_240min_before_seiz.mat


2 - Functions to include in PrepareFeatureDataset2Clustering folder to run feature_redundancy_assessment:
	- parsave.m
	- parload.m
	- getMetaDataPreIctal.m
	- prepareFeatureOverlapDataset.m


3 - If computeFromScratch = 1
	- Copy ResultsEPILEPSIAE folder containing ExtractedFeaturesNoEctopicBeatsOverlap folder
	- Add the path of ResultsEPILEPSIAE to outer_folder_path variable in prepareFeatureOverlapDataset2analyse.m
	- editInvalidWindows.m

4 - If computeFromScratch = 0
	- Copy ResultsOverlap folder containing 
feature_dataset_240min_before_seizure_all_feat.mat and structure2selectFeatPreIctalStudy.mat


5 - To run feature_redundancy_assessment.m it is necessary to copy to folder FunctionsFeatureSelection:
	
	- amiGrinstead.m 
	- average_mutual_information.m (from MATLAB folder)
	- corr_features.m (from MATLAB folder)
	- feature_redundancy_assessment.m
	- MutualInformation.m (from MATLAB folder)
	- Entropy (from MATLAB folder)

6 - When running code results will be added to ResultsOverlap from which they should be copied.

Toolboxes required to be installed:

- Predictive Maintenance

