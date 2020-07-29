# Heart-Rate-Variability-Preictal-Identification-Epilepsy
Heart Rate Variability Analysis for the Identiﬁcation of the Preictal Interval in Patients with Drug-resistant Epilepsy


1 - Open a folder e.g., mainFolderRunClustering that will contain the following folders and files:
	- PrepareFeatureDataset2Clustering
	- patients_info_final_preictal_study_240min_before_seiz.mat


2 - Functions to include in PrepareFeatureDataset2Clustering folder to run both clustering_by_seizure.m and clustering_solution_evaluation_by_seizure.m:
	- parsave.m
	- parload.m
	- initializeClusteringMethods.m
	- getClusterCentroids.m
	- getMetaDataPreIctal.m
	- prepareFeatureOverlapDataset2analyse.m

3 - If computeFromScratch = 1
	- Copy ResultsEPILEPSIAE folder containing ExtractedFeaturesNoEctopicBeatsOverlap folder
	- Add the path of ResultsEPILEPSIAE to outer_folder_path variable in prepareFeatureOverlapDataset2analyse.m
	- editInvalidWindows.m

4 - If computeFromScratch = 0
	- Copy ResultsOverlap folder containing 
feature_dataset_240min_before_seizure_all_feat.mat and structure2selectFeatPreIctalStudy.mat

5 - To run clustering_by_seizure.m it is necessary to copy to folder FunctionsClustering:
	- clustering_by_seizure.m
	- DBSCAN
	- getDBSCANClustering.m
	- GaussianMixtureClustering.m

6 - To run clustering_solution_evaluation_by_seizure.m it is necessary to copy to folder FunctionsClusterSolutionEvaluation:
	- clustering_solution_evaluation_by_seizure.m
	- findContinuousClusters.m
	- clusterSeparation.m 
	- compactness.m 
	- connectivity.m 
	- dunns.m 
	- nearest_neighbours.m 
	- plotClusterSolution.m 
	- add folder ResultsPreIctalStudyClusteringOverlap containing folder 
ResultsPreICtalStudyClusteringFeatComb3ByPatientSeizure with the clustering results

7 - When running code results will be added to ResultsPreIctalStudyClusteringOverlap folder from which they should be copied.

