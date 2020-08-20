# Heart-Rate-Variability-Preictal-Identification-Epilepsy
Heart Rate Variability Analysis for the Identification of the Preictal Interval in Patients with Drug-resistant Epilepsy

__CLUSTERING:__
Code to perform clustering on the three-by-three combination of HRV features, using the following clustering methods:

- K-means (KM)
- Agglomerative hierarchical clustering (AH)
- Density-based spatial clustering of applications with noise (DBSCAN)
- Gaussian mixture model clustering (GMM)


__How to run the code:__

1 - Run file main_clustering.m after defining the number of features to combine (in variable feat_comb).


__CLUSTER SOLUTION EVALUATION:__
Code to perform cluster solution evaluation of the solutions obtained after clustering, using the following methods:

(1) Validity indices based on cluster labels
    - Dunn's index (DI)
    - Connectivity (C)
    - Silhouette index (SI)

(2) Validity indices based on cluster prototypes
    - Cluster separation (CS)
    - Intra cluster variance (ICV)
    - Overall cluster deviation (OD)


__How to run the code:__

1 - Run file main_cluster_evaluation.m after performing clustering and also upon the definition of the number of features to combine (in variable feat_comb).


