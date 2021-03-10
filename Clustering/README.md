# Heart Rate Variability Analysis for the Identification of the Preictal Interval in Patients with Drug-resistant Epilepsy

__CLUSTERING:__

Code to perform clustering on the three-by-three combination of HRV features.


__How to run the code:__

Run main_clustering.m script. This is an example code to run clustering on the HRV features extracted for a total of four seizures from patient 21902. The chosen feature combinations (see variable feat_comb_seizures) are the ones depited in Figure 2 on the article.

The following clustering methods were used: 
- K-means (KM)
- Agglomerative hierarchical clustering (AH)
- Density-based spatial clustering of applications with noise (DBSCAN). The script was developed by Yarpiz (2021). DBSCAN Clustering Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm), MATLAB Central File Exchange. Retrieved March 9, 2021.
- Gaussian mixture model clustering (GMM)


The following cluster solution evaluation methods were used:

(1) Validity indices based on cluster labels
    - Dunn's index (DI) --> dunns.m

(2) Validity indices based on cluster prototypes 
    - Intra cluster variance (ICV) --> compactness.m
    - Overall cluster deviation (OD) --> compactness.m

The script numSuplots.m is used to automatically compute the number of subplots (corresponding to the number of seizures) for each patient. The code was developed by Rob Campbell (2021). numSubplots - neatly arrange subplots (https://www.mathworks.com/matlabcentral/fileexchange/26310-numsubplots-neatly-arrange-subplots), MATLAB Central File Exchange. Retrieved March 9, 2021.


Matlab toolboxes required:

- Statistics and Machine Learning

