# Heart-Rate-Variability-Preictal-Identification-Epilepsy
Heart Rate Variability Analysis for the Identiﬁcation of the Preictal Interval in Patients with Drug-resistant Epilepsy

Code to perform a feature redundancy study using the following two measures:

- Pearson's correlation coefficient
- Average mutual information using Thomas et al. Matlab code (https://www.sciencedirect.com/science/article/pii/S0022249614000571)

How to run the code:

1 - Run file main_get_feature_dataset.m to obtain the feature dataset.

2 - Define the dimension through which to evaluate redundancy in file main_feature_redundancy.m. There are two options: chosenDim = 'windows' or chosenDim = 'seizures'.

3 - Run file main_feature_redundancy.m.

