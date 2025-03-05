# FindMinimumPCs

## Description
The `FindMinimumPCs` function determines the minimum number of principal components (PCs) to retain for downstream analysis by evaluating the explained variance in a specified dimensionality reduction method (e.g., PCA). The function uses two criteria: 

1. The first component that explains more than 90% of the total variance, and 
2. the largest drop in explained variance between consecutive components greater than 0.1. The function returns the smaller of these two values to identify the optimal number of PCs.
