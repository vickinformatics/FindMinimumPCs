# FindMinimumPCs

## Description
This function determines the minimum number of significant principal components (PCs) required for a given dimensionality reduction (default: PCA) in a Seurat object. It calculates the minimum PCs by considering two criteria:

1. Cumulative Variance: The first component where the cumulative variance exceeds 90% while maintaining a variance contribution of less than 5% per component.
2. Variance Drop: The first component where the explained variance drops significantly (greater than 0.1) compared to the previous component.

## Note
The FindMinimumPCs function is a prerequisite for the SeuratWorkflow function. Ensure that FindMinimumPCs is loaded into your R session to allow SeuratWorkflow to function properly.
