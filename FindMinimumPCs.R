### Description: This function determines the minimum number of significant principal components (PCs) required for a given dimensionality reduction (default: PCA) in a Seurat object. It calculates the minimum PCs by considering two criteria: 1. Cumulative Variance: The first component where the cumulative variance exceeds 90% while maintaining a variance contribution of less than 5% per component. 2. Variance Drop: The first component where the explained variance drops significantly (greater than 0.1) compared to the previous component.

### Last updated: 2025-2-7

### Author: Vicki Do

FindMinimumPCs <- function(seurat, reduction_type = "pca") {
  # Check that the reduction_type exists in the Seurat object
  if (!(reduction_type %in% names(seurat@reductions))) {
    stop(paste("Reduction type '", reduction_type, "' not found in object.", sep = ""))
  }
  
  # Get the standard deviations of the specified reduction
  stdev <- seurat[[reduction_type]]@stdev
  sum_stdev <- sum(stdev)
  
  # Calculate the percentage of variance explained by each PC
  percent_stdev <- (stdev / sum_stdev) * 100
  
  # Calculate the cumulative sum of explained variance
  cumulative <- cumsum(percent_stdev)
  
  # Find the first significant component based on cumulative variance > 90%
  co1 <- which(cumulative > 90 & percent_stdev < 5)[1]
  
  # Find the largest drop in explained variance (greater than 0.1)
  co2 <- sort(which((percent_stdev[1:length(percent_stdev) - 1] - percent_stdev[2:length(percent_stdev)]) > 0.1), decreasing = TRUE)[1] + 1
  
  # Return the minimum of the two indices to determine the minimum number of significant PCs
  min_pc <- min(co1, co2)
  
  return(min_pc)
}