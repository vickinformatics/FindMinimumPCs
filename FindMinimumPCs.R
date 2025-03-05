#' FindMinimumPCs
#'
#' The `FindMinimumPCs` function determines the minimum number of principal components (PCs) to retain for downstream analysis by evaluating the explained variance in a specified dimensionality reduction method (e.g., PCA).
#' The function uses two criteria: 
#' (1) the first component that explains more than 90% of the total variance, and 
#' (2) the largest drop in explained variance between consecutive components greater than 0.1.
#' The function returns the smaller of these two values to identify the optimal number of PCs.
#'
#' @param seurat A Seurat object containing the results of a dimensionality reduction (e.g., PCA).
#' @param reduction_type Specifies the type of reduction to use (default is "pca").
#'
#' @return An integer representing the minimum number of principal components (PCs) to retain based on the explained variance.
#' 
#' @citations
#' Ianevski A, Giri AK, Aittokallio T (2022). "Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data." *Nature Communications*, 13, 1246. https://doi.org/10.1038/s41467-022-28803-w.
#' 
#' @author Vicki Do
#' @lastUpdated 2025-2-7

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