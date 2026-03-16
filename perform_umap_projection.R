#' Perform UMAP dimensionality reduction
#'
#' @description
#' Projects high-dimensional data to a lower-dimensional space using
#' Uniform Manifold Approximation and Projection (UMAP) algorithm.
#' Provides options for data preprocessing and UMAP parameter tuning.
#'
#' @param X Data frame or matrix containing input data
#' @param seed Integer, random seed for reproducibility (default: 42)
#' @param scaleX Logical, whether to scale input data (default: TRUE)
#' @param n_neighbors Integer, number of nearest neighbors used in UMAP (default: 15)
#'
#' @return A list containing:
#'   \item{Projected}{Data frame with projected coordinates}
#'   \item{UniqueData}{Original data with additional projection information}
#'   \item{umap_object}{The UMAP model object}
#'
#' @importFrom umap umap
#' @importFrom stats scale
perform_umap_projection <- function(data,
                                    seed = 42,
                                    scaleX = TRUE,
                                    n_neighbors = 15) {

  # Validate inputs
  if (!is.data.frame(data) && !is.matrix(data)) stop("Data must be a data frame or matrix.")
  if (is.null(data) || nrow(data) == 0) stop("Data must not be NULL or empty.")

  # Check for required library
  if (!requireNamespace("umap", quietly = TRUE)) {
    stop("Package 'umap' is required but not installed. Please install it.")
  }

  # Extract variables for projection (excluding Target and Label if present)
  if ("Target" %in% colnames(data) || "Label" %in% colnames(data)) {
    data_clean <- data[, !names(data) %in% c("Target", "Label")]
    labels <- if ("Label" %in% colnames(data)) data$Label else rownames(data)
    target <- if ("Target" %in% colnames(data)) data$Target else NULL
  } else {
    data_clean <- data
    labels <- rownames(data)
    target <- NULL
  }

  # Remove duplicates
  dupl <- which(duplicated(data_clean))
  if (length(dupl) > 0) {
    data_clean <- data_clean[!duplicated(data_clean),]
    if (!is.null(target)) target <- target[-dupl]
    if (!is.null(labels)) labels <- labels[-dupl]
  }

  # Scale data if requested
  if (scaleX) {
    data_clean <- scale(data_clean)
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Perform UMAP projection with specified n_neighbors
  message("Applying UMAP projection...")
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- n_neighbors
  res.umap <- umap::umap(data_clean, config = umap_config)

  # Extract projection coordinates
  projection <- res.umap$layout

  # Ensure projection has meaningful column names
  projection <- as.data.frame(projection)
  names(projection) <- paste0("Dim", seq_len(ncol(projection)))

  # Recombine with target and labels (if available)
  result <- list(Projected = projection)

  if (!is.null(target) || !is.null(labels)) {
    df_unique <- data.frame(projection)
    if (!is.null(target)) df_unique$Target <- target
    if (!is.null(labels)) df_unique$Label <- labels
    result$UniqueData <- df_unique
  }

  return(result)
}