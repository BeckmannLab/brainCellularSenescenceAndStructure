
#' Processing expression data from assay
#'
#' For raw counts, filter genes and samples, then estimate precision weights using linear mixed model weighting by number of cells observed for each sample.  For normalized data, only weight by number of cells
#'
#' @param y matrix of counts or log2 CPM
#' @param formula regression formula for differential expression analysis
#' @param data metadata used in regression formula
#' @param n.cells array of cell count for each sample
#' @param min.cells minimum number of observed cells for a sample to be included in the analysis
#' @param min.count minimum number of reads for a gene to be considered expressed in a sample.  Passed to \code{edgeR::filterByExpr}
#' @param min.samples minimum number of samples passing cutoffs for cell cluster to be retained
#' @param min.prop minimum proportion of retained samples with non-zero counts
#' @param min.total.count minimum total count required per gene for inclusion
#' @param isCounts logical, indicating if data is raw counts
#' @param normalize.method normalization method to be used by \code{calcNormFactors}
#' @param span Lowess smoothing parameter using by \code{variancePartition::voomWithDreamWeights()}
#' @param quiet show messages
#' @param weights matrix of precision weights
#' @param rescaleWeightsAfter default = FALSE, should the output weights be scaled by the input weights
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arguments passed to \code{dream}
#'
#' @return \code{EList} object storing log2 CPM and precision weights
#'
#' @seealso \code{processAssays()}
#' @importFrom BiocParallel SerialParam
#' @importClassesFrom limma EList
#' @importFrom variancePartition voomWithDreamWeights
#' @importFrom edgeR calcNormFactors filterByExpr DGEList
#' @importFrom methods is new
#' @importFrom stats model.matrix var
#' @importFrom SummarizedExperiment colData assays
#' @importFrom S4Vectors as.data.frame
#' @importFrom lme4 subbars
#' @importFrom MatrixGenerics colMeans2
#'


# y= countsMatrix; formula= form; data = data; n.cells= n.cells; min.cells = min.cells; min.count = min.count; min.samples = min.samples; min.prop = .4; min.total.count = 15; isCounts = isCounts; normalize.method = normalize.method; span = "auto"; quiet = quiet; weights = NULL; rescaleWeightsAfter = FALSE; BPPARAM = BPPARAM

library(rlang)
library(dreamlet)
library(edgeR)


processOneAssay_edgeR <- function(y, formula, data, n.cells, min.cells = 5, min.count = 2, min.samples = 4, min.prop = .4, min.total.count = 15, isCounts = TRUE, normalize.method = "TMM", span = "auto", quiet = TRUE, weights = NULL, rescaleWeightsAfter = FALSE, BPPARAM = SerialParam(), ...) {
  checkFormula(formula, data)

  if (is.null(n.cells)) {
    stop("n_cells must not be NULL")
  }
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  # samples to include of they have enough observed cells
  include <- (n.cells >= min.cells)

  # if no samples are retained
  if (sum(include) == 0) {
    return(NULL)
  }

  # subset expression and data
  y <- y[, include, drop = FALSE]
  data <- droplevels(data[include, , drop = FALSE])
  # n.cells <- n.cells[include]

  # if there are too few remaining samples
  if (nrow(data) < min.samples | nrow(y) == 0) {
    return(NULL)
  }

  if (!isCounts) {
    stop("isCounts = FALSE is not currently supported")
  }

  # Get count data and normalize
  y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
  # y <- calcNormFactors(y, method = normalize.method)

  # drop any constant terms from the formula
  formula <- removeConstantTerms(formula, data)

  # Drop variables in a redundant pair
  formula <- dropRedundantTerms(formula, data)

  # get samples with enough cells
  # filter genes
  # design: model.matrix( subbars(formula), data)
  # Design often includes batch and donor, which are very small
  #   this causes too many genes to be retained
  keep <- suppressWarnings(filterByExpr(y, 
      min.count = min.count, 
      min.prop = min.prop, 
      min.total.count = min.total.count)
    )

  # sample-level weights based on cell counts and mean library size
  # if (!is.null(weights) & !is(weights, "function")) {

  #   if( ! all(rownames(y)[keep] %in% rownames(weights)) ){
  #     msg <- "All genes retained in count matrix must be present in weights matrix.\n Make sure getExprGeneNames() and processAssays() use same parameter values."
  #     stop(msg)
  #   }

  #   precWeights <- weights[rownames(y)[keep], colnames(y)]
  # } else {
  #   precWeights <- rep(1, ncol(y))
  # }

  # if no genes are kept
  if (sum(keep) == 0) {
    return(NULL)
  }

  design = model.matrix(formula, data)
  y = y[keep, ,keep.lib.sizes=FALSE]
  y = (y)
  # browser()
  y = estimateDisp(y,design)
  geneExpr=y
  # geneExpr <- voomWithDreamWeights(y[keep, ], formula, data, 
  #   weights = precWeights, 
  #   rescaleWeightsAfter = rescaleWeightsAfter, BPPARAM = BPPARAM, 
  #   save.plot = TRUE, 
  #   quiet = quiet, 
  #   span = span, 
  #   hideErrorsInBackend = TRUE)

  # if no genes are succeed
  if (nrow(geneExpr) == 0) {
    return(NULL)
  }

  # save formula used after dropping constant terms
  if (!is.null(geneExpr)) geneExpr$formula <- formula
  if (!is.null(geneExpr)) geneExpr$isCounts <- isCounts

  # # geneExpr$E <- edgeR::cpm(y[rownames(geneExpr), ], log=TRUE, prior.count = 0.5)
  # tmp <- edgeR::cpm(y[rownames(geneExpr), ], log=TRUE, prior.count = 0.5)
  # tmp=tmp[,colnames(geneExpr)]
  # # identical(tmp,geneExpr$E)
  # # summary(c(tmp-geneExpr$E))
  # geneExpr$E=tmp

  geneExpr
}

environment(processOneAssay_edgeR) <- rlang::new_environment(parent = environment(dreamlet:::processOneAssay))
