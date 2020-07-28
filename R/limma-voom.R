#library(ggbio)
#library(edgeR)
#library(limma)

#' R6 class for limma voom analysis
#'
#' @param salmon salmon count data
#' @param groupcol category to use for DGE contrasts
#' @param outdir directory to write DGE results to
#'
#' @return
#' @export
#'
#' @examples
#' voom <- LimmaVoomClass$new(salmon, groupcol, outdir=outdir)

LimmaVoomClass <- R6::R6Class("LimmaVoomClass",

  public = list(

    salmon = NULL,
    samples = NULL,
    design = NULL,
    fit = NULL,

    initialize = function(salmon, groupcol, outdir)
    {
      samples <- salmon$project$getSamplesByGroup(groupcol)
      salmon <- salmon$subset(samples=samples$name)

      # https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
      dge <- edgeR::DGEList(counts=salmon$txi$counts, group=samples$group)

      # create a design matrix
      #design <- model.matrix(~ 0+samples$group)
      #colnames(design) <- levels(samples$group)
      design <- model.matrix(~ samples$group)
      colnames(design) <- levels(samples$group)

      # The next step is to remove rows that consistently have zero or very low counts
      keep <- edgeR::filterByExpr(dge, design)
      dge <- dge[keep, , keep.lib.sizes=FALSE]

      v <- limma::voom(dge, design, plot=FALSE)
      fit <- limma::lmFit(v, design)
      fit <- limma::eBayes(fit)

      tbl <- limma::topTable(fit, coef=ncol(design), sort.by='logFC', number=15000)
      groupname <- tolower(substring(colnames(design)[2], 14))
      outfile <- paste0(outdir, '/table-voom-group-',groupcol,'.txt')
      hlsgr::writeTable(tbl, outfile, row.names=TRUE)

      self$salmon <- salmon
      self$samples <- samples
      self$design <- design
      self$fit <- fit
      invisible(self)
    },

    getResults = function(p=NULL, padj=0.05, logfc=1.5)
    {
      results <- limma::topTable(self$fit, coef=ncol(self$design), sort.by='logFC')
      results <- na.omit(results)
      if (!is.null(p))
        results <- results[results$P.Value < p, ]
      if (!is.null(padj))
        results <- results[results$adj.P.Val < padj, ]
      if (!is.null(logfc))
        results <- results[abs(results$logFC) > logfc, ]
      return(results)
    }
  ),

  private = list(

  )
)
