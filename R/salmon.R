#' R6 class for importing salmon count data
#'
#' @param project NgsProjectClass object
#' @param quantdir folder containing salmon quant data
#'
#' @return
#' @export
#'
#' @examples
#' salmon <- SalmonClass$new(project, quantdir)

SalmonClass <- R6::R6Class("SalmonClass",

  public = list(

    project = NULL,
    samples = NULL,
    quantdir = NULL,
    files = NULL,
    txi = NULL,

    initialize = function(project, quantdir)
    {
      print(paste0('quantdir=', quantdir))
      self$project <- project
      self$samples <- project$samples
      self$quantdir <- quantdir

      self$files <- private$listQuantFiles(project, quantdir)
      tx_exp <- tximport::tximport(self$files, type='salmon', txOut=TRUE) # dtuScaledTPM
      self$txi <- tximport::summarizeToGene(tx_exp, tx2gene=hlsgr::txdb$tx2gene, ignoreTxVersion=TRUE, countsFromAbundance='scaledTPM')

      invisible(self)
    },

    save = function(filename)
    {
      saveRDS(self, file = filename)
      invisible(self)
    },

    subset = function(samples=NULL, genes=NULL)
    {
      #print(paste0('samples=', samples, ' genes=', genes))
      #https://adv-r.hadley.nz/r6.html
      salmon <- self$clone(deep=TRUE)

      txi <- salmon$txi
      if (is.null(samples))
        samples <- salmon$samples$name
      if (is.null(genes))
        genes <- rownames(txi$counts)

      genes <- rownames(txi$counts)[rownames(txi$counts) %in% genes]

      txi$abundance <- txi$abundance[genes, samples]
      txi$counts <- txi$counts[genes, samples]
      txi$length <- txi$length[genes, samples]

      salmon$txi <- txi
      salmon$samples <- salmon$samples[samples,]
      return(salmon)
    }
  ),

  private = list(

    # make a list of the quant.sf files in quantdir
    listQuantFiles = function(project, quantdir)
    {
      # make a list of the quant.sf files under the current directory
      files <- file.path(quantdir, project$samples$name, "quant.sf")

      # confirm that the files exist
      all(file.exists(files))

      # tag each file with its sample name
      names(files) <- project$samples$name

      # check to make sure there is a corresponding quant.sf for each sample
      for (sample in project$samples$name)
      {
        if (!file.exists(files[[sample]]))
          R.oo::throw(concat(files[[sample]],' not found'))
      }

      return(files)
    }
  )
)


getSalmonDeSeq2Results <- function()
{
  salmon <- hlsgr::readSalmonData(paste0(quantdir, '/salmon.rds'))
  levels <- c("Control", "Case")

  txi <- salmon$txi
  samples <- salmon$samples
  samples <- samples[!is.na(samples$group),]
  group <- samples$group
  group <- factor(group, levels=levels)

  result <- hlsgr::analyzeSalmonDataDeSeq2(salmon, outdir, protein_only=FALSE)
  return(result)
}

getDeSeq2ResultTable <- function(result, p=0.05, logfc=1)
{
  res_sig <- subset(result$res, padj<p)
  res_lfc_up <- subset(res_sig, log2FoldChange > logfc)
  res_lfc_down <- subset(res_sig, log2FoldChange < logfc*-1)

  sig <- c("up(logFC > 1)", "notsig", "down(logFC < -1)")
  sig_counts <- c(length(rownames(res_lfc_up)), length(rownames(result$res)), length(rownames(res_lfc_down)))
  sig_table <- data.frame(regulation = sig, counts = sig_counts)
  return(sig_table)
}
