#' R6 class for DESeq2 analysis counted by stringtie
#'
#' @param project NgsProjectClass object
#' @param groupcol category to use for DGE contrasts
#' @param stdir directory containing stringtie gene count matrix
#' @param outdir directory to write DGE results to
#'
#' @return
#' @export
#'
#' @examples
#' stringtietodeseq2 <- StringtieToDeSeq2Class$new(project=project, groupcol=groupcol, stdir=stdir, outdir=outdir)

StringtieToDeSeq2Class <- R6::R6Class("StringtieToDeSeq2Class",

  public = list(

    project = NULL,
    count = NULL,
    samples = NULL,
    dds = NULL,
    design = NULL,
    results = NULL,
    groupname = NULL,
    outfile = NULL,
    vsd = NULL,
    mds = NULL,
    annot_col = NULL,

    initialize = function(project, groupcol, stdir, outdir)
    {
      samples <- project$getSamplesByGroup(groupcol)
      countData <- as.matrix(read.csv(paste(stdir, "/gene_count_matrix.csv", sep=""), row.names="gene_id", check.names=F))
      countData <- countData[, rownames(samples)]

      # create the design matrix
      design <- model.matrix(~samples$group)
      print(design)

      # https://angus.readthedocs.io/en/2019/diff-ex-and-viz.html
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = samples, design = design)
      dds <- DESeq2::DESeq(dds)
      results <- DESeq2::results(dds)
      results <- results[order(results$log2FoldChange, decreasing=TRUE), ]

      #groupname <- tolower(substring(colnames(design)[2], 14))
      outfile <- paste0(outdir, '/table-stringtie-deseq2-group-',groupcol,'.txt')
      writeTable(as.data.frame(results), outfile, row.names=TRUE)

      vsd <- DESeq2::vst(dds)
      sample_dists <- SummarizedExperiment::assay(vsd) %>% t() %>% dist() %>% as.matrix()

      mdsData <- data.frame(cmdscale(sample_dists))
      mds <- cbind(mdsData, as.data.frame(SummarizedExperiment::colData(vsd))) # combine with sample data

      annot_col <- samples[,c('sample', 'group')] %>% dplyr::select(group) %>% as.data.frame()

      self$project <- project
      self$samples <- samples
      self$dds <- dds
      self$results <- results
      self$mds <- mds
      self$vsd <- vsd
      self$annot_col <- annot_col
      invisible(self)
      #return(list(salmon=salmon, samples=samples, dds=dds, res=res, res_lfc=res_lfc, mds=mds, vsd=vsd, annot_col=annot_col))

    },

    getResults = function(p=NULL, padj=0.05, logfc=1.5)
    {
      results <- self$results
      results <- na.omit(results)
      if (!is.null(p))
        results <- results[results$p < p, ]
      if (!is.null(padj))
        results <- results[results$padj < padj, ]
      if (!is.null(logfc))
        results <- results[abs(results$log2FoldChange) > logfc, ]
      return(results)
    },

    makeTable = function(num=50, padj=1, lfc=1.5)
    {
      tbl <- self$results
      tbl <- na.omit(tbl)
      tbl <- tbl[tbl$padj<padj, ]


      if (lfc > 0)
      {
        tbl <- tbl[tbl$log2FoldChange>lfc, ]
        tbl <- tbl[order(tbl$log2FoldChange, decreasing=TRUE), ]
      }
      if (lfc < 0)
      {
        tbl <- tbl[tbl$log2FoldChange<lfc, ]
        tbl <- tbl[order(tbl$log2FoldChange), ]
      }



      genetable <- read.csv('//10.37.17.251/out/test/multiomics_naka/ref/GRCh38/gencode.v34_gene_annotation_table.txt', sep="\t")
      tbl$ensembl_id <- rownames(tbl)

      tbl <- merge(as.data.frame(tbl), genetable[, c("ensembl_id", "GeneSymbol")], by.x='ensembl_id', by.y='ensembl_id', sort=FALSE)
      tbl <- tbl[, c("GeneSymbol", "log2FoldChange", "pvalue", "padj" )]

      tbl <- tbl[!duplicated(tbl), ]

      tbl <- head(tbl, n=num)
      return(tbl)
    }

    # getSigCounts <- function(padg=0.05, logfc=1)
    # {
    #   results <- results[results$padj<padj, ]
    #   res_lfc_up <- subset(res_sig, log2FoldChange > logfc)
    #   res_lfc_down <- subset(res_sig, log2FoldChange < logfc*-1)
    #
    #   sig <- c("up(logFC > 1)", "notsig", "down(logFC < -1)")
    #   sig_counts <- c(length(rownames(res_lfc_up)), length(rownames(result$res)), length(rownames(res_lfc_down)))
    #   sig_table <- data.frame(regulation = sig, counts = sig_counts)
    #   return(sig_table)
    # }
  ),

  private = list(

  )
)
