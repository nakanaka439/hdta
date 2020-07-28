#library(R6)
#library(mixOmics)

#' R6 class for mixomics analysis
#'
#' @param result DESeq2 analysis result
#' @param datTraits clinicaldata
#'
#' @return
#' @export
#'
#' @examples
#' mixomixs <- MixomicsClass$new(result, datTraits)

MixomicsClass <- R6::R6Class("MixomicsClass",
  public = list(

    data = NULL,
    samples = NULL,
    design = NULL,
    ncomp = NULL,
    sgccda.res = NULL,
    perf.diablo = NULL,
    list.keepX = NULL,

    initialize = function(result, datTraits, genecode=FALSE) {
      count_data <- DESeq2::counts(result$dds, normalized=TRUE)
      keep <- rowSums(counts(result$dds)) >= 70
      count_data <- count_data[keep,]

      if(genecode){
        fixed_ids <- genecode2ensembl(rownames(count_data))
        rownames(count_data) <- fixed_ids
      }else{
        count_data <- count_data
      }

      samples <- datTraits
      samples <- samples[!is.na(samples$group),]

      max_genes <- 10000
      mRNA <- t(count_data[order(apply(count_data, 1, mad), decreasing = T)[1:max_genes],])

      metabolomics <- read.csv(paste0(metabolomicsdir, "/metabolome_normalization.txt"), header=T, check.names=F, row.names=1)

      multi <- list(
        mRNA = mRNA,
        metabolomics = t(metabolomics))

      lapply(multi, dim)

      ids <- base::Reduce(intersect, list(rownames(multi$mRNA),
                                          rownames(multi$metabolomics)))

      multi$mRNA <- multi$mRNA[ids, ]
      multi$metabolomics <- multi$metabolomics[ids, ]

      lapply(multi, dim)

      samples <- samples[ids,]

      samples$group <- factor(samples$group, levels=c("Control", "Case"))

      self$data <- multi
      self$samples <- samples
      self$design <- private$createDesign(multi)

      invisible(self)
    },

    evaluatePerformance = function()
    {
      self$sgccda.res <- mixOmics::block.splsda(X = self$data,
                                                Y = self$samples$group,
                                                ncomp = 5, design = self$design)

      # this code takes a few minutes to run
      self$perf.diablo <- mixOmics::perf(self$sgccda.res, validation = 'Mfold',
                                         folds = 10, nrepeat = 10)
    },

    plotPerformance = function()
    {
      plot(self$perf.diablo)
    },

    tune = function()
    {
      # $choice.ncomp indicates an optimal number of components for the final DIABLO model.
      #self$perf.diablo$choice.ncomp$WeightedVote

      # optimal number of components
      self$ncomp <- self$perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

      # model tuning
      #test.keepX = list(mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
      #		metabolomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))

      test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
      		metabolomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))


      tune.TCGA = mixOmics::tune.block.splsda(X = self$data, Y = self$samples$group,
              ncomp = self$ncomp, test.keepX = test.keepX, design = self$design,
              validation = 'Mfold', folds = 10, nrepeat = 3, cpus = 2,
              dist = "centroids.dist")

      self$list.keepX <- tune.TCGA$choice.keepX
      #list.keepX = list(mRNA = c(9, 6, 6), metabolomics = c(30, 7, 5))
      print(self$list.keepX)
      #mRNA: 15  5 20 20, metabolomics: 10 20 10
      invisible(self)
    },

    getFinalModel = function()
    {
      self$sgccda.res <- mixOmics::block.splsda(X = self$data,
        Y = self$samples$group, ncomp = self$ncomp,
        keepX = self$list.keepX, design = self$design)
      #self$sgccda.res$design
      invisible(self)
    },

    plotDiablo = function(ncomp = 1)
    {
      mixOmics::plotDiablo(self$sgccda.res, ncomp = ncomp)
      invisible(self)
    },

    plotIndiv = function()
    {
      if(self$ncomp > 1){
        mixOmics::plotIndiv(self$sgccda.res, ind.names = FALSE,
                            legend = TRUE, title = 'DIABLO')
      }else{
        print("Cannot plot because ncomp=1.")
      }
      invisible(self)
    },

    plotArrow = function()
    {
      if(self$ncomp > 1){
        mixOmics::plotArrow(self$sgccda.res, ind.names = FALSE,
                            legend = TRUE, title = 'DIABLO')
      }else{
        print("Cannot plot because ncomp=1.")
      }
      invisible(self)
    },

    plotVar = function()
    {
      if(self$ncomp > 1){
        mixOmics::plotVar(self$sgccda.res, var.names = FALSE, style = 'graphics',
                          legend = TRUE, pch = c(16, 17), cex = c(2,2),
                          col = c('darkorchid', 'lightgreen'))
      }else{
        print("Cannot plot because ncomp=1.")
      }
      invisible(self)
    },

    plotCircos = function()
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      mixOmics::circosPlot(sgccda.res, cutoff = 0.4, line = TRUE,
                           color.blocks= c('darkorchid','lightgreen'),
                           color.cor = c("chocolate3","grey20"),
                           size.labels = 1.5, size.variables = 0.5)
      invisible(self)
    },

    plotLoadings = function()
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      for(i in 1:self$ncomp)
      {
        mixOmics::plotLoadings(sgccda.res, comp = i, contrib = 'max',
          method = 'median', size.name = 1.0)#, col=self$col
      }
      invisible(self)
    },

    plotNetwork = function(cutoff=0.4, output=FALSE)
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      net <- mixOmics::network(sgccda.res, blocks = c(1,2),
                               color.node = c('darkorchid', 'lightgreen'), cutoff = cutoff)
      if(output){
        library(igraph)
        write.graph(net$gR, file = paste(outdir, "myNetwork.gml", sep=""), format = "gml")
      }else{
        print("No output")
      }
      #
      invisible(self)
    },

    plotHeatmap = function()
    {
      sgccda.res <- private$annotateResults(self$sgccda.res)
      mixOmics::cimDiablo(sgccda.res, margins=c(20,20), color.blocks=c('pink', 'green'))
      invisible(self)
    },

    cytoscapeImportTable = function()
    {
      if (self$ncomp > 1){
      	loading_score_comp1_mrna <- data.frame(comp1 = subset(self$sgccda.res$loading$mRNA[, 1], abs(self$sgccda.res$loading$mRNA[, 1]) > 0))
      	loading_score_comp1_metabo <- data.frame(comp1 = subset(self$sgccda.res$loading$metabolomics[, 1], abs(self$sgccda.res$loading$metabolomics[, 1]) > 0))
      	loading_score_comp1 <- rbind(loading_score_comp1_mrna, loading_score_comp1_metabo)
      	loading_score_comp1 <- data.frame(comp1=loading_score_comp1, abs_comp1=abs(loading_score_comp1), comp_id=rownames(loading_score_comp1))
      	colnames(loading_score_comp1) <- c("comp1", "abs_comp1", "comp_id")

      	component_name <- read.csv(paste(datadir, "/component_name.txt", sep=""), sep="\t")
      	colnames(component_name) <- c("comp_id", "comp_name")

      	loading_score_comp1_name <- merge(loading_score_comp1, unique(component_name), by="comp_id", sort=FALSE)
      	write.csv(loading_score_comp1_name, paste(outdir, "/loading_score_comp1_name.csv", sep=""))
      }else{
      	loading_score_comp1_mrna <- subset(self$sgccda.res$loading$mRNA, abs(self$sgccda.res$loading$mRNA) > 0)
      	loading_score_comp1_metabo <- subset(self$sgccda.res$loading$metabolomics, abs(self$sgccda.res$loading$metabolomics) > 0)
      	loading_score_comp1 <- rbind(loading_score_comp1_mrna, loading_score_comp1_metabo)
      	loading_score_comp1 <- data.frame(comp1=loading_score_comp1, abs_comp1=abs(loading_score_comp1), comp_id=rownames(loading_score_comp1))
      	colnames(loading_score_comp1) <- c("comp1", "abs_comp1", "comp_id")

      	component_name <- read.csv(paste(datadir, "/component_name.txt", sep=""), sep="\t")
      	colnames(component_name) <- c("comp_id", "comp_name")

      	loading_score_comp1_name <- merge(loading_score_comp1, unique(component_name), by="comp_id", sort=FALSE)
      	write.csv(loading_score_comp1_name, paste(outdir, "/loading_score_comp1_name.csv", sep=""))
      }
    }
  ),

  private = list(
    createDesign = function(data)
    {
      design <- matrix(0.1, ncol = length(data), nrow = length(data),
                       dimnames = list(names(data), names(data)))
      diag(design) = 0
      return(design)
    },

    annotateResults = function(sgccda.res.orig)
    {
      sgccda.res <- sgccda.res.orig

      ensembl_ids <- rownames(sgccda.res$loadings$mRNA)
      genes <-  lookupEnsemblIds(ensembl_ids)
      rownames(sgccda.res$loadings$mRNA) <- genes
      colnames(sgccda.res$X$mRNA) <- genes
      sgccda.res$names$colnames$mRNA <- genes

      hmt_ids <- rownames(sgccda.res$loading$metabolomics)
      hmt_names <- lookupHmtIds(hmt_ids)
      rownames(sgccda.res$loadings$metabolomics) <- hmt_names
      colnames(sgccda.res$X$metabolomics) <- hmt_names
      sgccda.res$names$colnames$metabolomics <- hmt_names

      return(sgccda.res)
    }
  )
)
