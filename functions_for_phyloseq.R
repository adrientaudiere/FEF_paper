################################################################################
#' Transform the otu_table of a \code{\link{phyloseq-class}} object into a binary otu_table. Usefull to test if the results are not biaised by sequences bias that appended during PCR or NGS pipeline.  
#' 
#' @param physeq (Required): a \code{\link{phyloseq-class}} object.
#' @param minNumber (default = 1): the minimum number of sequences to put 
#' a 1 in the otu table.
#' @author Adrien Taudiere
#'
#' @return A \code{physeq} object with only 0/1 in the OTU table
#'
#' @examples 
#' data(enterotype)
#' #enterotype.bin <- as.binaryOtuTable(enterotype)
as.binaryOtuTable <- function(physeq, minNumber = 1) {
  if (!inherits(physeq, "phyloseq")) {
    stop("physeq must be a phyloseq object")
  }
  res <- physeq
  res@otu_table[res@otu_table >= minNumber] <- 1
  res@otu_table[res@otu_table < minNumber] <- 0
  return(res)
}
################################################################################

################################################################################
#' Plot edgeR results for a phyloseq or a edgeR object.
#' @param data (Required): a \code{\link{phyloseq-class}} or a \code{\link{edgeR}} object.
#'#' @param tax_table: Require if data is a \code{\link{edgeR}} object. 
#' The taxonomic table used to find the \code{taxa} and \code{color_taxa} arguments. If data is a \code{\link{phyloseq-class}} object, data@tax_table is used.
#' @param contrast (Required):This argument specifies what comparison to extract from the object to build a results table. See \code{\link[DESeq2]{results}} man page for more details.
#' @param alpha (default = 0.01): the significance cutoff used for optimizing the independent filtering. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param taxa (default = 'Genus'): taxonomic level of interest
#' @param color_tax (default = 'Phylum'): taxonomic level used for color assignation
#' @param verbose : whether the function print some information during the computation
#' @param ... Additional arguments passed on to \code{\link[edgeR]{exactTest}} or \code{\link[ggplot2]{ggplot}}
#' 
#' @export
#' 
#' @examples 
#'\dontrun{
#'data(GlobalPatterns)
#'plot_edgeR_phyloseq(GlobalPatterns, c('SampleType', 'Soil', 'Feces'), tax_table = GlobalPatterns@tax_table, color_tax = 'Kingdom')
#'plot_edgeR_phyloseq(GlobalPatterns, c('SampleType', 'Soil', 'Feces'), taxa = 'Class', tax_table = GlobalPatterns@tax_table, color_tax = 'Kingdom')
#'}
#' @author Adrien Taudiere
#'
#' @return A \code{\link{ggplot}}2 plot representing edgeR results
#'
#' @seealso \code{\link[edgeR]{exactTest}}
#' @seealso \code{\link{plot_deseq2_phyloseq}}

plot_edgeR_phyloseq <- function(data, contrast = NULL, alpha = 0.01, taxa = "Genus", color_tax = "Phylum", 
                                verbose = TRUE, ...) {
  
  if (!inherits(data, "phyloseq")) {
    stop("data must be an object of class 'phyloseq'")
  }
  
  if (verbose) {
    message("Conversion to edgeR format")
  }
  data_edgeR <- phyloseq_to_edgeR(data, group = contrast[1])
  if (verbose) {
    message("Perform edgeR binary test")
  }
  et <- exactTest(data_edgeR, pair = c(contrast[2], contrast[3]), ...)
  
  tt <- topTags(et, n = nrow(et$table), adjust.method = "BH", sort.by = "PValue")
  res <- tt@.Data[[1]]
  sigtab <- res[(res$FDR < alpha), ]
  sigtab <- cbind(as(sigtab, "data.frame"))
  
  sigtabgen <- subset(sigtab, !is.na(taxa))
  
  d <- tapply(sigtabgen$logFC, sigtabgen[, color_tax], function(x) max(x))
  d <- sort(d, TRUE)
  sigtabgen$col_tax <- factor(as.character(sigtabgen[, color_tax]), levels = names(d))
  
  d <- tapply(sigtabgen$logFC, sigtabgen[, taxa], function(x) max(x))
  d <- sort(d, TRUE)
  sigtabgen$tax <- factor(as.character(sigtabgen[, taxa]), levels = names(d))
  
  p <- ggplot(sigtabgen, aes(x = tax, y = logFC, color = col_tax), ...) + geom_point(size = 6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
    labs(title = paste("Change in abundance for ", contrast[1], " (", contrast[2], " vs ", contrast[3], ")", sep = ""))
  
  return(p)
}
################################################################################

################################################################################
#' Convert phyloseq data to DESeq2 dds object
#'
#' No testing is performed by this function. The phyloseq data is converted
#' to the relevant \code{\link[DESeq2]{DESeqDataSet}} object, which can then be
#' tested in the negative binomial generalized linear model framework
#' of the \code{\link[DESeq2]{DESeq}} function in DESeq2 package.
#' See the
#' \href{http://joey711.github.io/phyloseq-extensions}{phyloseq-extensions}
#' tutorials for more details.
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}.
#'  Must have a \code{\link{sample_data}} component.
#'
#' @param design (Required). A \code{\link{formula}} which specifies the design of the experiment,
#'  taking the form \code{formula(~ x + y + z)}. That is, a formula with right-hand side only.
#'  By default, the functions in this package and DESeq2
#'  will use the last variable in the formula (e.g. \code{z})
#'  for presenting results (fold changes, etc.) and plotting.
#'  When considering your specification of experimental design, you will want to 
#'  re-order the levels so that the \code{NULL} set is first.
#'  For example, the following line of code would ensure that Enterotype 1 is used as the 
#'  reference sample class in tests by setting it to the first of the factor levels
#'  using the \code{\link{relevel}} function:
#'  
#'  \code{sample_data(entill)$Enterotype <- relevel(sample_data(entill)$Enterotype, "1")}
#'  
#' @param ... (Optional). Additional named arguments passed to \code{\link[DESeq2]{DESeqDataSetFromMatrix}}.
#'  Most users will not need to pass any additional arguments here.
#'  Most testing-related options should be provided in 
#'  a following call to \code{\link[DESeq2]{DESeq}}.
#'  
#' @return A \code{\link[DESeq2]{DESeqDataSet}} object.
#' 
#' @seealso
#' 
#' \code{vignette("phyloseq-mixture-models")}
#' 
#' The 
#' \href{http://joey711.github.io/phyloseq-extensions}{phyloseq-extensions}
#' tutorials.
#' 
#'  \code{\link[DESeq2]{DESeq}}
#'  
#'  \code{\link[DESeq2]{results}}
#'
#' @export
#'  
#' @examples
#'  # Check out the vignette phyloseq-mixture-models for more details.
#'  # vignette("phyloseq-mixture-models")
#'  data(soilrep)
#'  phyloseq_to_deseq2(soilrep, ~warmed)
phyloseq_to_deseq2 = function(physeq, design, ...){
  # Need to add check here for missing sample_data
  if( is.null(sample_data(physeq, FALSE)) ){
    stop("There must be sample_data present, for specifying experimental design. See ?phyloseq_to_deseq2")
  }
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  colData = data.frame(sample_data(physeq))
  # Create the DESeq data set, dds.
  if(requireNamespace("DESeq2")){
    dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData, design, ...)
    return(dds)
  }
}
################################################################################




####################################################################################
#Plot the result of a DESeq2 test
####################################################################################
#' Plot DESeq2 results for a phyloseq or a DESeq2 object.
#' @param data (Required): a \code{\link{phyloseq-class}} or a \code{\link{DESeqDataSet-class}} object.
#' @param tax_table : Required if data is a \code{\link{DESeqDataSet-class}} object. 
#' The taxonomic table used to find the \code{taxa} and \code{color_taxa} arguments. If data is a \code{\link{phyloseq-class}} object, data@tax_table is used.
#' @param contrast (Required):This argument specifies what comparison to extract from the object to build a results table. See \code{\link[DESeq2]{results}} man page for more details.
#' @param alpha (default = 0.01): the significance cutoff used for optimizing the independent filtering. If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
#' @param taxa (default = 'Genus'): taxonomic level of interest
#' @param color_tax (default = 'Phylum'): taxonomic level used for color or a color vector.
#' @param taxDepth (default = NULL): Taxonomic depth to test for differential distribution among contrast. If Null the analysis is done at the OTU (i.e. Species) level. If not Null data need to be a \code{\link{phyloseq-class}} object.
#' @param verbose : whether the function print some information during the computation
#' @param ... Additional arguments passed on to \code{\link[DESeq2]{DESeq}} or \code{\link[ggplot2]{ggplot}}
#'
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 results
#' @export
#'
#' @examples
#'\dontrun{
#'data("GlobalPatterns")
#'GlobalPatterns_Proteobacteria <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 2] == 'Proteobacteria')
#'res_deseq2 <- DESeq(phyloseq_to_deseq2(GlobalPatterns_Proteobacteria, ~ SampleType), test = 'Wald', fitType = 'local')
#'plot_deseq2_phyloseq(res_deseq2, c('SampleType', 'Soil', 'Feces'), tax_table = GlobalPatterns_Proteobacteria@tax_table, color_tax = 'Kingdom')
#'plot_deseq2_phyloseq(res_deseq2, c('SampleType', 'Soil', 'Feces'), tax_table = GlobalPatterns_Proteobacteria@tax_table, color_tax = c("red", "black"), alpha = 0.7))
#'plot_deseq2_phyloseq(GlobalPatterns_Proteobacteria, c('SampleType', 'Soil', 'Feces'), taxDepth = 'Family', taxa = 'Family',  color_tax = 'Class')
#'}
#' @author Adrien Taudiere
#' 
#' @return A \code{\link{ggplot}}2 plot representing DESeq2 results
#'
#' @seealso \code{\link[DESeq2]{DESeq}}
#' @seealso \code{\link[DESeq2]{results}}
#' @seealso \code{\link{plot_edgeR_phyloseq}}

plot_deseq2_phyloseq <- function(data, contrast = NULL, tax_table = NULL, alpha = 0.01, 
                                 taxa = "Genus", color_tax = "Phylum", taxDepth = NULL, verbose = TRUE, ...) {
  
  if (!inherits(data, "phyloseq")) {
    if (!inherits(data, "DESeqDataSet")) {
      stop("data must be an object of class 'phyloseq' or 'DESeqDataSet'")
    }
  } else {
    # Calculate new dataset given the Taxa depth if taxDepth is not null
    if (!is.null(taxDepth)) {
      data_TAX <- data
      data_TAX@otu_table <- otu_table(apply(data@otu_table, 2, function(x) tapply(x, 
                                                                                  data@tax_table[, taxDepth], sum)), taxa_are_rows = T)
      data_TAX@tax_table <- tax_table(apply(data@tax_table[, 1:match(taxDepth, colnames(data@tax_table))], 
                                            2, function(x) X <- tapply(x, data@tax_table[, taxDepth], function(xx) xx[1])))
      data_TAX@refseq <- NULL
      data <- data_TAX
      if (is.na(match(taxa, colnames(data@tax_table)))) {
        taxa <- taxDepth
      }
    }
    
    if (is.null(tax_table) & inherits(data, "phyloseq")) {
      tax_table <- data@tax_table
    }
    
    if (verbose) {
      message("Conversion to Deseq2 format.")
    }
    data_deseq2 <- phyloseq_to_deseq2(data, as.formula(paste("~", contrast[1])))
    
    if (verbose) {
      message("Calculation of Deseq2 results.")
    }
    data <- DESeq(data_deseq2, test = "Wald", fitType = "parametric", quiet = !verbose, 
                  ...)
  }
  
  # Calcul deseq2 results
  res <- results(data, contrast = contrast)
  
  d <- res[which(res$padj < alpha), ]
  
  if (dim(d)[1] == 0) {
    message("None taxa present significant distribution pattern through contrast.")
    return("None taxa present significant distribution pattern through contrast.")
  }
  d <- cbind(as(d, "data.frame"), as(tax_table[rownames(d), ], "matrix"))
  
  # Compute colors
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
    })
  }
  
  if (!sum(areColors(color_tax)) > 0) {
    x <- tapply(d$log2FoldChange, d[, color_tax], function(x) max(x))
    x <- sort(x, TRUE)
    d$col_tax <- factor(as.character(d[, color_tax]), levels = names(x))
  } else {
    d$col_tax <- rep(color_tax, length = dim(d)[1])
  }
  
  # Compute log2FoldChange values
  x <- tapply(d$log2FoldChange, d[, taxa], function(x) max(x))
  x <- sort(x, TRUE)
  d$tax <- factor(as.character(d[, taxa]), levels = names(x))
  
  if (!sum(areColors(color_tax)) > 0) {
    p <- ggplot(d, aes(x = tax, y = log2FoldChange, color = col_tax), ...) + geom_point(size = 6) + 
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + labs(title = paste("Change in abundance for ", 
                                                                                                  contrast[1], " (", contrast[2], " vs ", contrast[3], ")", sep = ""))
  } else {
    p <- ggplot(d, aes(x = tax, y = log2FoldChange), ...) + geom_point(size = 6, color = d$col_tax) + 
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + labs(title = paste("Change in abundance for ", 
                                                                                                  contrast[1], " (", contrast[2], " vs ", contrast[3], ")", sep = ""))
  }
  
  return(p)
}
################################################################################



####################################################################################
#Plot the result of a mt test
####################################################################################
#' Plot the result of a mt test (\code{\link[phyloseq]{mt}})
#' @param mt (Required): result of a mt test
#' @param alpha (default = 0.05): Choose the cut off p-value to plot taxa
#' @param color_tax : A taxonomic level to color the points
#' @param taxa : The taxonomic level choose for x-positioning
#' @examples
#' \dontrun{
#' data("GlobalPatterns")
#' res = mt(GlobalPatterns, "SampleType", test="f")
#' plot_mt(res, color_tax = "Phylum") + scale_color_hue()
#' }
#' @author Adrien Taudiere
#' 
#' @return a \code{\link{ggplot}}2 plot of result of a mt test
#' 
#' @seealso \code{\link[phyloseq]{mt}}

plot_mt <- function(mt = NULL, alpha = 0.05, color_tax = "Class", taxa = "Species"){
  d <- mt[mt$plower < alpha,]
  d$tax_col <- factor(as.character(d[, color_tax]))
  d$tax_col[is.na(d$tax_col)] <- "unidentified"
  d$tax <- as.character(d[, taxa])
  d$tax[is.na(d$tax)] <- "unidentified"
  d$tax <- factor(d$tax, levels = unique(factor(as.character(d[,"Species"]))[rev(order(d$teststat))]))
  
  p <- ggplot(d, aes(x = tax, y = teststat, color = tax_col)) + geom_point(size = 6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) 
  p
}
####################################################################################


################################################################################
#' Plot accumulation curves for \code{\link{phyloseq-class}} object
#' 
#' @param physeq (Required): a \code{\link{phyloseq-class}} object.
#' @param fact (Required): Name of the factor in physeq@sam_data used to plot different lines
#' @param nbSeq (logical): Either plot accumulation curves using sequences or using samples 
#' @param step (integer): distance among points calculated to plot lines. A
#'  low value give better plot but is more time consuming. Only use if nbSeq = TRUE
#' @param by.fact (logical): First merge the OTU table by factor to plot only one line by factor 
#' @param ci.col : Color vector for confidence intervall. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to cange plot.
#' @param col : Color vector for lines. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to change plot. 
#' @param lwd  (default = 3): thickness for lines. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to cange plot. 
#' @param leg (logical): Plot legend or not. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to cange plot.
#' @param printSamplesNames (logical): Print samples names or not? Only use if nbSeq = TRUE.
#' @param CI (default = 2) : Confidence intervall value used to multiply the standard error to plot confidence intervall
#' @param ... Additional arguments passed on to \code{\link{ggplot}} if nbSeq = TRUE 
#' or \code{\link{plot}} if nbSeq = FALSE
#'
#' @examples 
#' data("GlobalPatterns")
#' GlobalPatterns_Archaea <- subset_taxa(GlobalPatterns,
#'  GlobalPatterns@tax_table[, 1] == 'Archaea')
#' \dontrun{
#' accu_plot(GlobalPatterns_Archaea, 'SampleType', nbSeq = TRUE, by.fact = TRUE)
#' }
#' 
#' @return A \code{\link{ggplot}}2 plot representing the richness 
#' accumulation plot if nbSeq = TRUE, else, if nbSeq = FALSE
#' return a base plot.
#'  
#' @importFrom vegan specaccum
#' @importFrom vegan rarefy
#' @importFrom plyr ldply
#'
#' @author Adrien Taudiere
#' @seealso \code{\link[vegan]{specaccum}}
accu_plot <- function(physeq, fact = NULL, nbSeq = TRUE, step = NULL, by.fact = FALSE,
                      ci.col = NULL, col = NULL, lwd = 3, leg = TRUE, 
                      printSamplesNames = FALSE, CI = 2, ...) {
  
  if (!inherits(physeq, "phyloseq")) {
    stop("physeq must be a phyloseq object")
  }
  
  if (!nbSeq){
    factor.interm <- eval(parse(text = paste("physeq@sam_data$", fact, sep = "")))
    factor.interm <- as.factor(factor.interm)
    
    physeq_accu <- as.matrix(t(physeq@otu_table))
    physeq_accu[physeq_accu > 0] <- 1
    accu_all <- specaccum(physeq_accu)
    
    accu <- list()
    for (i in 1:nlevels(factor.interm)) {
      accu[[i]] <- specaccum(physeq_accu[factor.interm == levels(factor.interm)[i], ])
      #print(paste(round(i/nlevels(factor.interm) * 100), "%"))
    }
    
    funky.color <- colorRampPalette(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
                                      "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"))
    
    if (is.null(col)) {
      col <- funky.color(nlevels(factor.interm) + 1)
    }
    if (is.null(ci.col)) {
      transp <- function (col, alpha = 0.5) {
        res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
        return(res)
      }    
      ci.col <- transp(funky.color(nlevels(factor.interm) + 1), 0.3)
    }
    
    plot(accu_all, ci.type = "poly", ci.col = ci.col[1], col = col[1], lwd = lwd, ci.lty = 0, 
         xlab = "Sample", ...)
    
    for (i in 1:nlevels(factor.interm)) {
      lines(accu[[i]], ci.type = "poly", ci.col = ci.col[i + 1], col = col[i + 1], lwd = lwd, 
            ci.lty = 0)
    }
    if (leg) {
      legend("bottomright", c("all", levels(factor.interm)), col = col, lty = 1, lwd = 3)
    }
  }
  
  if (nbSeq) {
    FACT <- as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))
    
    if (!by.fact) {
      x <- t(physeq@otu_table)
    } else {
      x <- apply(physeq@otu_table, 1, function(x) tapply(x, FACT, sum))
    }
    
    tot <- rowSums(x)
    nr <- nrow(x)
    
    if (is.null(step)) {step = round(max(tot) / 30, 0)}
    if (is.null(step)) {step = 1}
    
    res <- list()
    n_max <- seq(1, max(tot), by = step)
    out <- lapply(seq_len(nr), function(i) {
      n <- seq(1, tot[i], by = step)
      if (n[length(n)] != tot[i]) {
        n <- c(n, tot[i])
      }
      res_interm <- rarefy(x[i, ], n, se = TRUE)
      res <- cbind(as.matrix(res_interm)[1,], as.matrix(res_interm)[2,])
      return(res)
    })
    
    names(out) <- names(tot)
    
    df <- ldply(out, data.frame)
    
    cond <- c()
    for (i in 1:nlevels(as.factor(df$.id))) {
      cond <- c(cond, 1:table(df$.id)[i])
    }
    
    df$x <- n_max[cond]
    
    if (!by.fact) {
      df$fact <- as.factor(unlist(unclass(physeq@sam_data[match(df$.id, sample_names(physeq)), fact])[fact]))
    } else {
      df$fact <- df$.id
    }
    
    df$ymin <- df$X1 - df$X2 * CI
    df$ymin[is.na(df$ymin)] <- df$X1[is.na(df$ymin)]
    df$ymax <- df$X1 + df$X2 * CI
    df$ymax[is.na(df$ymax)] <- df$X1[is.na(df$ymax)]
    dff <- data.frame(matrix(nrow = length(tot)))
    dff$xlab <- tapply(df$x, df$.id, max) 
    dff$xlab <- dff$xlab + max(dff$xlab, na.rm = TRUE)/20
    dff$ylab <- tapply(df$X1, df$.id, max)
    dff$.id <-  names(dff$ylab)
    p <- ggplot(data = df, aes(x = x, y = X1, group = .id, col = fact)) +  
      geom_ribbon(aes(ymin = ymin, ymax = ymax, col = NULL, fill = fact), alpha = 0.2) +
      geom_line() + xlab("Number of sequences") + ylab("Number of OTUs (with standard error)") 
    
    if (printSamplesNames) {
      p + geom_text(data = dff, aes(x = xlab , y = ylab, label = .id, col = NULL))
    } else {p}
    return(p)
  }
}
################################################################################

################################################################################
#' Plot OTU circle for \code{\link{phyloseq-class}} object
#' @param physeq (Required): a \code{\link{phyloseq-class}} object.
#' @param fact (Required): Name of the factor to cluster samples by modalities. Need to be in \code{physeq@sam_data}.
#' @param taxa (Default:'Order'): Name of the taxonomic rank of interest
#' @param nbSeq (Default: TRUE): Represent the number of sequences or the number of OTUs (nbSeq = FALSE)
#' @param rarefy (logical): Does each samples modalities need to be rarefy in order to compare them with the same amount of sequences?
#' @param min.prop.tax (Default: 0.01): The minimum proportion for taxon to be ploted
#' @param min.prop.mod (Default: 0.1) : The minimum proportion for modalities to be ploted
#' @param gap.degree : Gap between two neighbour sectors. It can be a single value or a vector. If it is a vector, the first value corresponds to the gap after the first sector.
#' @param start.degree : The starting degree from which the circle begins to draw. Note this degree is measured in the standard polar coordinate which means it is always reverse-clockwise.
#' @param row.col : Color vector for row
#' @param grid.col : Grid colors which correspond to sectors. The length of the vector should be either 1 or the number of sectors. It's preferred that grid.col is a named vector of which names correspond to sectors. If it is not a named vector, the order of grid.col corresponds to order of sectors.
#' @param log10trans (logical): Should sequence be log10 transformed (more precisely by log10(1+x))?
#' @param ... Additional arguments passed on to \code{\link[circlize]{chordDiagram}} or \code{\link[circlize]{circos.par}}
#'
#' @examples
#' data("GlobalPatterns")
#' # GlobalPatterns_Archaea <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == 'Archaea')
#' # otu_circle(GlobalPatterns_Archaea, 'SampleType')
#' @author Adrien Taudiere
#' 
#' @return A \code{\link{chordDiagram}} plot representing the distribution 
#' of OTUs or sequences in the different modalities of the factor fact
#' 
#' @seealso \code{\link[circlize]{chordDiagram}}
#' @seealso \code{\link[circlize]{circos.par}}

otu_circle <- function(physeq = NULL, fact = NULL, taxa = "Order", nbSeq = TRUE, rarefy = FALSE, 
                       min.prop.tax = 0.01, min.prop.mod = 0.1, gap.degree = NULL, start.degree = NULL, row.col = NULL, 
                       grid.col = NULL, log10trans = F,...) {
  
  if (!inherits(physeq, "phyloseq")) {
    stop("physeq must be an object of class 'phyloseq'")
  }
  
  if (!nbSeq) {
    physeq@otu_table[physeq@otu_table > 0] <- 1
  }
  
  taxcol <- match(taxa, colnames(physeq@tax_table))
  if (is.na(taxcol)) {
    stop("The taxa argument do not match any taxa rank in physeq@tax_table")
  }
  
  taxsamp <- match(fact, colnames(physeq@sam_data))
  if (is.na(taxsamp)) {
    stop("The samples argument do not match any sample attributes in physeq@sam_data")
  }
  
    otu_table_tax <- apply(physeq@otu_table, 2, function(x) tapply(x, physeq@tax_table[, taxcol], 
                                                                 function(xx) sum(xx, na.rm = T)))
    otu_table_ech <- apply(otu_table_tax, 1, function(x) tapply(x, physeq@sam_data[, taxsamp], 
                                                              function(xx) sum(xx, na.rm = T)))
  if (rarefy) {
    otu_table_ech_interm <- rrarefy(otu_table_ech, min(rowSums(otu_table_ech)))
    print(paste("Rarefaction by modalities deletes ", sum(otu_table_ech) - sum(otu_table_ech_interm), 
                " (", round(100 * (sum(otu_table_ech) - sum(otu_table_ech_interm))/sum(otu_table_ech), 
                            2), "%) sequences.", sep = ""))
    otu_table_ech <- otu_table_ech_interm
  }
  
  otu_table_ech <- otu_table_ech[, colSums(otu_table_ech) > 0]
  
  # Keep only taxa and modalities with a sufficient proportion (min.prop.tax,
  # min.prop.mod) to plot
  o_t_e_interm <- otu_table_ech[(rowSums(otu_table_ech)/sum(otu_table_ech)) > min.prop.mod, 
                                (colSums(otu_table_ech)/sum(otu_table_ech)) > min.prop.tax]
  if (nrow(o_t_e_interm) != nrow(otu_table_ech)) {
    print(paste("Only ", nrow(o_t_e_interm), " modalities are plot (", round(100 * 
                                                                               nrow(o_t_e_interm)/nrow(otu_table_ech), 2), "%). Use 'min.prop.mod' to plot more samples.", 
                sep = ""))
  }
  
  if (ncol(o_t_e_interm) != ncol(otu_table_ech)) {
    print(paste("Only ", ncol(o_t_e_interm), " taxa are plot (", round(100 * ncol(o_t_e_interm)/ncol(otu_table_ech), 
                                                                       2), "%). Use 'min.prop.tax' to plot more taxa", sep = ""))
  }
  otu_table_ech <- o_t_e_interm
  
  if (log10trans) {
    otu_table_ech <- apply(otu_table_ech, 2, function(x) log10(1+x))
  }
  
  
  if (is.null(gap.degree)) {
    col2keep <- rep(1, ncol(otu_table_ech) - 1)
    row2keep <- rep(1, nrow(otu_table_ech) - 1)
    gap.degree <- c(row2keep, 10, col2keep, 10)
    
  }
  if (is.null(start.degree)) {
    start.degree <- 170
  }
  
  funky.color <- colorRampPalette(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
                                    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"))
  
  if (is.null(grid.col)) {
    grid.col <- c(funky.color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
  }
  
  if (is.null(row.col)) {
    row.col <- c(funky.color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
  }
  
  circos.par(gap.degree = gap.degree, start.degree = start.degree, ...)
  chordDiagram(otu_table_ech, row.col = row.col, grid.col = grid.col, ...)
  circos.clear()
}
################################################################################

################################################################################
#' Sankey plot of \code{\link{phyloseq-class}} object
#' @param physeq (Required): a \code{\link{phyloseq-class}} object.
#' @param fact (Optional): Name of the factor to cluster samples by modalities. 
#' Need to be in \code{physeq@sam_data}.
#' @param taxa (Default: c(1:4)): a vector of taxonomic rank to plot
#' @param nbSeq (Default: FALSE): Represent the number of sequences or the number of OTUs (nbSeq = FALSE). Note that ploting the number of sequences is slower.
#' @param min.prop.tax (Default: 0): The minimum proportion for taxon to be ploted. EXPERIMENTAL. For the moment each links below the min.prop. tax is discard from the sankey network resulting in sometimes weird plot.
#' @param tax2remove : a vector of taxonomic groups to remove from the analysis (e.g. \code{c('Incertae sedis', 'unidentified')})
#' @param units : character string describing physical units (if any) for Value
#' @param Symbol2sub (default = c('\\.', '-')): vector of symbol to delete in the taxonomy 
#' @param ... Additional arguments passed on to \code{\link[networkD3]{sankeyNetwork}}
#'
#' @examples
#' data("GlobalPatterns")
#' #GlobalPatterns_Archaea <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == 'Archaea')
#' #sankey_phyloseq(GlobalPatterns_Archaea)
#' #sankey_phyloseq(GlobalPatterns_Archaea, fact = 'SampleType')
#'\dontrun{
#' sankey_phyloseq(GlobalPatterns, taxa = c(1:5), min.prop.tax = 0.01)
#' sankey_phyloseq(GlobalPatterns, taxa = c(2:6), min.prop.tax = 0.01, nbSeq = TRUE)
#'}
#' @importFrom networkD3 sankeyNetwork
#' @author Adrien Taudiere
#' 
#' @return A \code{\link{sankeyNetwork}} plot representing the taxonomic distribution
#' of OTUs or sequences. If \code{fact} is set, represent the distribution of
#' the last taxonomic level in the modalities of \code{fact}
#'  
#' @seealso \code{\link[networkD3]{sankeyNetwork}}

sankey_phyloseq <- function(physeq = NULL, fact = NULL, taxa = c(1:4), nbSeq = FALSE, 
                            min.prop.tax = 0, tax2remove = NULL, units = NULL, Symbol2sub = c("\\.", "-"), ...) {
  
  if (!inherits(physeq, "phyloseq")) {
    stop("physeq must be an object of class 'phyloseq'")
  }
  
  if (!nbSeq) {
    physeq@otu_table[physeq@otu_table > 0] <- 1
    mat.interm <- matrix()
    mat <- matrix(ncol = 3)
    colnames (mat) <- c("Var1", "Var2", "value")
    for (i in 1:(length(taxa) - 1)) {
      res.interm <- table(physeq@tax_table[, taxa[i]], physeq@tax_table[, taxa[i + 1]])
      mat.interm <- reshape2::melt(res.interm)
      mat.interm <- mat.interm[mat.interm[,3]>0,]
      mat <- rbind(mat, mat.interm)
    }
  } else if (nbSeq) {
    mat.interm <- matrix()
    mat <- matrix(ncol = 3)
    colnames (mat) <- c("Var1", "Var2", "value")
    tax_table.interm <- physeq@tax_table[rep(1:dim(physeq@tax_table)[1], times = taxa_sums(physeq))]
    
    for (i in 1:(length(taxa) - 1)) {
      res.interm <- table(tax_table.interm[, taxa[i]], tax_table.interm[, taxa[i + 1]])
      mat.interm <- reshape2::melt(res.interm)
      mat.interm <- mat.interm[mat.interm[,3]>0,]
      mat <- rbind(mat, mat.interm)
    }
  }
  
  if (!is.null(fact)) {
    NetMatrix2Links <- function(m = NULL) {
      res <- matrix(ncol = 3)
      for (i in 1:dim(m)[1]) {
        for (j in 1:dim(m)[2]) {
          if (m[i, j] > 0) {
            res <- rbind(res, c(rownames(m)[i], colnames(m)[j], m[i, j]))
          }
        }
      }
      return(res)
    }
    
    mat.interm <- apply(physeq@otu_table, 1, function(x) tapply(x, physeq@sam_data[, fact], 
                                                                sum))
    
    if (!nbSeq) {
      mat.interm <- apply(mat.interm, 1, function(x) tapply(x, physeq@tax_table[, 
                                                                                taxa[length(taxa)]], function(x) sum(x > 0)))
    } else if (nbSeq) {
      mat.interm <- apply(mat.interm, 1, function(x) tapply(x, physeq@tax_table[, 
                                                                                taxa[length(taxa)]], sum))
    }
    
    sampLinks <- NetMatrix2Links(mat.interm)
    sampLinks[, 2] <- toupper(sampLinks[, 2])
    colnames(sampLinks) <- colnames(mat)
    mat <- rbind(mat, sampLinks)
  }
  
  mat <- as.data.frame(mat[rowSums(is.na(mat)) == 0, ])
  mat[, 3] <- as.numeric(as.vector(mat[, 3]))
  mat <- mat[rowSums(is.na(mat)) == 0, ]
  
  
  if (!is.null(tax2remove)) {
    mat <- mat[!mat[, 1] %in% tax2remove, ]
    mat <- mat[!mat[, 2] %in% tax2remove, ]
  }
  
  if (min.prop.tax != 0) {
    min.nb.tax <- min.prop.tax * sum(mat[, 3])/length(taxa)
    mat <- mat[mat[, 3] >= min.nb.tax, ]
  }
  
  for (i in 1:length(Symbol2sub)) {
    mat <- apply(mat, 2, function(x) gsub(Symbol2sub[i], "", x))
  }
  
  taxSank <- list()
  namesNodes <- unique(c(as.vector(mat[, 1]), as.vector(mat[, 2])))
  namesNodes <- namesNodes[!is.na(namesNodes)]
  taxSank$nodes <- data.frame((1:length(namesNodes)) - 1, namesNodes)
  names(taxSank$nodes) <- c("code", "name")
  mat2 <- mat
  for (i in 1:nrow(taxSank$nodes)) {
    mat2[, 1] <- gsub(paste("\\<", taxSank$nodes[i, 2], "\\>", sep = ""), taxSank$nodes[i, 
                                                                                        1], mat2[, 1])
    mat2[, 2] <- gsub(paste("\\<", taxSank$nodes[i, 2], "\\>", sep = ""), taxSank$nodes[i, 
                                                                                        1], mat2[, 2])
  }
  
  taxSank$links <- apply(mat2, 2, as.numeric)
  taxSank$links <- data.frame(taxSank$links[rowSums(is.na(taxSank$links)) == 0, ])
  taxSank$nodes <- as.data.frame(as.character(taxSank$nodes[, 2]))
  names(taxSank$nodes) <- c("name")
  names(taxSank$links) <- c("source", "target", "value")
  if (is.null(units)) {
    if (!nbSeq) {
      units <- "OTUs"
    } else if (nbSeq) {
      units <- "Sequences"
    }
  }
  sankeyNetwork(Links = taxSank$links, Nodes = taxSank$nodes, Source = "source", 
                Target = "target", Value = "value", NodeID = "name", units = units, ...)
}
################################################################################

################################################################################
#' Venn diagram of \code{\link{phyloseq-class}} object
#' @param physeq (Required): a \code{\link{phyloseq-class}} object.
#' @param fact (Required): Name of the factor to cluster samples by modalities. 
#' Need to be in \code{physeq@sam_data}.
#' @param min.nb.seq (Default: 0)): minimum number of sequences by OTUs by samples 
#' to take into count this OTUs in this sample
#' @param printValues (logical) : Print (or not) the table of number of OTUs 
#' for each combination. If printValues is TRUE the object is not a ggplot object. 
#' Please use printValues = FALSE if you want to add ggplot function (cf example). 
#'
#' @importFrom grid grid.newpage
#' @importFrom grid viewport
#' @importFrom grid upViewport
#' @importFrom grid grid.draw
#' @importFrom gridExtra tableGrob
#' @importFrom grid pushViewport
#' @importFrom venneuler venneuler
#'
#' @examples
#' data("enterotype")
#'\dontrun{
#' venn_phyloseq(enterotype, fact = 'SeqTech')
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus')
#' venn_phyloseq(enterotype, fact = 'Nationality', printValues = F)
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus', printValues = F) + scale_fill_hue()
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus', printValues = F) + scale_fill_hue()
#'}
#'
#' @return A \code{\link{ggplot}}2 plot representing Venn diagramm of 
#' modalities of the argument \code{factor}
#'
#' @author Adrien Taudiere
#' @seealso \code{\link[venneuler]{venneuler}}

venn_phyloseq <- function(physeq, fact, min.nb.seq = 0, printValues = TRUE) {
  
  if (!inherits(physeq, "phyloseq")) {
    stop("physeq must be an object of class 'phyloseq'")
  }
  
  moda <- as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))
  data_Venn <- t(apply(physeq@otu_table, 1, function(x) by(x, moda, max)))
  combinations <- data_Venn > min.nb.seq
  
  e <- new.env(TRUE, emptyenv())
  cn <- colnames(combinations)
  for (i in seq.int(dim(combinations)[1])) if (any(combinations[i, ])) {
    ec <- paste(cn[combinations[i, ]], collapse = "&")
    e[[ec]] <- if (is.null(e[[ec]])) 
      1L else (e[[ec]] + 1L)
  }
  
  en <- ls(e, all.names = TRUE)
  weights <- as.numeric(unlist(lapply(en, get, e)))
  combinations <- as.character(en)
  
  table_value <- data.frame(combinations = as.character(combinations), 
                            weights = as.double(weights))
  
  VENN <- venneuler(data_Venn > min.nb.seq)
  venn_res <- data.frame(x = VENN$centers[, 1], y = VENN$centers[, 2],
                         radius = VENN$diameters/2)
  
  nmod <- nrow(venn_res)
  x1 <- list()
  for (i in 1:nmod) {
    x1[[i]] <- grep(rownames(venn_res)[i], table_value$combinations)
  }
  
  for (i in 1:nrow(table_value)) {
    table_value$x[i] <- mean(VENN$centers[, "x"][unlist(lapply(x1,
                                                               function(x) sum(x %in% i) > 0))])
    table_value$y[i] <- mean(VENN$centers[, "y"][unlist(lapply(x1, 
                                                               function(x) sum(x %in% i) > 0))])
  }
  
  df <- venn_res
  df$xlab <- df$x + (df$x - mean(df$x))
  df$ylab <- df$y + (df$y - mean(df$y))
  
  circularise <- function(d, n = 360) {
    angle <- seq(-pi, pi, length = n)
    make_circle <- function(x, y, r, Modality) {
      data.frame(x = x + r * cos(angle), y = y + r * sin(angle), Modality)
    }
    lmat <- mapply(make_circle, Modality = rownames(d), 
                   x = d[, 1], y = d[, 2], r = d[, 3], SIMPLIFY = FALSE)
    do.call(rbind, lmat)
  }
  
  circles <- circularise(df)
  
  p <- ggplot() + geom_polygon(data = circles, aes(x, y, group = Modality, fill = Modality), 
                               alpha = 0.5) + theme_void()
  
  if (printValues) {
    g_legend <- function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    legend <- g_legend(p)
    
    
    grid.newpage()
    vp1 <- viewport(width = 0.75, height = 1, x = 0.375, y = .5)
    vpleg <- viewport(width = 0.25, height = 0.5, x = 0.85, y = 0.75)
    subvp <- viewport(width = 0.3, height = 0.3, x = 0.85, y = 0.25)
    print(p + theme(legend.position = "none"), vp = vp1)
    upViewport(0)
    pushViewport(vpleg)
    grid.draw(legend)
    #Make the new viewport active and draw
    upViewport(0)
    pushViewport(subvp)
    grid.draw(tableGrob(table_value[, c(1, 2)], rows = NULL))
  }
  else{return(p)}
}
################################################################################


################################################################################
#' Make Krona files 
#' @aliases merge_krona
#' @param physeq (Required): a \code{\link{phyloseq-class}} object.
#' @param file: the location of the html file to save
#' @param nbSeq (logical, default to TRUE): If true, Krona set the distribution 
#' of sequences in the taxonomy. If False, Krona set the distribution of OTUs 
#' in the taxonomy.
#' @param ranks: Number of the taxonomic ranks to plot (num of the column in tax_table 
#' of your physeq object). Default setting plot all the ranks (argument 'All').
#' @param add_unassigned_rank (number; default = 0). Add unassigned for rank 
#' inferior to 'add_unassigned_rank' when necessary
#' @param name = A name for intermediary files, usefull to name your krona 
#' dataset when merge using merge_krona
#' 
#' @examples
#' data("GlobalPatterns")
#'\dontrun{
#' GlobalPatterns_Acidobacter <- subset_taxa(GlobalPatterns, Phylum=="Acidobacteria")
#' krona(GlobalPatterns_Acidobacter, "Number.of.sequences.html")
#' krona(GlobalPatterns_Acidobacter, "Number.of.OTUs.html", nbSeq = F)
#' merge_krona(c("Number.of.sequences.html", "Number.of.OTUs.html"))
#'}
#'
#' @return A \code{\link{html}} file
#'
#' @author Adrien Taudiere
krona <- function(physeq, file = "krona.html", nbSeq = TRUE, ranks = "All", add_unassigned_rank = 0, name = NULL){
  
  df <- data.frame(unclass(physeq@tax_table[,ranks]))
  df$OTUs <- rownames(physeq@tax_table)
  
  if (ranks[1] == "All") {
    ranks <- seq_along(physeq@tax_table[1,])
  }
  
  if (is.null(name)) {
    if (nbSeq) {name <- "Number.of.sequences"}
    else {name <- "Number.of.OTUs"}
  }
  
  if (nbSeq) {
    df$nb_seq <- taxa_sums(physeq)
  } else {
    df$nb_Otu <- rep(1, length(taxa_sums(physeq)))
  }
  
  df <- df[c(ncol(df), 2:ncol(df) - 1)]             
  res <- lapply(split(df, seq_along(physeq@tax_table[,1])), function(x) as.vector(as.matrix(x))[!is.na(unlist(x))])
  
  res <- lapply(res, function(x) if (length(x) < add_unassigned_rank) {
    x <- c(x, "unassigned")[c(1:length(x) - 1 , length(x) + 1, length(x))]} else {x})

  interm.txt <- paste(tempdir(),"/", name, ".html", sep = "")
  #tempfile(pattern = "file", tmpdir = tempdir(), fileext = "")
  
  lapply(res, cat, "\n", file = interm.txt, append = TRUE, sep = "\t")
  
  cmd <- paste("ktImportText ", interm.txt, " -o ", file, sep = "")
  system(command = cmd)
  system(command = paste("rm", interm.txt))
}

merge_krona <- function(files = NULL, output="mergeKrona.html"){
  
  cmd <- paste("ktImportKrona ", paste(files, collapse = " "), " -o ", output, sep = "")
  
  system(command = cmd)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}