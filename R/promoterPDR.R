#' @title Calculate promoter PDR
#' @description Promoter proportion of discordant reads (PDR) was defined as the average PDR of all the CpGs inside.
#' @param file_name A tab-separated values input file. The input file contains details of single CpG PDR with
#'     following columns: chr, start, strand, ConMethReadCount, ConUMethReadCount, DisReadCount, NAReadCount.
#'     Details of PDR were described in Landau \emph{et al.}, \emph{Cancer Cell}, 2014.
#'     extdata/pdrCall_from_Bismark.py can be used to call PDR of single CpG from
#'     Bismark (Krueger \emph{et al.}, \emph{Bioinformatics}, 2011) outputs (files starting with CpG_OB or CpG_OT).
#' @param pro A GRanges object containing promoter annotation.
#' @param min_cpgs An integer value defining the minimum number of available CpGs for each promoter.
#' @return A data frame summarizing promoter PDR levels.
#'     \itemize{
#'       \item Hugo: Hugo symbol
#'       \item PDR: PDR level
#'     }
#' @examples
#' promoterPDR(file_name = system.file("extdata", "PDR.SRR2069925.txt", package = "MethSig"),
#'             pro = makeHG19Promoters())
#'
promoterPDR <- function(file_name, pro, min_cpgs = 3) {
  data <- read.table(file_name, header=T, sep='\t')
  data$totalReadCount <- data$ConMethReadCount + data$ConUMethReadCount + data$DisReadCount
  data <- data[data$totalReadCount >= 10,]
  data$pdr <- data$DisReadCount / data$totalReadCount
  anno <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(data$chr),
                                 ranges=IRanges::IRanges(start=data$start, end=data$start))
  index <- as.list(GenomicRanges::findOverlaps(pro, anno, ignore.strand=T))
  pdrCal <- function(x) {
    if (length(x) < min_cpgs) {
      return(c(NA, NA))
    } else {
      return(c(length(x), mean(data[x, 'pdr'])))
    }
  }
  pdr <- as.data.frame(do.call(rbind, lapply(index, pdrCal)))
  colnames(pdr) <- c('num', 'pdrtumor')
  pdr$Hugo <- names(pro)
  pdr <- na.omit(pdr)
  pdr$num <- NULL
  colnames(pdr) <- c('PDR', 'Hugo')
  pdr <- pdr[,c('Hugo', 'PDR')]
  return(pdr)
}
