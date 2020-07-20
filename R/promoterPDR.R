#' @title Calculate promoter PDR
#' @description Promoter proportion of discordant reads (PDR) was defined as the average PDR of all the CpGs inside.
#' @details If all the CpGs on a specific read are methylated, or all of the CpGs on a read are unmethylated,
#'     the read is classified as concordant; otherwise it is classified as discordant. At each CpG,
#'     the PDR is equal to the number of discordant reads that cover that location divided by
#'     the total number of read that cover that location. Only CpGs with read depth greater than 10 reads
#'     and covered by reads that contain at least 4 CpGs were included into the analysis.
#' @param file_name A tab-separated values input file. The input file contains PDR of single CpG with
#'     following columns: chr, start, strand, ConMethReadCount, ConUMethReadCount, DisReadCount, NAReadCount.
#'     Details of PDR were described in Landau et al., Cancer Cell, 2014. extdata/pdrCall_from_Bismark.py
#'     can be used to call PDR of single CpG from Bismark outputs (files starting with CpG_OB or CpG OT).
#' @param pro A GRanges object containing promoter annotation.
#' @param min_cpgs Single integer value >= 0L. The minimum number available CpGs for each promoter.
#' @return A data frame summarizing promoter PDR levels.
#'     \itemize{
#'       \item Hugo: Hugo symbol
#'       \item PDR: PDR level
#'     }
#' @examples
#' promoterDHcR('PDR.SRR2069925.txt', pro)
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
