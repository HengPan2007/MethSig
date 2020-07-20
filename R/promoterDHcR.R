#' @title Calculate promoter DHcR
#' @description Promoter differentially hypermethylated cytosine ratio (DHcR) was defined as the ratio of
#'     hypermethylated cytosines (HCs) to the total number of promoter CpGs profiled.
#' @details HCs of each sample were defined as CpGs at which DNAme is statistically higher than
#'     the average DNAme of control samples (FDR = 20\%, Chi-squared test). Only CpGs with read depth
#'     greater than 10 reads were included in the analysis.
#' @param file_name A tab-separated values input file. The input file contains details of differentially methylated
#'     cytosines with following columns (V1 to V11): chr, pos, numC in control, numC+numT in control,
#'     numC in tumor, numC + numT in tumor, CpG methylation ratio (tumor methylation / control methylation),
#'     chi-squared p-value, adjusted p-value, significance, hyper or hypo in tumor.
#'     Details of generating this type of files were described in Pan et al, Cancer Systems Biology, 2018.
#' @param pro A GRanges object containing promoter annotation.
#' @param min_cpgs Single integer value >= 0L. The minimum number of covered CpGs for each promoter.
#' @return A data frame summarizing promoter DHcR levels.
#'     \itemize{
#'       \item Hugo: Hugo symbol
#'       \item Depth: Average sequencing depth
#'       \item HCs: Number of HCs
#'       \item CpGs: Number of covered CpGs
#'       \item DHcR: DHcR level
#'     }
#' @examples
#' promoterDHcR('DMC.SRR2069925.txt', pro)
#'
promoterDHcR <- function(file_name, pro, min_cpgs=5) {
  example <- read.table(file_name, sep='\t')
  anno <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(example$V1),
                                 ranges=IRanges::IRanges(start=example$V2, end=example$V2))
  index <- as.list(GenomicRanges::findOverlaps(pro, anno, ignore.strand=T))
  tumor <- as.data.frame(t(sapply(index, function(x) c(table(example[x, 'V11']), mean(example[x, 'V6'])))))
  colnames(tumor)[4] <- 'cov'
  tumor$total <- tumor$UP + tumor$DOWN + tumor$'-'
  tumor$hratio <- tumor$UP/tumor$total
  tumor$Hugo <- names(pro)
  tumor$'-' <- NULL
  tumor$DOWN <- NULL
  tumor <- tumor[tumor$total >= min_cpgs, ]
  colnames(tumor)[1:5] <- c('HCs', 'Depth', 'CpGs', 'DHcR', 'Hugo')
  tumor <- tumor[, c('Hugo', 'Depth', 'HCs', 'CpGs', 'DHcR')]
  return(tumor)
}
