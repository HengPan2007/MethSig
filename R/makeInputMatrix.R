#' @title Make input matrix
#' @description Add DHcR, PDR and sequencing depth information of samples to be infered into input matrix.
#' @param names_list A list of sample names.
#' @param matCV A matrix of covarites which are included into analysis. Ideally, users can include any covariate
#'     that will be correlated with promoter DHcR in tumor samples and a hugo symbol column must be included.
#'     An example can be loaded by using invisible(matCV). This covarite matrix contains following columns:
#'     Hugo, DHcR_normal, PDR_normal, GEXP_Normal, Reptime.
#' @param pro A GRanges object containing promoter annotation.
#' @param input_dir A string indicates where DHcR and PDR files of all the samples are stored.
#'     DHcR files should be named as DMC.sample_name.txt format.
#'     PDR files should be named as PDR.sample_name.txt format.
#' @param min_cpgs A list of two integer value >= 0L. (Minimum number of CpGs requied by promoterDHcR,
#'     Minimum number of CpGs required by promoterPDR).
#' @return A data frame used as input for hypermethylation inference.
#'     \itemize{
#'       \item Id: Sample Id
#'       \item Hugo: Hugo symbol
#'       \item Covariates included in matCV
#'       \item PDR_Tumor: PDR level in tumor samples
#'       \item Depth_Tumor: Average sequencing depth in tumor samples
#'       \item CpGs_Tumor: Number of covered CpGs in tumor samples
#'       \item DHcR_Tumor: DHcR level in tumor samples
#'     }
#' @examples
#' makeInputMatrix(names_list, matCV, pro, input.dir)
#'
makeInputMatrix <- function(names_list, matCV, pro, input_dir, min_cpgs=c(5, 3)) {
  col_name <- colnames(matCV)
  makeMatrixPerSample <- function(name) {
    dhcr <- promoterDHcR(paste0(input_dir, '/DMC.', name, '.txt'), pro, min_cpgs[1])
    colnames(dhcr)[2:5] <- paste0(colnames(dhcr)[2:5], '_Tumor')
    input.mat <- merge(matCV, dhcr)

    pdr <- promoterPDR(paste0(input_dir, 'PDR.', name, '.txt'), pro, min_cpgs[2])
    colnames(pdr)[2] <- paste0(colnames(pdr)[2], '_Tumor')

    input.mat <- merge(input.mat, pdr)
    input.mat$Id <- name
    return(input.mat)
  }
  foo <- do.call(rbind, lapply(names_list, makeMatrixPerSample))
  foo <- foo[,c('Id', col_name, 'PDR_Tumor', 'Depth_Tumor', 'CpGs_Tumor', 'DHcR_Tumor')]
  return(foo)
}
