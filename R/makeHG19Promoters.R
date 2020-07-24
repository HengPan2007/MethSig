#' Generate promoter annotation of hg19 refseq genes
#'
#' @param up An integer value defining the number of nucleotides toward the 5' end relative to
#'     the transcription start site.
#' @param down An integer value defining the number of nucleotides toward the 3' end relative to
#'     the transcription start site.
#' @return A GRanges object containing promoter annotation.
#' @examples
#' makeHG19Promoters(2000, 2000)
#' makeHG19Promoters(1000, 1000)

makeHG19Promoters <- function(up=2000, down=2000) {
  genes <- invisible(geneAnnoHG19)
  pro <- GenomicRanges::promoters(genes, upstream=up, downstream=down)
  return(pro)
}
