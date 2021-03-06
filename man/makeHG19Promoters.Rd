% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeHG19Promoters.R
\name{makeHG19Promoters}
\alias{makeHG19Promoters}
\title{Generate promoter annotation of hg19 RefSeq genes}
\usage{
makeHG19Promoters(up = 2000, down = 2000)
}
\arguments{
\item{up}{An integer value defining the number of nucleotides towards the 5' end relative to
the transcription start site.}

\item{down}{An integer value defining the number of nucleotides towards the 3' end relative to
the transcription start site.}
}
\value{
A GRanges object containing promoter annotation.
}
\description{
Generate promoter annotation of hg19 RefSeq genes
}
\examples{
makeHG19Promoters(up = 2000, down = 2000)

}
