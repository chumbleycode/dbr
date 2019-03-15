#'A gene by transcription factor binding motif matrix.
#'
#' Broadly, this contains binding loci counts for each transcription factor motif in each gene. In particular, it counts the number of times any part of the transcription factor binding motif overlaps with the 1000 base pairs upstream of the transcription start site.
#'
#' @format A gene by transcription factor binding motif matrix.
#' @source biomart
"utr1"

#'A gene by transcription factor binding motif matrix.
#'
#' Broadly, this contains binding loci counts for each transcription factor motif in each gene. In particular, it counts the number of times any part of the transcription factor binding motif overlaps with an exon.
#'
#' @format A gene by transcription factor binding motif matrix.
#'
#' @source biomart
"exonic1"


#'A gene by transcription factor binding motif matrix.
#'
#' Broadly, this contains binding loci counts for each transcription factor motif in each gene. In particular, it counts the number of times any part of the transcription factor binding motif overlaps with an exon AND simultaneously overlaps with the 1000 base pairs upstream of the transcription start site.
#'
#' @format A gene by transcription factor binding motif matrix.
#'
#' @source biomart
"exonic1_utr1"
