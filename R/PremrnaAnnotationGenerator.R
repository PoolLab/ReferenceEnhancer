#' @title PremrnaAnnotationGenerator
#'
#' @description It supplements original normal gene annotation entries by
#' traditional pre-mRNA entries where transcripts have been redefined as exons
#' and map in the --include-introns mode to retrieve most of available intronic reads.
#'
#' @param genome_annotation ENSEMBL/10x Genomics default genome annotation file (.gtf).
#'
#' @return Generates a basic pre-mRNA reference and saves in working directory as premrna.gtf
#' @export
#'
#' @examples
#' genome_annotation <- LoadGtf("test_genes.gtf")
#' PremrnaAnnotationGenerator(genome_annotation)
PremrnaAnnotationGenerator <- function(genome_annotation){

  exonic_df <- genome_annotation
  premrna_df = exonic_df[exonic_df$type == "transcript",] # Extract all "transcript" entries in the genome annotation to a new variable
  premrna_df$feature = rep("exon", nrow(premrna_df)) # Rename all "feature"

  premrna_df = GenomicRanges::makeGRangesFromDataFrame(premrna_df, keep.extra.columns=TRUE)
  write_gtf(premrna_df, "premrna.gtf")
  print("Pre-mRNA reference has been saved in working directory as premrna.gtf")

}
