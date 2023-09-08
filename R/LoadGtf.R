#' @title LoadGtf
#'
#' @description Use to import the Ensembl/10x Genomics default genome annotation
#' or other desired genome annotation file in GTF format for optimization for scRNA-seq
#' analysis. Note: This file can be downloaded from 10x Genomics provided reference
#' transcriptome "gene" folder at
#' "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"
#' or Ensembl.org if wish to customize more.
#'
#' @param unoptimized_annotation_path Path to the unoptimized genome annotion GTF file.
#'
#' @return Resulting object contains the genome annotation entries from the genome annotation GTF file.
#' @export
#'
#' @examples
#' LoadGtf(unoptimized_annotation_path = "test_genes.gtf")
LoadGtf <- function(unoptimized_annotation_path){

  #Access test data
  if(unoptimized_annotation_path == "test_genes.gtf"){
    unoptimized_annotation_path <- system.file("extdata", "test_genes.gtf", package = "ReferenceEnhancer")
  }

  genome_annotation <- rtracklayer::import(con = unoptimized_annotation_path, format = "gtf") # Import the original exonic genome annotation file
  genome_annotation <- as.data.frame(genome_annotation)
  return(genome_annotation)
}
