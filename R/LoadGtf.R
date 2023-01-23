#' @title LoadGtf
#'
#' @description Use to import the ENSEMBL/10x Genomics default genome annotation
#' file (GTF).
#' Note: This file can be downloaded from 10x Genomics provided reference
#' transcriptome "gene" folder at "https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"
#' or Ensembl.org if wish to customize more.
#'
#' @param genes_gtf_path Path to annotion .gtf file.
#'
#' @return Genome annotation file (genome_annotation).
#' @export
#'
#' @examples
#' LoadGtf("test_genes.gtf")
LoadGtf <- function(genes_gtf_path){

  #Access test data
  if(genes_gtf_path == "test_genes.gtf"){
    genes_gtf_path <- system.file("extdata", "test_genes.gtf", package = "ReferenceEnhancer")
  }

  genome_annotation <- rtracklayer::import(con = genes_gtf_path, format = "gtf") # Import the original exonic genome annotation file
  genome_annotation <- as.data.frame(genome_annotation)
  return(genome_annotation)
}
