#' @title IdentifyOverlappers
#'
#' @description Identifies all overlapping genes based on the ENSEMBL/10x Genomics
#' default genome annotation file (GTF), rank-order them according to the number of gene overlaps.
#' Prioritize this gene list for manual curation focusing on exonically overlapping genes.
#' Saves the list of overlapping genes in working directory (overlapping_gene_list.csv).
#'
#' @param genome_annotation Genome annotation file in .gtf format.
#'
#' @return Rank-ordered gene list of same-strand overlapping genes (gene_overlaps).
#' @export
#'
#' @examples
#' genome_annotation <- LoadGtf("test_genes.gtf")
#' IdentifyOverlappers(genome_annotation)
IdentifyOverlappers <- function(genome_annotation){

  genes_df = genome_annotation[genome_annotation$type == "gene",1:13] # Extract all "gene" entries in the genome annotation to a new variable
  row.names(genes_df) = 1:nrow(genes_df)
  gene_names = genes_df$gene_name
  genes_df = GenomicRanges::makeGRangesFromDataFrame(genes_df, keep.extra.columns=T) # convert into granges object


  overlapper = rep(FALSE, length(gene_names))
  number_of_overlaps = rep(0, length(gene_names))
  overlapping_genes = rep("", length(gene_names))

  for (i in 1:length(gene_names)){
    a = sum(GenomicRanges::countOverlaps(genes_df, genes_df[i]))
    if (a>1){
      overlapper[i] = TRUE
      number_of_overlaps[i] = a-1
      conflict_genes = gene_names[as.logical(GenomicRanges::countOverlaps(genes_df, genes_df[i]))]
      conflict_genes = setdiff(conflict_genes, gene_names[i])
      overlapping_genes[i] = paste(conflict_genes, collapse = ', ')
    }
  }

  overlapping_gene_list = as.data.frame(cbind(gene_names, number_of_overlaps, overlapping_genes))[overlapper,]
  colnames(overlapping_gene_list) = c("gene", "number_of_gene_overlaps", "overlapping_genes")
  overlapping_gene_list$number_of_gene_overlaps = as.integer(overlapping_gene_list$number_of_gene_overlaps)

  o = order(overlapping_gene_list$number_of_gene_overlaps, decreasing = TRUE) # Rank order genes by the number of gene overlaps
  overlapping_gene_list = overlapping_gene_list[o,]

  if(dim(overlapping_gene_list)[1] > 0){
    row.names(overlapping_gene_list) = 1:nrow(overlapping_gene_list)
  }

  overlapping_gene_list["automatic_classification"] <- ""
  overlapping_gene_list["final_classification"] <- ""
  overlapping_gene_list["transcripts_for_deletion"] <- ""

  write_csv(overlapping_gene_list, 'overlapping_gene_list.csv')
  print("A list of overlapping genes has been saved in your working directory (overlapping_gene_list.csv) for manual curation.")
  return(overlapping_gene_list)

}
