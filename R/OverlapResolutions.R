#' @title  OverlapResolutions
#'
#' @description Based on original genome annotation .gtf file and a list of overlapping genes,
#' generates recommended actions for overlapping genes.
#'
#' @param overlap_data A list of overlapping genes generated from IdentifyOverlappers.
#' @param genome_annotation ENSEMBL/10x Genomics default genome annotation file (GTF).
#'
#' @return overlapping_gene_list.csv with added recommendations.
#' @export
#'
#' @examples
#' genome_annotation <- LoadGtf("test_genes.gtf")
#' gene_overlaps <- IdentifyOverlappers(genome_annotation)
#' OverlapResolutions(genome_annotation, gene_overlaps, gene_pattern)
OverlapResolutions <- function(genome_annotation, overlap_data, gene_pattern){
  gene_list <- unique(overlap_data$gene)
  gene_address <- rep(0, length(gene_list))
  for (i in 1:length(gene_list)){
    gene_address[i] <- which(overlap_data$gene==gene_list[i])[1]
  }

  overlap_data <- overlap_data[gene_address,]

  rownames(overlap_data) <- overlap_data[,'gene']

  overlap_data['automatic_classification'] <- NA

  for(key in (rownames(overlap_data))){

    # Check that the gene is not classified already
    if(is.na(overlap_data[key,'automatic_classification'])){
      gene_A <- subset(genome_annotation, gene_name == key)

      if(overlap_data[key,'number_of_gene_overlaps'] > 1){
        overlaps <- as.list(strsplit(overlap_data[key,'overlapping_genes'], ", "))

        for(item in overlaps[[1]]){
          gene_B = genome_annotation[genome_annotation['gene_name'] == item,]

          gene_A_exons = return_exons(gene_A)
          gene_B_exons = return_exons(gene_B)

          if(exon_overlap(gene_A_exons, gene_B_exons) == TRUE){

            overlap_data[item, 'automatic_classification'] = 'Manual inspection'

            if(is.na(overlap_data[key,'automatic_classification']) | overlap_data[key,'automatic_classification'] != 'Manual inspection'){
              overlap_data[key,'automatic_classification'] = 'Manual inspection'
            }
          }
          else{

            if(is.na(overlap_data[key,'automatic_classification'])){
              overlap_data[key,'automatic_classification'] = 'Keep as is'

              if(overlap_data[item,'number_of_gene_overlaps'] > 1){
                overlap_data[item,'automatic_classification'] = 'Manual inspection'
              }
              else{
                overlap_data[item,'automatic_classification'] = 'Keep as is'
              }
            }
          }
        }
      }

      if(overlap_data[key,'number_of_gene_overlaps'] == 1){
        overlapping <- overlap_data[key,'overlapping_genes'][[1]]
        gene_B <- subset(genome_annotation, gene_name == overlapping)
        strand <- gene_A[1,'strand']

        gene_A_exons = return_exons(gene_A)
        gene_B_exons = return_exons(gene_B)

        # Check if both - key and overlapping gene - are pseudogenes
        if(both_pseudo(key, overlapping, gene_pattern) == TRUE){
          overlap_data[key, 'automatic_classification'] = 'Manual inspection'
          overlap_data[overlapping[[1]], 'automatic_classification'] = 'Manual inspection'
        }

        # Check for pseudogene
        else if(pseudo_overlap(key, overlapping, gene_A_exons, gene_B_exons, gene_pattern) == key){
          overlap_data[key, 'automatic_classification'] = 'Delete'
          overlap_data[overlapping[[1]], 'automatic_classification'] = 'Keep as is'
        }

        else if(pseudo_overlap(key, overlapping, gene_A_exons, gene_B_exons, gene_pattern) == overlapping){
          overlap_data[key, 'automatic_classification'] = 'Keep as is'
          overlap_data[overlapping[[1]], 'automatic_classification'] = 'Delete'
        }

        else if(pseudo_overlap(key, overlapping, gene_A_exons, gene_B_exons, gene_pattern) == 'exonic'){
          overlap_data[key, 'automatic_classification'] = 'Keep as is'
          overlap_data[overlapping[[1]], 'automatic_classification'] = 'Keep as is'
        }

        # Check for readthrough
        else if(exon_overlap(gene_A_exons, gene_B_exons) == TRUE){
          if(strand == '+'){
            name_A = key
            name_B = overlapping
            result = readthrough_or_premature_plus(name_A, gene_A, name_B, gene_B, gene_A_exons, gene_B_exons)

            if(result[[3]] == 'readthrough'){
              overlap_data[result[[1]],'automatic_classification'] = 'Readthrough transcript deletion'
              overlap_data[result[[2]],'automatic_classification'] = 'Keep as is'
            }
            else if(result[[3]] == 'premature'){
              overlap_data[result[[1]],'automatic_classification'] = 'Keep as is'
              overlap_data[result[[2]],'automatic_classification'] = 'Premature transcript deletion'
            }
            else if(result[[3]] == 'manual'){
              overlap_data[result[[1]],'automatic_classification'] = 'Manual inspection'
              overlap_data[result[[2]],'automatic_classification'] = 'Manual inspection'
            }
          }

          else if(strand == '-'){
            name_A = key
            name_B = overlapping
            result = readthrough_or_premature_min(name_A, gene_A, name_B, gene_B, gene_A_exons, gene_B_exons)

            if(result[[3]] == 'readthrough'){
              overlap_data[result[[1]],'automatic_classification'] = 'Readthrough transcript deletion'
              overlap_data[result[[2]],'automatic_classification'] = 'Keep as is'
            }
            else if(result[[3]] == 'premature'){
              overlap_data[result[[1]],'automatic_classification'] = 'Keep as is'
              overlap_data[result[[2]],'automatic_classification'] = 'Premature transcript deletion'
            }
            else if(result[[3]] == 'manual'){
              overlap_data[result[[1]],'automatic_classification'] = 'Manual inspection'
              overlap_data[result[[2]],'automatic_classification'] = 'Manual inspection'
            }
          }
        }

        else if(exon_overlap(gene_A_exons, gene_B_exons) == FALSE){
          overlap_data[key,'automatic_classification'] = 'Keep as is'
          overlap_data[overlapping,'automatic_classification'] = 'Keep as is'
        }
      }
    }
  }

  print("Overlapping genes list (overlapping_gene_list.csv) has been updated with recommended action categories and the file has been saved in your working directory")
  write_csv(overlap_data, "overlapping_gene_list.csv")

}
