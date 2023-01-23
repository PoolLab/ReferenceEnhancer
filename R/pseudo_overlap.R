pseudo_overlap <- function(key, overlapping, gene_A_exons, gene_B_exons){
  # Check for exon overlap
  if(stringr::str_sub(key, 1, 2) == 'Gm' | stringr::str_sub(key, - 3, - 1) == 'Rik' | stringr::str_sub(overlapping, 1, 2) == 'Gm' | stringr::str_sub(overlapping, - 3, - 1) == 'Rik'){
    if(exon_overlap(gene_A_exons, gene_B_exons) == TRUE){
      # Check if gene_A is a pseudogene
      if(stringr::str_sub(key, 1, 2) == 'Gm' | stringr::str_sub(key, - 3, - 1) == 'Rik'){
        return(key)
      }
      else{
        return(overlapping)
      }
    }
    else{
      return('exonic')
    }
  }

  else{
    return('empty')
  }
}
