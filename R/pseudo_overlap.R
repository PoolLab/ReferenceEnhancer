pseudo_overlap <- function(key, overlapping, gene_A_exons, gene_B_exons, gene_pattern){

  if(missing(gene_pattern)){
    return('empty')
  }

  else{
    # Check for exon overlap
    if(sum(stringr::str_detect(key, gene_pattern)) > 0 | sum(stringr::str_detect(overlapping, gene_pattern)) > 0 ){
      if(exon_overlap(gene_A_exons, gene_B_exons) == TRUE){
        # Check if gene_A is a pseudogene
        if(sum(stringr::str_detect(key, gene_pattern)) > 0){
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

}
