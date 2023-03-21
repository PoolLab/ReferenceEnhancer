both_pseudo <- function(key, overlapping, gene_pattern){

  if (missing(gene_pattern)){
    return(FALSE)
  }

  else{
    key_pseudo = sum(stringr::str_detect(key, gene_pattern))
    overlapping_pseudo = sum(stringr::str_detect(overlapping, gene_pattern))
    return((key_pseudo + overlapping_pseudo) > 1)

  }

  }
