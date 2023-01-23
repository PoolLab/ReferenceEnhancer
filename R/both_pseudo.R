both_pseudo <- function(key, overlapping){
  return((stringr::str_sub(key, 1, 2) == 'Gm' | stringr::str_sub(key, - 3, - 1) == 'Rik') & (stringr::str_sub(overlapping, 1, 2) == 'Gm' | stringr::str_sub(overlapping, - 3, - 1) == 'Rik'))
}
