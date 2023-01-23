readthrough_or_premature <- function(upstream_name, downstream_name, upstream_trx, downstream_trx){
  max_u = 0
  for(trx_u in 1:length(upstream_trx[[1]])){
    count_u = 0
    for (trx_d in 1:length(downstream_trx[[1]])){
      x = seq(downstream_trx[[1]][trx_d], downstream_trx[[2]][trx_d]-1)
      y = seq(upstream_trx[[1]][trx_u], upstream_trx[[2]][trx_u]-1)

      if(length(intersect(x,y)) > 0){
        count_u = count_u + 1
      }

      if(count_u > max_u){
        max_u = count_u
      }
    }
  }

  max_d = 0
  for(trx_d in 1:length(downstream_trx[[1]])){
    count_d = 0
    for(trx_u in 1:length(upstream_trx[[1]])){
      x = seq(downstream_trx[[1]][trx_d], downstream_trx[[2]][trx_d]-1)
      y = seq(upstream_trx[[1]][trx_u], upstream_trx[[2]][trx_u]-1)
      if(length(intersect(x,y)) > 0){
        count_d = count_d + 1
      }
    }
    if(count_d > max_d){
      max_d = count_d
    }
  }

  if(max_u > max_d){
    result <- list(upstream_name, downstream_name, "readthrough")
    return(result)
  }
  else if(max_d > max_u){
    result <- list(upstream_name, downstream_name, "premature")
    return(result)
  }
  else if(min(length(upstream_trx), length(downstream_trx)) == 1 & max(length(upstream_trx), length(downstream_trx)) != 1){
    print("MAX & MIN")
    result <- list(upstream_name, downstream_name, "manual")
    return(result)
  }
  else if(max_u == 1 & max_d == 1){
    result <- list(upstream_name, downstream_name, "manual")
    return(result)
  }
  else if(max_u == max_d){
    if(length(upstream_trx[[1]]) > length(downstream_trx[[1]])){
      result <- list(upstream_name, downstream_name, "readthrough")
      return(result)
    }
    else{
      result <- list(upstream_name, downstream_name, "premature")
      return(result)
    }
  }
  else{
    result <- list(upstream_name, downstream_name, "manual")
    return(result)
  }
}
