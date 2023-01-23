readthrough_or_premature_plus <- function(name_A, gene_A, name_B, gene_B, gene_A_exons, gene_B_exons){
  select_gene_A <- gene_A[gene_A$type=='gene', ]
  gene_A_start <-  select_gene_A[,'start']
  gene_A_end <- select_gene_A[,'end']
  select_gene_B <- gene_B[gene_B$type=='gene', ]
  gene_B_start <-  select_gene_B[,'start']
  gene_B_end <- select_gene_B[,'end']

  if(gene_A_start < gene_B_start){
    upstream_name <- name_A
    upstream <- gene_A
    downstream_name <- name_B
    downstream <- gene_B
  }

  else if(gene_B_start < gene_A_start){
    upstream_name <- name_B
    upstream <- gene_B
    downstream_name <- name_A
    downstream <- gene_A
  }

  else if(gene_A_end < gene_B_end){
    upstream_name <- name_A
    upstream <- gene_A
    downstream_name <- name_B
    downstream <- gene_B
  }

  else{
    upstream_name <- name_B
    upstream <- gene_B
    downstream_name <- name_A
    downstream <- gene_A
  }

  upstream_trx <- list(upstream[upstream$type == 'transcript',][,'start'], upstream[upstream$type == 'transcript',][,'end'])
  downstream_trx = list(downstream[downstream$type == 'transcript',][,'start'], downstream[downstream$type == 'transcript',][,'end'])

  return(readthrough_or_premature(upstream_name, downstream_name, upstream_trx, downstream_trx))
}
