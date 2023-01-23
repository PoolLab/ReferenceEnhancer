return_exons <- function(gene_name){
  exon_subset <- subset(gene_name, type == 'exon')
  return(data.frame(exon_subset['start'], exon_subset['end']))
}
