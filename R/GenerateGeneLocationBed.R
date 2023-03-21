#' @title GenerateGeneLocationBed
#'
#' @description Makes a bed file with gene boundaries, which is required for
#' assigning intergenic reads to a specific gene.
#'
#' Note: This step is partially run in linux Terminal in Bash and requires
#' bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html).
#' Put bedtools in PATH after installing.
#'
#' In linux terminal, navigate to folder with the genome annotation of interest.
#' Assuming that it is named "genes.gtf" per 10x Genomics convention.
#'
#' @param genome_annotation Genome annotation file in .gtf format.
#'
#' @return Saves gene_ranges.bed in working directory.
#' @export
#'
#' @examples
#' genome_annotation <- LoadGtf("test_genes.gtf")
#' GenerateGeneLocationBed(genome_annotation, bedops_loc)
GenerateGeneLocationBed <- function(genome_annotation, bedops_loc = NULL){
  gene_ranges_df <- genome_annotation
  gene_ranges_df <- gene_ranges_df[gene_ranges_df$type == "gene",] # Extract all "gene" entries in the genome annotation ot a new variable
  gene_ranges_df <- GenomicRanges::makeGRangesFromDataFrame(gene_ranges_df, keep.extra.columns=TRUE)
  rtracklayer::export(gene_ranges_df, "gene_ranges.gtf", format = "gtf")

  ## Add "transcript_id """ column to the gtf file to make it compatible with bedtools format (through terminal)
  system('awk \'{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }\' gene_ranges.gtf > gene_ranges1.gtf')

  ## Convert reference gtf into bed with bedops (make sure bedops is in the PATH variable for linux or MAcOS)
  if(length(bedops_loc) != 0){
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(old_path, bedops_loc, sep = ":"))
  }

  else{
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(old_path, "/usr/bin/bedops", sep = ":"))
  }

  system('gtf2bed < gene_ranges1.gtf > gene_ranges.bed') # Creates a bed file with gene boundaries

  ## The following code in R replaces final column with gene name. Make sure you navigate to same folder in R.

  gene_ranges = read.table("gene_ranges.bed", sep = "\t")

  if(dim(gene_ranges)[1] > 0){
    for (i in 1:dim(gene_ranges)[1])
    {
      a = gene_ranges[i,10]
      res <- stringr::str_match(a, "gene_name\\s*(.*?)\\s*;")
      b = res[,2]
      gene_ranges[i, 10] = b
    }
  }

  ## Remove gene_ranges.gtf
  file.remove("./gene_ranges.gtf")
  file.remove("./gene_ranges1.gtf")

  ## Save outcome
  write.table(gene_ranges, "gene_ranges.bed", sep="\t",row.names=FALSE, col.names=FALSE, quote = FALSE)
  print("Gene ranges file (gene_ranges.bed) has been saved in working directory.")

}
