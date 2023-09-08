#' @title GenerateExtensionCandidates
#'
#' @description Identifies candidate genes for extension with excess 3' intergenic
#' reads and creates a rank ordered list of genes as a function of 3' intergenic
#' read mapping within 10kb of known gene end. You can use this as a prioritized
#' gene list for gene extension to examine in Integrated Genomics Viewer.
#'
#' Note: It runs partially in Bash/Linux terminal. Make sure bedtools is installed
#' and provide a path in the function if you get an error message.
#'
#' @param bedops_loc Optional. Location of bedtools in file system.
#'
#' @return Rank ordered list of gene extension candidates saved to working directory
#' as “gene_extension_candidates.csv”.
#' @export
#'
#' @examples
#' GenerateExtensionCandidates(bedtools_loc = NULL)
GenerateExtensionCandidates <- function(bedtools_loc = NULL){

  ## In bash/linux terminal: Make sure bedtools is in PATH (make sure bedtools is installed and in the PATH variable in Linux or MacOS)

  system("sort -k 1,1 -k2,2n gene_ranges.bed > gene_ranges1.bed")
  system("sort -k 1,1 -k2,2n intergenic_reads.bed > intergenic_reads1.bed")

  # Checks and adds bedtools to path
  if(is.null(bedtools_loc)){
    if(is.na(unlist(strsplit(system("whereis bedtools", intern = TRUE),": "))[2])){
      print("Didn't find bedtools. Please install bedtools or provide a path to bedtools.")
    }
    else{
      old_path <- Sys.getenv("PATH")
      bedtools_loc = unlist(strsplit(system("whereis bedtools", intern = TRUE),": "))[2]
      Sys.setenv(PATH = paste(old_path, bedtools_loc, sep = ":"))
    }
  }
  else{
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(old_path, bedtools_loc, sep = ":"))
  }


  system("bedtools closest -a intergenic_reads1.bed -b gene_ranges1.bed -s -D a -fu > results.txt") # resulting file contains sequencing reads with distance data from closest 3' gene identity and end

  ## In R: Save a rank ordered list of genes with highest-to-lowest number of intergenic reads within 10kb of its known gene end.

  summary_data = read.table("results.txt", sep = "\t")

  summary_data = summary_data[summary_data$V23>-10000,] # retain only sequencing reads within 10kb of known gene ends. Change to more or less stringent as desired.
  summary_data = summary_data[summary_data$V23<0,] # retain only sequencing reads within 10kb of known gene ends. Change to more or less stringent as desired.

  hist(summary_data$V23) # plot histogram of intergenic sequencing reads as a function of distance from 3' gene ends.

  summary_data_genes = table(summary_data$V22) # Summarizes # of intergenic reads within 10kb of known gene ends for each gene.
  o = order(summary_data_genes, decreasing = TRUE) # Rank order the gene list
  length(summary_data_genes)
  summary_data_genes = summary_data_genes[o]
  length(summary_data_genes[summary_data_genes>10]) # Determine number of genes with more than 10 intergenic reads within 10kb of known gene end
  summary_data_genes = summary_data_genes[summary_data_genes>10] # Threshold gene list based on the amount of intergenic gene loading.
  summary_data_genes = data.frame(summary_data_genes)
  dim(summary_data_genes)
  summary_data_genes[1:40,]
  summary_data_genes["update_start"] <- ""
  summary_data_genes["update_end"] <- ""


  write_csv(summary_data_genes, "gene_extension_candidates.csv")
  print("A rank ordered list of gene extension candidates has been saved to working directory as gene_extension_candidates.csv")

}
