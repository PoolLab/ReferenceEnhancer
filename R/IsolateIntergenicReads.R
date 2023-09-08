#' @title IsolateIntergenicReads
#'
#' @description Intergenic reads are extracted from Cell Ranger aligned bam file.
#' Use a scRNA-seq dataset of interest that has been aligned to the unoptimized
#' genome reference with the Cell Ranger count pipeline. Intergenic reads can be
#' identified by two features: their read identity tag RE = "I" (for intergenic)
#' OR their RE=E (for exonic) with AN = <some gene>. The latter reads are in fact
#' intergenic reads since Cell Ranger wrongly classifies reads mapping antisense
#' to an exon as exonic (i.e. RE="E"). The false exonic reads can be recognized
#' and captured as proper intergenic reads by extracting two kinds of reads
#' (RE=I and RE=E & AN=<something else than NA). Also, removing duplicates command
#' in GenomicAlignments package does not work for intergenic (nor for intronic) reads.
#' Duplicate and corrupt read removal has to be done manually (i.e. make sure cellular
#' and molecular barcodes have specified lengths and duplicate barcodes removed).
#'
#' Note that bam files can often be many tens of gigabytes and thus this step is
#' highly memory intensive.
#'
#' @param bam_file_name Path to Cell Ranger generated bam file (run Cell Ranger
#' count pipeline on sequencing data of interest and aligning it to the unoptimized
#' transcriptomic reference).
#' @param index_file_name Path to Cell Ranger generated bam.bai file (run Cell Ranger
#' count pipeline on sequencing data of interest and aligning it to the unoptimized
#' transcriptomic reference).
#'
#' @param barcode_length Optional. Specifies the length of barcode needed. If not specified, defaults to 26.
#'
#' @return Saves extracted intergenic reads as a separate file (“intergenic_reads.bed”)
#' @export
#'
#' @examples
#' IsolateIntergenicReads(
#' bam_file_name = "test_bam.bam",
#' index_file_name = "test_index.bam.bai")
IsolateIntergenicReads <- function(bam_file_name, index_file_name, barcode_length = NULL){

  bamfile = bam_file_name
  indexfile = index_file_name

  if(bamfile == "test_bam.bam" & indexfile == "test_index.bam.bai"){
    bamfile <- system.file("extdata", "test_bam.bam", package = "ReferenceEnhancer")
    indexfile <- system.file("extdata", "test_index.bam.bai", package = "ReferenceEnhancer")
  }

  seq_data = GenomicAlignments::readGAlignments(bamfile, index=indexfile, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE), tag = c("GN", "RE", "CB", "UB", "AN"), what = "flag", tagFilter = list("RE"=c("I", "E"))))
  seq_data = data.frame(seq_data)

  ## Keep only intergenic reads by removing all true exonic reads (i.e. remove exonic reads that lack antisense gene mapping: AN tag =NA)
  intergenic_reads = seq_data$RE=="I"
  false_exonic_reads = !is.na(seq_data$AN)
  all_intergenic_reads = as.logical(intergenic_reads + false_exonic_reads)
  seq_data = seq_data[all_intergenic_reads,] # remaining dataframe contains only intergenic reads. Note, that we are assuming that all false exonic reads are intergenic, which slightly overestimates intergenic read count. This is since some will likely also end up being intronic.

  ## Remove all duplicate reads and reads with corrupt barcodes (i.e. keep reads with 16 nucleotide cellular barcodes and 10 nucleotide molecular barcodes). Note that duplicate removal is required since Cell Ranger does not automatically flag duplicates for intronically and intergenically classified reads.
  seq_data$CB = stringr::str_sub(seq_data$CB, end=-3) # Remove last two elements of the cell barcode. This is an artifact ("-1") added by Cell Ranger software.
  seq_data$barcodes = paste(seq_data$CB, seq_data$UB, sep="") # Assemble the cell barcode / molecular barcode list. Each read included in the gene_cell matrix will have a unique index comprised of the two.

  if(is.null(barcode_length)){
    a = nchar(seq_data$barcodes)==26 # logical vector for selecting reads with non-corrupt barcodes
  }
  else{
    a = nchar(seq_data$barcodes)==barcode_length
  }

  seq_data = seq_data[a,] # exclude all reads that don't have an intact full cellular and molecular barcodes
  length(unique(seq_data$barcodes)) # Determine # of unique intergenic reads
  seq_data = seq_data[!duplicated(seq_data$barcodes),] # exclude all duplicated intergenic reads

  ## Save extracted intergenic reads as a separate file
  gr_seq_data = GenomicRanges::makeGRangesFromDataFrame(seq_data) # coerce to granges object as that makes it possible to save it as a bedfile that bedtools can parse

  ga_seq_data = as(gr_seq_data, "GAlignments") # Coerces GRanges object into a GAlignments object, that can be saved as a bed file. Required for bedtools to link reads to closest 3' gene end.
  m = rtracklayer::asBED(ga_seq_data)# converts GAlignments object into the bed format
  write_bed(m, "intergenic_reads.bed")
  print("Extracted intergenic reads have been saved to your working directory as intergenic_reads.bed.")

}
