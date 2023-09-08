#' @title  OptimizedAnnotationAssembler
#'
#' @description
#' OptimizedAnnotationAssembler generates the scRNA-seq optimized genome annotation.
#' The resulting optimized genome annotation can be used to generate the transcriptomic
#' reference for mapping single-cell sequencing data (e.g. with cellranger mkref
#' or STAR --runMode genomeGenerate). Note that completing this step is time intensive
#' and can sometimes take 12-24 hours depending on the length of the annotation
#' to be optimized.
#' This function goes through the following steps:
#' 0. Load data and libraries:
#'    - genome annotation file to be optimized in GTF.
#'    - "overlapping_gene_list.csv" file specifying how to resolve gene overlap
#'    derived issues. "Delete" entries in $final_classification field mark genes
#'    for deletion. Transcript names in $transcripts_for_deletion mark specific
#'    transcripts for deletion.
#'    - "gene_extension_candidates.csv" specifying updated gene boundaries for
#'    incorporating intergenic reads.
#'    - "rename_genes.csv" specifying gene names to be replaced and new names
#'    (under $old_names and $new_names fields, respectively).
#'  1. Resolve "self-overlapping" gene (duplicate gene_ids) derived issues.
#'  Required for making references compatible with multiome workflows.
#'  2. Creates pre-mRNA genome annotation from input genome annotation. This step
#'  extracts all transcript entries from the genome annotation and defines them
#'  as full length exons with new transcript IDs and corresponding transcripts.
#'  This allows to capture many intronically mapped reads that otherwise get discarded.
#'  3. Gene deletion step: Deletes all annotation entries for genes destined for
#'  deletion (has "Delete" entry in $final_classification field of
#'  "overlapping_gene_list.csv".
#'  4. Transcript deletion step: Deletes all transcripts destined for deletion
#'  (transcript names listed in the "transcripts_for_deletion" column in
#'  "overlapping_gene_list.csv".
#'  5. Gene coordinate adjustment step: Replaces the left most or right most
#'  coordinate of the first exon of a gene in genome annotation if there is a
#'  coordinate in columns $new_left or $new_right in the
#'  "gene_extension_candidates.csv".
#'  6. Adds pre-mRNA reads to all genes not in the gene overlap list.
#'  7. Renames genes to avoid discarding expression data with near perfect terminal
#'  exon overlap.
#'  8. Saves the optimized genome annotation in a new GTF file.
#'
#' @param unoptimized_annotation_path path to unoptimized genome annotation file in GTF.
#' @param gene_overlaps overlapping genes list generated with IdentifyOverlappers function.
#' @param gene_extension list of gene extension candidates generated with GenerateExtensionCandidates function.
#' @param gene_replacement manually generated list of gene names to be replaced in .csv format. Column names: old_name, new_name. Optional.
#'
#' @return Single-cell RNA-seq optimized genome annotation that can be used to
#' generate the transcriptomic reference (e.g. with cellranger mkref or
#' STAR --runMode genomeGenerate pipelines) for mapping single-cell sequencing data.
#' @export
#'
#' @examples
#' OptimizedAnnotationAssembler(
#' unoptimized_annotation_path = "test_genes.gtf",
#' gene_overlaps = "test_overlapping_gene_list.csv",
#' gene_extension = "./gene_extension_candidates.csv",
#' gene_replacement = "test_gene_replacement.csv")
OptimizedAnnotationAssembler <- function(unoptimized_annotation_path, gene_overlaps, gene_extension, gene_replacement){

  if(gene_overlaps == "test_overlapping_gene_list.csv"){
    gene_overlaps <- system.file("extdata", "test_overlapping_gene_list.csv", package = "ReferenceEnhancer")
  }


  unoptimized_df <- LoadGtf(unoptimized_annotation_path)

  overlap_df = read.csv(gene_overlaps, header=T)

  new_df = unoptimized_df


  ####  1. Create premRNA genome annotation from input gtf that defines transcripts as exons ####
  ###############################################################################################
  transcripts_df = unoptimized_df[unoptimized_df$type == "transcript",]
  exons_df = transcripts_df # Create new dataframe to contain premrna exons
  exons_df$type = rep("exon", nrow(exons_df)) # rename "type" from transcripts to exon

  premrna_df = gdata::interleave(transcripts_df, exons_df) # interleave transript entries with exon entries
  premrna_df$transcript_id = gsub("000000", "100000", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000001", "110001", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000002", "110002", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000003", "110003", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000004", "110004", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000005", "110005", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000006", "110006", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000007", "110007", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000008", "110008", premrna_df$transcript_id)
  premrna_df$transcript_id = gsub("000009", "110009", premrna_df$transcript_id)

  rm(unoptimized_df)

  ####  2. Delete select genes ####
  #################################
  genes_to_delete = overlap_df$genes[overlap_df$final_classification == "Delete"]
  new_df = new_df[!new_df$gene_name %in% genes_to_delete,]

  ####  3. Delete select transcripts ####
  #######################################
  transcripts_to_delete = overlap_df$transcripts_for_deletion
  transcripts_to_delete <- transcripts_to_delete[transcripts_to_delete!=""]

  transcripts_to_delete_final = transcripts_to_delete[!stringr::str_detect(transcripts_to_delete, ", ")]

  if(length(transcripts_to_delete) != 0){
    for (i in 1:length(transcripts_to_delete)){
      a = transcripts_to_delete[i]
      if (stringr::str_detect(a, ", ")){
        split_elements <- unlist(stringr::str_split(a, ", "))
        transcripts_to_delete_final = c(transcripts_to_delete_final, split_elements)
      }
    }
  }

  transcripts_to_delete = transcripts_to_delete_final

  new_df = new_df[!new_df$transcript_name %in% transcripts_to_delete,]

  ####  4. Adjust gene coordinates ####
  #####################################
  boundary_fix = read.csv(gene_extension, header=T)

  left_genes = as.data.frame(cbind(boundary_fix$genes[!is.na(boundary_fix$update_start)], boundary_fix$update_start[!is.na(boundary_fix$update_start)]))

  colnames(left_genes) = c("genes", "update_start")
  left_genes$update_start = as.numeric(left_genes$update_start)
  right_genes = as.data.frame(cbind(boundary_fix$genes[!is.na(boundary_fix$update_end)], boundary_fix$update_end[!is.na(boundary_fix$update_end)]))
  colnames(right_genes) = c("genes", "update_end")
  right_genes$update_end = as.numeric(right_genes$update_end)

  left_exon_difs = rep(0, length(left_genes)) # for troubleshooting
  right_exon_difs = rep(0, length(right_genes))

  for (i in 1:dim(left_genes)[1]){
    gene_entries = which(new_df$gene_name == left_genes[i, 1])
    type_entries = new_df$type[gene_entries]
    first_gene_exon = head(gene_entries[type_entries == "exon"], 1)
    new_df[first_gene_exon, 2] = left_genes[i, 2]

    if(identical(new_df[first_gene_exon, 3], integer(0)) & identical(new_df[first_gene_exon, 2], integer(0))){

    }
    else{
      left_exon_difs[i] = new_df[first_gene_exon, 3] - new_df[first_gene_exon, 2]
    }


  }

  for (i in 1:dim(right_genes)[1]){
    gene_entries = which(new_df$gene_name == right_genes[i, 1])
    type_entries = new_df$type[gene_entries]
    last_gene_exon = tail(gene_entries[type_entries == "exon"], 1)
    new_df[last_gene_exon, 3] = right_genes[i, 2]

    if(identical(new_df[last_gene_exon, 3], integer(0)) & identical(new_df[last_gene_exon, 2], integer(0))){

    }
    else{
      right_exon_difs[i] = new_df[last_gene_exon, 3] - new_df[last_gene_exon, 2]
    }

  }

  #### 5. Add pre-mRNA transcripts to genes not in the gene overlap list ####
  ############################################################################
  # Explanation: Cellranger --include-introns mode unfortunately does not pick up on many intronic reads (unclear why despite lengthy correspondence with their support). I can pick those up however if I add the pre-mRNA transcripts to respective genes as exons with new transcript_id values.

  ## Genes to modify

  #overlap_df$genes # genes to exclude from premrna reference appending

  genes_to_append = unique(new_df$gene_name)
  genes_to_append = setdiff(genes_to_append, overlap_df$genes)

  ## Give new transcript_ids to everything in the pre-mRNA gtf

  for (i in 1:dim(premrna_df)[1]){
    premrna_df$transcript_id[i] = as.character(i)
  }

  ## Reformat the gtf dataframes such that we can add premrna entries to the original unoptimized entries and thus compile a hybrid reference for capturing intronic reads

  final_colnames = intersect(colnames(new_df), colnames(premrna_df))

  new_df = new_df[, final_colnames]
  premrna_df = premrna_df[, final_colnames]

  ## Append premrna transcript to the end of the gene

  genes_to_append = genes_to_append[1:(length(genes_to_append)-1)]

  for (i in genes_to_append){
    insert = premrna_df[premrna_df$gene_name %in% i,]
    first_section = new_df[0:tail(which(new_df$gene_name == i), 1),]
    last_section = new_df[(tail(which(new_df$gene_name == i), 1)+1):dim(new_df)[1],]
    new_df = rbind(first_section, insert, last_section)
  }

  #### 6. Rename desired genes ####
  #################################
  # Rename desired genes (example from mouse genome): "Cers1"==>"Cers1_Gdf1" // "Chtf8" ==> "Chtf8_Derpc" // "Insl3" ==> "Insl3_Jak3" // "Pcdhga1" ==> "Pcdhg_all" // "Pcdha1" ==> "Pcdha_all" // "Ugt1a10" ==> "Ugt1a_all" // "4933427D14Rik" ==> "4933427D14Rik_Gm43951" // "Mkks" ==> "Mkks_plus"
  if(missing(gene_replacement)){

  }
  else{
    if(gene_replacement == "test_gene_replacement.csv"){
      gene_replacement <- system.file("extdata", "test_gene_replacement.csv", package = "ReferenceEnhancer")
          }

    gene_replacement <- read.csv(gene_replacement, header=T)

    old_names <- gene_replacement[,'old_name']
    new_names <- gene_replacement[,'new_name']


    for (i in 1:length(old_names)){
      new_df$transcript_name = stringr::str_replace_all(new_df$transcript_name, old_names[i], new_names[i])
      new_df$gene_name = stringr::str_replace_all(new_df$gene_name, old_names[i], new_names[i])
    }
  }

  #### 7. Export object to gtf file ####
  ######################################
  new_gtf = GenomicRanges::makeGRangesFromDataFrame(new_df, keep.extra.columns=TRUE)

  write_gtf(new_gtf, "optimized_reference.gtf")
  print("Optimized annotation reference has been saved in working directory as optimized_reference.gtf")

}
