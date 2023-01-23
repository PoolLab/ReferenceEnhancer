write_bed <- function(output_data, file_name){
  path = "."
  out_path <- file.path(path, file_name)
  rtracklayer::export.bed(output_data, con = out_path)
}
