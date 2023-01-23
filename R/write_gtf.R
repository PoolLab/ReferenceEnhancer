write_gtf <- function(output_data, file_name){
  #path = "."
  #out_path <- file.path(path, file_name)
  rtracklayer::export(output_data, file_name, format = "gtf")
}
