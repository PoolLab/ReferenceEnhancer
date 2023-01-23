write_csv <- function(output_data, file_name){
  path = "."
  out_path <- file.path(path, file_name)
  write.csv(output_data, out_path)
}
