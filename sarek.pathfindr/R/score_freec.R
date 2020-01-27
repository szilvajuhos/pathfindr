score_freec <- function() {
  cat(" ------------ Score Control-FREEC function -------------\n")
  cat("Reference is: ",reference_genome,'\n')
  freec_result_files <- freec_files()
  
  cat('Selected files for Control-FREEC in',getwd(),':\n\n')
  cat("Control-FREEC Tumor Ratio file:  ", freec_result_files["freec_Tratio_file"], '\n')
  cat("Control-FREEC Normal Ratio file: ", freec_result_files["freec_Nratio_file"], '\n')
  cat("Control-FREEC Tumor BAF file:    ", freec_result_files["freec_Tbaf_file"], '\n')
  cat("Control-FREEC Normal BAF file:   ", freec_result_files["freec_Nbaf_file"], '\n')
  cat("Control-FREEC Info file:      ", freec_result_files["freec_info_file"], '\n')
  cat("Control-FREEC CNV file:      ", freec_result_files["freec_cnv_file"], '\n')
  load_freec(freec_result_files)
}
