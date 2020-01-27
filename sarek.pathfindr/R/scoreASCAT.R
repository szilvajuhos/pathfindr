ascat_files <- function() {
  ascat_result_files <- dir(recursive = T,full.names = T)
  file_names_to_remove = unique(c(grep(pattern = '.png', x = ascat_result_files),
                  grep(pattern = 'Annotation.old',x = ascat_result_files),
                  grep(pattern = '/work/',x = ascat_result_files),
                  grep(pattern = 'VEP.CADD/Manta',x = ascat_result_files),
                  grep(pattern = 'VEP/haplo',x = ascat_result_files),
                  grep(pattern = 'VEP/mute',x = ascat_result_files),
                  grep(pattern = 'VEP/Strelk',x = ascat_result_files),
                  grep(pattern = 'old.ControlFreec',x = ascat_result_files),
                  grep(pattern = 'old.BTB_ControlFreec',x = ascat_result_files)
  ))
  
  if (length(file_names_to_remove)>0) ascat_result_files <- ascat_result_files[-file_names_to_remove]
  
  # Ascat files
  ascat_Tratio_file <- ascat_result_files[grep(pattern = "^.*Ascat.*[TR2]\\.LogR$",ascat_result_files)]
  ascat_Nratio_file <- ascat_result_files[grep(pattern = "^.*Ascat.*[BN1]\\.LogR$",ascat_result_files)][1]
  ascat_Tbaf_file <- ascat_result_files[grep(pattern = "^.*Ascat.*[TR2]\\.BAF$",ascat_result_files)]
  ascat_Nbaf_file <- ascat_result_files[grep(pattern = "^.*Ascat.*[BN1]\\.BAF$",ascat_result_files)][1]
  ascat_segment_file <- ascat_result_files[grep(pattern = "^.*Ascat.*[TR2]\\.cnvs\\.txt$",ascat_result_files)]  
  c( "ascat_Tratio_file" = ascat_Tratio_file,
     "ascat_Nratio_file" = ascat_Nratio_file,
     "ascat_Tbaf_file" = ascat_Tbaf_file,
     "ascat_Nbaf_file" = ascat_Nbaf_file,
     "ascat_segment_file" = ascat_segment_file
     )
}

score_ascat <- function() {
  cat(" ------------ Score ASCAT function -------------\n")
  cat("Reference is: ",reference_genome,'\n')
  ascat_result_files <- ascat_files()
  cat('Selected files for ASCAT in',getwd(),':\n\n')
  cat("ASCAT Tumor Ratio file:  ", ascat_result_files["ascat_Tratio_file"], '\n')
  cat("ASCAT Normal Ratio file: ", ascat_result_files["ascat_Nratio_file"], '\n')
  cat("ASCAT Tumor BAF file:    ", ascat_result_files["ascat_Tbaf_file"], '\n')
  cat("ASCAT Normal BAF file:   ", ascat_result_files["ascat_Nbaf_file"], '\n')
  cat("ASCAT Segment file:      ", ascat_result_files["ascat_segment_file"], '\n')
  load_ascat(ascat_result_files)
}

