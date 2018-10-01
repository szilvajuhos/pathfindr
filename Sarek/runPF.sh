#!/bin/bash -x
R -e "rmarkdown::render('Sarek-view.Rmd',output_file='Sarek-output.html')"
#R -e "rmarkdown::render('Sarek-view.Rmd',knit_root_dir='~/reports/SAMPLES/P2233_104T_P2233_123N',output_file='Sarek-output.html')"


