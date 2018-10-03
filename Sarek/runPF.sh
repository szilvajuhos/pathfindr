#!/bin/bash -x
#R -e "rmarkdown::render('Sarek-view.Rmd',output_file='Sarek-output.html')"
#R -e "rmarkdown::render('Sarek-view.Rmd',knit_root_dir='SAMPLES/PFtest',output_file='Sarek-output.html') " 

usage() { echo "Usage: $0 [-r reference_data ] [-s sample_dir] [-o output.html] " 1>&2; exit 1; }

while getopts "s:o:r:" p; do
	case "${p}" in
		r)
			reference_dir=${OPTARG}
			;;
		s)
			sample_dir=${OPTARG}
			;;
		o)
			outputFile=${OPTARG}
			;;
		*)
			usage
			;;
	esac
done
shift $((OPTIND-1))

if [ -z ${sample_dir} ] || [ ! -d ${sample_dir} ]; then echo "Incorrect sample directory "; usage; fi
if [ -z ${reference_dir} ] || [ ! -d ${reference_dir} ]; then echo "Incorrect reference directory "; usage; fi
if [ -z ${outputFile} ]; then echo "Missing output file name"; usage; fi


R -e "rmarkdown::render('Sarek-view.Rmd', knit_root_dir='${sample_dir}', output_file='${outputFile}', params=list(reference_data='${reference_dir}'))" 


