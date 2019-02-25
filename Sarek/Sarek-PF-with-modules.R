#!/usr/bin/env Rscript
# call like:
# Rscript ~/dev/pathfindr/Sarek/Sarek-PF.R -r ~/dev/pathfindr/Sarek/reference_data/ -o bugger -s /data0/btb/P2233/test/
suppressWarnings(suppressMessages(library("optparse")))
write("\n *** Pathfindr: filtering and prioritizing somatic mutations \n",stdout())
library(tryCatchLog)
library(futile.logger)

# Command-line arguments to read
option_list = list(
  make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="reference directory name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=getwd(), 
              help="Sample directory [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$reference)) {
  print_help(opt_parser)
  stop("Reference dir must be supplied.\n", call.=FALSE)
}

library(pathfindr)
# we are using tictoc to benchmark runtimes all over the code
library(tictoc)
tic("Pathfindr")
# Print out sample name and other handy stuff
printList(list('Sample:',basename(getwd()),'Directory:',getwd(),'Date:',date()))


# reference is a command-line parameter
ref_data <- opt$reference
write(paste("Reference directory: ",ref_data),stdout())

sample_dir <- opt$sample
write(paste("Sample dir: ",sample_dir),stdout())

# we do want to print out tables
write_tables=TRUE

# by default and right now we are using only GRCh38
reference_genome='GRCh38'
write(paste("Reference genome: ",reference_genome),stdout())

# read all filenames we have in the current directory
files <- dir(path=sample_dir,recursive = T,full.names = T)
# get rid filenames that we are not interested in 
# TODO actually would be better to list only file that we are interested
remove=unique(c(grep(pattern = '.png',x = files),
	  grep(pattern = 'Annotation.old',x = files),
	  grep(pattern = '/work/',x = files)))
if (length(remove)>0) files <- files[-remove]

##################################### Sample data: we are storing results here #####################################
library(data.table)
sampleData <- data.table(
  date=date(),
  directory=sample_dir,
  name=basename(sample_dir )#,
  # patient='',
  # sample='',
  # freec_cnv_file,
  # freec_Tratio_file,
  # freec_Nratio_file,
  # freec_Tbaf_file,
  # ascat_Tratio_file
)

##################################### make chromosomes #####################################
# TODO: read it from the FASTA index file
if (reference_genome=='GRCh38') chrsz <- data.table(
  chr = paste0('chr',c("1", "2", "3", "4", "5", "6", "7", "8", 
          "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
          "20", "21", "22", "X", "Y")), 
  label = c("1", "2", "3", "4", "5", 
          "6", "7", "8", "9", "10", "11", "12", "13", 
          "14", "15", "16", "17", "18", "19", "20", 
          "21", "22", "X", "Y"),
  starts = c(0, 253902921, 500669562, 
             703689723, 898025604, 1084280925, 1259870526, 1424144247, 1574191848, 
             1717288929, 1855896810, 1995944331, 2134165212, 2253485013, 2363996694, 
             2470839615, 2565902256, 2654091777, 2739309138, 2802869979, 2871941700, 
             2923598421, 2979326382, 3140008743), 
  length = c(248902921, 241766641, 
             198020161, 189335881, 181255321, 170589601, 159273721, 145047601, 
             138097081, 133607881, 135047521, 133220881, 114319801, 105511681, 
             101842921, 90062641, 83189521, 80217361, 58560841, 64071721, 
             46656721, 50727961, 155682361, 56827081)
) 

##################################### Tumor Genes - General Ones #####################################
library(stringr)
tumorgenes_data = paste(ref_data,'cancer_gene_census.csv',sep='')
write(paste("Tumor genes data ",tumorgenes_data),stderr());
tumorgenes       = fread(tumorgenes_data,key='Gene Symbol')

tumorgenes$chr=paste0('chr',str_replace(string=tumorgenes$`Genome Location`,pattern = ':.*',replacement = ''))
temp=str_replace(string=tumorgenes$`Genome Location`,pattern = '.*:',replacement = '')
tumorgenes$start=as.numeric(str_replace(string=temp,pattern = '.*-',replacement = ''))
tumorgenes$end=as.numeric(str_replace(string=temp,pattern = '-.*',replacement = ''))
tumorgenes$cumstart=NA; tumorgenes$cumend=NA

for (i in 1:nrow(chrsz)) {
    ix <- tumorgenes$chr==chrsz$chr[i]
    tumorgenes$cumstart[ix] <- tumorgenes$start[ix]+chrsz$starts[i]
    tumorgenes$cumend[ix] <- tumorgenes$end[ix]+chrsz$starts[i]
}

##################################### Tumor Genes - Locals: these are usually cancer-specific extensions  #####################################
local_tumorgenes_data = paste(ref_data,'2018_gene_list_tere_ref.csv',sep='')
write(paste("Local tumor genes data ",local_tumorgenes_data),stderr());
local_tumorgenes=fread(local_tumorgenes_data,key='Gene')[Gene!='',1:2]

##################################### Create Tiers #####################################
alltumorgenes=unique(c(tumorgenes$`Gene Symbol`,local_tumorgenes$Gene))
alltier1=union(tumorgenes[Tier==1,`Gene Symbol`],local_tumorgenes[`Tier 1 and 2 for pediatric cancers final`==1,Gene])
alltier2=union(tumorgenes[Tier==2,`Gene Symbol`],local_tumorgenes[`Tier 1 and 2 for pediatric cancers final`==2,Gene])



# uncomment to enable modules 
tic("ControlFREEC")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             ControlFREEC                                                      x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
# read control freec result files
freec_Tratio_file <- files[grep(pattern = '[TR].hg38.pileup.gz_ratio.txt',files,perl = T)]
freec_Nratio_file <- files[grep(pattern = '[TR].hg38.pileup.gz_normal_ratio.txt',files,perl = T)][1]
freec_Tbaf_file <- files[grep(pattern = '[TR].hg38.pileup.gz_BAF.txt',files,perl = T)]
freec_Nbaf_file <- files[grep(pattern = '[BN].hg38.pileup.gz_BAF.txt',files,perl = T)][1]
freec_info_file <- files[grep(pattern = 'hg38.pileup.gz_info.txt',files,perl = T)]
freec_cnv_file <- files[grep(pattern = "^.*[TR]\\.hg38\\.pileup\\.gz_CNVs$",files,perl = T)]

write("ControlFREEC files:",stdout())
printList( list(freec_Tratio_file,freec_Nratio_file, freec_Tbaf_file, freec_Nbaf_file, freec_info_file, freec_cnv_file) )
tryCatchLog(loadControlFREEC(freec_Tratio_file,freec_Nratio_file, freec_Tbaf_file, freec_Nbaf_file, freec_info_file, freec_cnv_file))
toc()

#tic("ASCAT")
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#write("x                                             ASCAT                                                             x",stderr());
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#ascat_Tratio_file <- files[grep(pattern = "^.*Ascat.*[TR]\\.LogR$",files)]
#ascat_Nratio_file <- files[grep(pattern = "^.*Ascat.*[BN]\\.LogR$",files)][1]
#ascat_Tbaf_file <- files[grep(pattern = "^.*Ascat.*[TR]\\.BAF$",files)]
#ascat_Nbaf_file <- files[grep(pattern = "^.*Ascat.*[BN]\\.BAF$",files)][1]
#ascat_segment_file <- files[grep(pattern = "^.*Ascat.*[TR]\\.cnvs\\.txt$",files)]
#write("ASCAT files:",stdout())
#printList( list(ascat_Tratio_file,ascat_Nratio_file,ascat_Tbaf_file,ascat_Nbaf_file,ascat_segment_file) )
##loadASCAT(ascat_Tratio_file,ascat_Nratio_file,ascat_Tbaf_file,ascat_Nbaf_file,ascat_segment_file)
#toc()
#
#tic("Manta")
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#write("x                                             Manta                                                             x",stderr());
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#manta_tumor_file <- grep(pattern = ".*Manta_.*vs.*somaticSV.vcf.snpEff.ann.vep.ann.vcf$",files,value = T)
#manta_normal_file <- grep(pattern = ".*Manta_.*vs.*diploidSV.vcf.snpEff.ann.vep.ann.vcf$",files,value = T)[1]
#cosmic_fusions = fread(paste0(ref_data,'cosmic_fusions_table.csv'),key = 'name')
#swegen_manta_all=fread(paste0(ref_data,'swegen_sv_counts.csv'),key='name')
#write("Files for Manta structural variants: ", stdout())
#printList( list(manta_tumor_file,manta_normal_file) )
## TODO: sort it out to be local or a general on-demand object
#allfusionpairs=NULL
#allfusion=tumorgenes[grep('fusion',`Role in Cancer`),`Gene Symbol`]
#
#for (i in 1:length(allfusion)) {
#  t=trimws(strsplit(tumorgenes[allfusion[i],`Translocation Partner`],', ')[[1]])
#  if (length(t)>0) 
#		for (j in 1:length(t))
# 		allfusionpairs=c(allfusionpairs,paste(sort(c(allfusion[i],t[j])),collapse = ' '))
#}
##loadMantaTumor(manta_tumor_file)
##loadMantaNormal(manta_normal_file)
#toc()
#
#tic("Database read")
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#write("x                                     Reading databases for SNV and small indel prioritisation                  x",stderr());
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
## These are also needed for Strelka
#write("Reading SweGen SNP counts ...",stdout())
#snptable         = fread(paste0(ref_data,'swegen_snp_counts.csv'),key='name')
#write("Reading COSMIC tables ...",stdout())
#cosmic_coding    = fread(paste0(ref_data,'cosmic_coding_table.csv'),key = 'name')
#cosmic_noncoding = fread(paste0(ref_data,'cosmic_noncoding_table.csv'),key = 'name')
#cosmic_fusions   = fread(paste0(ref_data,'cosmic_fusions_table.csv'),key = 'name')
#alltsg = tumorgenes[grep('TSG',`Role in Cancer`),`Gene Symbol`]
#write("Reading hotspots ...",stdout())
#hotspots_snv = unique( fread(paste0(ref_data,'hotspots_v2_snv.csv'))[,.(Hugo_Symbol,Amino_Acid_Position)])[-grep('splice',Amino_Acid_Position)]
#hotspots_snv$pos <- as.numeric(hotspots_snv$Amino_Acid_Position)
#hotspots_inframe = unique(fread(paste0(ref_data,'hotspots_v2_inframe.csv'))[,.(Hugo_Symbol,Amino_Acid_Position)])
#hotspots_inframe$start=as.numeric(str_replace(string = hotspots_inframe$Amino_Acid_Position,pattern = '-[0-9]*',replacement = ''))
#hotspots_inframe$end=as.numeric(str_replace(string = hotspots_inframe$Amino_Acid_Position,pattern = '[0-9]*-',replacement = ''))
#near_hotspots=NULL
#for (i in 1:nrow(hotspots_inframe)) 
#  near_hotspots = c(near_hotspots,
#                    paste(hotspots_inframe$Hugo_Symbol[i],
#                    c(seq(hotspots_inframe$start[i]-2,hotspots_inframe$start[i]+2),
#                    hotspots_inframe$end[i]-2,hotspots_inframe$end[i]+2))
#                  )
#for (i in 1:nrow(hotspots_snv))
#  near_hotspots = c(near_hotspots,
#                    paste(hotspots_snv$Hugo_Symbol[i],
#                    seq(hotspots_snv$pos[i]-2,hotspots_snv$pos[i]+2))
#                  )
#
#toc()
#tic("Mutect2")
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#write("x                                             Mutect2 with GATK 3.8                                             x",stderr());
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#write("Reading MuTect2 file ...",stdout())
#mutect2_file <- grep(pattern = ".*mutect2_.*AF.*vep.ann.vcf$",files,value = T)
#
##loadMutect2(mutect2_file)
#toc()
#
#tic("Strelka")
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#write("x                                             Strelka                                                           x",stderr());
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
## snpEff + VEP Strelka snv file
#strelka_snv_file <- grep(pattern = ".*Strelka_.*_somatic_snvs.*AF.*vep.ann.vcf$",files,value = T)
## snpEff + VEP Strelka indel file
#strelka_indel_file <- grep(pattern = ".*Strelka_.*_somatic_indels.*AF.*vep.ann.vcf$",files,value = T)
#write(paste("Strelka SNV",strelka_snv_file),stderr());
#write(paste("Strelka indel",strelka_indel_file),stderr());
#
##loadStrelka()
#toc()
#
#tic("HaplotypeCaller")
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
#write("x                                             HaplotypeCaller                                                   x",stderr());
#write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
## Haplotype caller variant files
#haplotypecaller_T_file <- grep(pattern = ".*haplotypecaller.*[TR].*AF.*vep.ann.vcf$",files,value = T)
#haplotypecaller_N_file <- grep(pattern = ".*haplotypecaller.*[NB].*AF.*vep.ann.vcf$",files,value = T)
#write("Tumour and normal germline calls:",stdout())
#printList(list(haplotypecaller_T_file,haplotypecaller_N_file))
##loadHaplotypecaller()
toc()

toc()	# measuring the whole stuff
write("\nReady.",stdout())
