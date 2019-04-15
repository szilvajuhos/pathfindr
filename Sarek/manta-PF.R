#!/usr/bin/env Rscript
# call like:
# Rscript ~/dev/pathfindr/Sarek/manta-PF.R -r ~/dev/pathfindr/Sarek/reference_data/ -o bugger -s /data0/btb/P2233/test/
options(warn=0)
suppressWarnings(suppressMessages(library("optparse")))

write("\n *** Pathfindr: filtering and prioritizing somatic mutations \n",stdout())

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
opt$reference <- '~/dev/pathfindr/Sarek/reference_data/'
opt$sample <- '../TEST'
opt$out <- 'bugger'

if (is.null(opt$reference)) {
  print_help(opt_parser)
  stop("Reference dir must be supplied.\n", call.=FALSE)
}



# we are using tictoc to benchmark runtimes all over the code
library(tictoc)
tic("Pathfindr")
# Utility finction to print out lists
printList <- function(aList) {
	for (item in 1:length(aList)) { write(head(aList[[item]]),stdout()) }
}

# Print out sample name and other handy stuff
printList(list('Sample:',basename(getwd()),'Directory:',getwd(),'Date:',date()))


# reference is a command-line parameter
ref_data <- opt$reference
#ref_data <- "reference_data/"
write(paste("Reference directory: ",ref_data),stdout())

sample_dir <- opt$sample
#sample_dir <- "../TEST"
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
    grep(pattern = 'Annotation/Manta/VEP.CADD',x = files),
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
##################################### ControlFREEC #####################################

loadMantaTumor <- function(manta_tumor_file) {
  suppressWarnings(suppressMessages(library(VariantAnnotation)))

	# get SweGen AFs: TODO make it optional
	swegen_manta_all=fread(paste0(ref_data,'swegen_sv_counts.csv'),key='name')
	write("Number of entries in SweGen SV swegen_manta_all",stderr())
	write(str(dim(swegen_manta_all)),stderr())
	# first collect PASS ids from all samples
	allpass=NULL
	if (length(manta_tumor_file)>0) for (s in 1:length(manta_tumor_file)) {
		vcf=readVcf(file = manta_tumor_file[s],genome = reference_genome)	# @TODO save the VCF object, we want to read only once
		write(paste("Number of VCF entries:",dim(vcf)),stdout())
		# for entries with PASS flag it as TRUE 
		pass=rowRanges(vcf)$FILTER=='PASS'
		write(paste("Number of PASS-ed calls:",sum(pass)),stdout())
		allpass=c(allpass,names(vcf)[pass])
	}
	write(paste("Number of ALL PASS-ed Manta calls to process:",summary(allpass)),stdout())

	# then collect variants...
	manta_tumor_table=NULL
	if (length(manta_tumor_file)>0) for (s in 1:length(manta_tumor_file)) {
		sample=strsplit(basename(manta_tumor_file[s]),'[.]')[[1]][1]
		cat('Parsing tumor file\n')
		cat(manta_tumor_file[s],'\n')
	
		vcf=readVcf(file = manta_tumor_file[s],genome = reference_genome)
		vcf=vcf[names(vcf) %in% allpass]
	
		if (length(vcf)>0) {
			g=geno(vcf)
			inf=info(vcf)
			rr=rowRanges(vcf)
		
			# Read counts
			srtable=as.data.table(g$SR)
			colnames(srtable)=paste0('SplitReadSupport_',c('N','T')) # Verified on a test sample
			prtable=as.data.table(g$PR)
			colnames(prtable)=paste0('PairedReadSupport_',c('N','T'))
		
			# Calculate allele ratio
			temp=data.table(
				pr.ref=sapply(prtable$PairedReadSupport_T,"[",1),
				sr.ref=sapply(srtable$SplitReadSupport_T,"[",1),
				pr.alt=sapply(prtable$PairedReadSupport_T,"[",2),
				sr.alt=sapply(srtable$SplitReadSupport_T,"[",2)
			)
			AFreq=apply(temp[,3:4],1,sum,na.rm=T)/
				(apply(temp[,3:4],1,sum,na.rm=T)+apply(temp[,1:2],1,sum,na.rm=T))
		
		
			# make table of variants
			sv_table=cbind(
				data.table(ID=as.character(names(vcf)),
				sample,
				chr=as.character(seqnames(rr))),
				as.data.table(rr)[,-1],
				AFreq,prtable,srtable,
				data.table(Swegen_count=rep(NA,length(vcf)),
					rank_score=0,
					rank_terms='',
					LOH=''
				),
				as.data.table(inf)
			)
			sv_table$cumstart=sv_table$start
			sv_table$cumend=sv_table$end
			sv_table$altchr=sv_table$chr
			sv_table$altpos=sv_table$END
			sv_table$plot=F
			sv_table$arc= -1
	
			# now we have a table of SVs with all the data
			# Filter out variants seen 2+ (?) times in reference data
			cat("Now get rid of duplicates ...\n")
			## Key has only chr,start,end
			key=sv_table[,c('chr','start','end')]
			key$chr=substr(key$chr,4,6)
			key$imprecise='(pr)'
			## If imprecise, round the pos to 10
			ix=sv_table$IMPRECISE==T
			key$imprecise[ix]='(impr)'
			key$start[ix]=round(key$start[ix]/10)*10
			key$end[ix]=round(key$end[ix]/10)*10
			key=paste(key$chr,key$start,key$end,key$imprecise)
		
			# put in Swegen count
			sv_table$Swegen_count=swegen_manta_all[key,value]
			sv_table$Swegen_count[is.na(sv_table$Swegen_count)]=0
			## do the filter
			sv_table <- sv_table[Swegen_count<2]
		}

		cat("Collected SVs so far:\n ")	
		cat(paste(sum(sv_table$SVTYPE=='INV'),"Inversions\n"))
		cat(paste(sum(sv_table$SVTYPE=='DEL'),"Deletions\n"))
		cat(paste(sum(sv_table$SVTYPE=='INS'),"Inserts\n"))
		cat(paste(sum(sv_table$SVTYPE=='BND'),"Translocations\n"))
		cat(paste(sum(sv_table$SVTYPE=='DUP'),"Duplications\n"))
		
		if (exists('sv_table')) if (nrow(sv_table)>0) {
			# loop through all and extract endpoint chr and pos   <------ TODO: This one must be remade.. (refactor)
				for (i in 1:nrow(sv_table)) try( {                    # <----  sometimes error here, so try..
					t=strsplit(x = sv_table$ALT[[i]],split = ':')[[1]]
					if (length(t)>1 & t[1]!="<DUP") {
						tchr=str_extract(t[1],'[0-9,X,Y]*$')
						sv_table$altchr[i] <- paste0('chr',tchr) #else sv_table$altchr[i] <- tchr for non-BTB thingy
						tt=str_extract(t[2],'^[0-9]*')
						sv_table$altpos[i]=as.numeric(tt)
					}
				}, silent=T)
			# for each chromosome get cumulative pos
			sv_table$altcumpos=sv_table$altpos
			for (i in 1:nrow(chrsz)) {
				ix=sv_table$chr==chrsz$chr[i]
				if (sum(ix)>0) {
					sv_table$cumstart[ix]=sv_table$start[ix]+chrsz$starts[i]
					sv_table$cumend[ix]=sv_table$end[ix]+chrsz$starts[i]
				}
				ix=sv_table$altchr==chrsz$chr[i]
				if (sum(ix)>0) {
					sv_table$altcumpos[ix]=sv_table$altpos[ix]+chrsz$starts[i]
				}
			}
			# decide how it is to be plotted (not represented elsewhere, up or down arc)
			for (i in 1:nrow(sv_table)) {
				if (sv_table$chr[i]==sv_table$altchr[i]) {
					# intrachromosomal: plot always, with "positive" arc
					sv_table$plot[i]=T
					sv_table$arc[i]=1 
				} else if (sv_table$altcumpos[i] > sv_table$cumstart[i]) {
					# interchromosomal: plot if mate is to right (else mate will be plotted) with negative arc
					sv_table$plot[i]=T 
					sv_table$arc[i]=-1 
				}
			}
		
			# snpEff annotation is in the ANN column.
			h=strsplit(info(header(vcf))['ANN',][,3],'annotations: \'')[[1]][2]
			snpEff_header=trimws(strsplit(h,'\\|')[[1]])
	
			# snpEff annotation is put in snpEff_table
			snpEff_table=matrix(data = NA,nrow = length(unlist(sv_table$ANN)),ncol = length(snpEff_header)+1)
			colnames(snpEff_table)=c('ID',snpEff_header)
			row=1
			for (i in 1:nrow(sv_table)) { # for each variant
				# for each VEP annotation:
				for (j in 1:length(sv_table$ANN[[i]]))  if (length(sv_table$ANN[[i]])>0) {
					line=strsplit(sv_table$ANN[[i]][j],'\\|')[[1]]
					snpEff_table[row,1]=sv_table$ID[i]
					snpEff_table[row,1+(1:length(line))]=line
					row=row+1
				}
			}
			snpEff_table=unique(as.data.table(snpEff_table))
			
			# Filter out annotations of certain types where the ID has >N annotations of that type
			ids=unique(snpEff_table$ID)
			for (id in ids) {
				common=c('protein_protein_contact', 'duplication', 'structural_interaction_variant', 'inversion',
						'transcript_ablation', 'feature_ablation', 'sequence_feature', 
						'intergenic_region', 'downstream_gene_variant', 'upstream_gene_variant')
				annotations=snpEff_table[ID==id,Annotation] # the annotations of this variant
				table=table(annotations[annotations %in% common])
				table=table[table>20] # the common annotations that appear >N times for this variant
				if (length(table)>0) {
					remove=which(snpEff_table$ID==id & snpEff_table$Annotation %in% names(table))
					snpEff_table=snpEff_table[-remove[-1],] # saves one to make sure the variant has some annotation
				}
			}
				# Add to data (all samples, one table)
				manta_tumor_table=rbind(manta_tumor_table,merge(sv_table,snpEff_table,by='ID',all=F))
				setkey(manta_tumor_table,'sample')
			}
		}# done parsing each sample
	cat(paste("Parsing done, collected ",length(manta_tumor_table),"entries for further processing\n"))
	
		# Prepare ranking and (some more) filtering
		if (!is.null(manta_tumor_table)) if (nrow(manta_tumor_table)>0) {
			# ## Make table with most important info for ranking, and report
			# selected <- unique(manta_tumor_table[,.(ID,sample,SVTYPE,chr,start,REF,ALT,AFreq,PairedReadSupport_T,SplitReadSupport_T,
			#                        Swegen_count,Rank_score='',Rank_terms='',Gene_Name,Annotation,Annotation_Impact)])
			selection=manta_tumor_table[!is.na(Annotation)]
			
			if (nrow(selection)>0) {
				# Known cancer genes affect ranking by +1
				cat("Increase rank by being a known cancer gene\n")
				for (gene in unique(selection$Gene_Name)) if (!is.na(gene)) if (gene!='') {
					ix=selection$Gene_Name==gene
					gene=sort(strsplit(gene,'&')[[1]])
					# cosmic/local Tier2:
					if (any(gene %in% alltier2)) {
						selection$rank_score[ix]=2
						selection$rank_terms[ix]='T2_gene'
					}
					# tier 1 top priority:
					if (any(gene %in% alltier1)) {
						selection$rank_score[ix]=2
						selection$rank_terms[ix]='T1_gene'
					}
				}
				
			cat("Add high impact\n")
			ix=selection$Annotation_Impact=='HIGH'
			if (any(ix)) {
				selection$rank_score[ix]=selection$rank_score[ix]+2
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high_impact')
			}
			cat("Add moderate impact\n")
			ix=selection$Annotation_Impact=='MODERATE'
			if (any(ix)) {
				selection$rank_score[ix]=selection$rank_score[ix]+1
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'moderate_impact')
			}
			
			
			cat(" +1 for focal\n")
			ix=selection$chr==selection$altchr & selection$end-selection$start < 3e6 & selection$altpos-selection$start < 3e6
			if (any(ix)) {
				selection$rank_score[ix]=selection$rank_score[ix]+1
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'focal')
			}
			
			cat("cosmic_>xx, fusion_gene or just fusion\n")
			ix=grep('fusion',selection$Annotation)
			if (length(ix)>0) for (i in ix) if (selection$Gene_Name[i]!='') {
				gene=sort(strsplit(selection$Gene_Name[i],'&')[[1]])
				if (paste(gene,collapse = ' ') %in% cosmic_fusions[value>5,name]) {
					selection$rank_score[i]=selection$rank_score[i]+3
					selection$rank_terms[i]=paste(selection$rank_terms[i],'cosmic_>5')
				} else if (paste(gene,collapse = ' ') %in% cosmic_fusions$name) {
					selection$rank_score[i]=selection$rank_score[i]+2
					selection$rank_terms[i]=paste(selection$rank_terms[i],'cosmic_>1')
				}
				if (paste(gene,collapse = ' ') %in% allfusionpairs) {
					selection$rank_score[i]=selection$rank_score[i]+2
					selection$rank_terms[i]=paste(selection$rank_terms[i],'CGC_fusion')
				} else if (any(gene %in% allfusion)) {
					selection$rank_score[i]=selection$rank_score[i]+1
					selection$rank_terms[i]=paste(selection$rank_terms[i],'partial_CGC_fusion')
				} else {
					selection$rank_score[i]=selection$rank_score[i]+0
					selection$rank_terms[i]=paste(selection$rank_terms[i],'fusion')
				}
			}
			
			cat("-1 for ablation if long del\n")
			ix=intersect(
				which(selection$SVTYPE=='DEL' & selection$chr==selection$altchr & abs(selection$altpos-selection$start) > 3e6),
				grep('ablation',selection$Annotation)
			)
			if (length(ix)>0) {
				selection$rank_score[ix]=selection$rank_score[ix]-1
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'long_del')
			}
			
			cat("-1 for duplication if long dup\n")
			ix=intersect(
				which(selection$SVTYPE=='DUP' & selection$chr==selection$altchr & abs(selection$altpos-selection$start) > 10e6),
				grep('duplication',selection$Annotation)
			)
			if (length(ix)>0) {
				selection$rank_score[ix]=selection$rank_score[ix]-1
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'long_dup')
			}
	 
			firstcols=c('ID','sample','Gene_Name','rank_score','rank_terms','LOH','AFreq','Annotation','Annotation_Impact','Swegen_count')
			cols=colnames(selection)
			setcolorder(x = selection,neworder = c(firstcols,cols[!cols %in% firstcols]))
			
			manta_tumor_selected <- selection[order(Feature_ID)][order(rank_score,decreasing = T)]
			cat("Selected Manta SVs:")
			cat(length(manta_tumor_selected))
			#write(str(head(manta_tumor_selected)),stderr())
      outfile <- paste0(sampleData$name,'_manta_tumor.csv')
			# exclude ANN and Gene_ID
			manta_tumor_selected_ranked <- subset(manta_tumor_selected[,-c('ANN','Gene_ID')], rank_score > 1)
			#manta_tumor_selected_ranked <- subset(manta_tumor_selected, rank_score > 1)
			#View(head(manta_tumor_selected_ranked))
			#a <- readLines("stdin",n=1);
		browser()	
      write(paste0(" Writing to ",outfile),stdout())
			#write.table(manta_tumor_selected_ranked,file=outfile,quote=TRUE,sep=",",na="NA")
    	#write.csv(manta_tumor_selected_ranked,file=outfile,na="")
      fwrite(manta_tumor_selected[,-c('ANN','Gene_ID')][rank_score>1],file=outfile)
      write(paste0(" *** Manta results written to ",outfile),stdout())
			#tableWrapper(manta_tumor_selected[,-c('ANN','Gene_ID')][rank_score>1])
		} 
	}
}

tic("Manta")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             Manta                                                             x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
manta_tumor_file <- grep(pattern = ".*VEP\\.CADD.*Manta_.*vs.*somaticSV.*ann.vcf$",files,value = T)
manta_normal_file <- grep(pattern = ".*VEP\\.CADD.*Manta_.*vs.*diploidSV.*ann.vcf$",files,value = T)[1]
cosmic_fusions = fread(paste0(ref_data,'cosmic_fusions_table.csv'),key = 'name')
swegen_manta_all=fread(paste0(ref_data,'swegen_sv_counts.csv'),key='name')
write("Files for Manta structural variants: ", stdout())
write(paste("manta_tumor_file: ",manta_tumor_file), stdout())
write(paste("manta_normal_file: ",manta_normal_file), stdout())
#printList( list(manta_tumor_file,manta_normal_file) )
# TODO: sort it out to be local or a general on-demand object
allfusionpairs=NULL
allfusion=tumorgenes[grep('fusion',`Role in Cancer`),`Gene Symbol`]

for (i in 1:length(allfusion)) {
  t=trimws(strsplit(tumorgenes[allfusion[i],`Translocation Partner`],', ')[[1]])
  if (length(t)>0) 
		for (j in 1:length(t))
 		allfusionpairs=c(allfusionpairs,paste(sort(c(allfusion[i],t[j])),collapse = ' '))
}
toc()

tic("Database read")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                     Reading databases for SNV and small indel prioritisation                  x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
# These are also needed for Strelka
write("Reading COSMIC tables ...",stdout())
cosmic_fusions   = fread(paste0(ref_data,'cosmic_fusions_table.csv'),key = 'name')

toc()


write("Starting Manta tumor processing",stdout())
loadMantaTumor(manta_tumor_file)
#loadMantaNormal(manta_normal_file)

write("\nReady.",stdout())
