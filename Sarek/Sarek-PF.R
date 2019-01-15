#!/usr/bin/env Rscript
# call like:
# Rscript ~/dev/pathfindr/Sarek/Sarek-PF.R -r ~/dev/pathfindr/Sarek/reference_data/ -o bugger -s /data0/btb/P2233/test/

library("optparse")

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


# we are using tictoc to benchmark runtimes all over the code
library(tictoc)

# Utility finction to print out lists
printList <- function(aList) {
	for (item in 1:length(aList)) { write(head(aList[[item]]),stdout()) }
}

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

# we are going to read quite a few files, and are expecting them to be there at 
# their place
checkFileExistence <- function(searchPattern, files) {
	fileToSearch <- grep(searchPattern, files, value=T)
	if(!length(fileToSearch)) {
		write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
		write(paste(" ERROR: File is missing. Search pattern is:",searchPattern),stderr());
		write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
		q(status=1)
	}
	fileToSearch
}

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

##################################### ControlFREEC #####################################
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

loadControlFREEC <- function(freec_Tratio_file,freec_Nratio_file, freec_Tbaf_file, freec_Nbaf_file, freec_info_file, freec_cnv_file) {
  cnvs=NULL
  freec_cnv=NULL
  samplename=NULL
  tratio=NULL; tbaf=NULL
  if (!is.na(freec_Tratio_file[1])) {
    # extract sample names
    for (s in 1:length(freec_Tratio_file)) samplename[s]=strsplit(basename(freec_Tratio_file[s]),'[.]')[[1]][1]
		write("Sample name, Tratio, Tbaf")
		printList(samplename)
    for (s in 1:length(freec_Tratio_file)) {
      # tumor log ratio
      temp <- fread(file = freec_Tratio_file[s])
      temp$sample=samplename[s]
      temp=temp[,c(10,1:6)]
      temp$CopyNumber[1:2] <- c(0,4) # to secure range of colors
      tratio=rbind(tratio,temp)
      setkey(tratio,'sample')
 
      # tumor BAF
      temp <- fread(freec_Tbaf_file[s])[,1:3]
      temp$sample=samplename[s]
      temp=temp[,c(4,1:3)]
      tbaf=rbind(tbaf,temp)
      setkey(tbaf,'sample')
    }
    # normal log ratio
    write("Normal log ratio and BAF",stdout());
    nratio <- fread(file = freec_Nratio_file)[,1:6]
    # normal BAF
    nbaf <- fread(freec_Nbaf_file)[,1:3]
 
    # manipulate for plot:
    write("Reshape data for plotting",stdout());
    tratio$cumstart <- tratio$Start
    nratio$cumstart <- nratio$Start
    tbaf$cumstart <- tbaf$Position
    nbaf$cumstart <- nbaf$Position
    tempchr=substr(chrsz$chr,4,6)
    for (i in 2:nrow(chrsz)) {
      ix <- tratio$Chromosome==tempchr[i]
      tratio$cumstart[ix] <- tratio$Start[ix]+chrsz$starts[i]
      ix <- nratio$Chromosome==tempchr[i]
      nratio$cumstart[ix] <- nratio$Start[ix]+chrsz$starts[i]
      ix <- tbaf$Chromosome==tempchr[i]
      tbaf$cumstart[ix] <- tbaf$Position[ix]+chrsz$starts[i]
      ix <- nbaf$Chromosome==tempchr[i]
      nbaf$cumstart[ix] <- nbaf$Position[ix]+chrsz$starts[i]
    }
 
    # modify for plot
    tratio$CopyNumber[tratio$CopyNumber>4] <- 4
    nratio$CopyNumber[nratio$CopyNumber>4] <- 4
    tratio$BAF[tratio$BAF<0.5 | tratio$BAF>1] <- NA
    nratio$BAF[nratio$BAF<0.5 | nratio$BAF>1] <- NA
    tratio$Chromosome=paste0('chr',tratio$Chromosome)
    nratio$Chromosome=paste0('chr',nratio$Chromosome)
    tbaf$Chromosome=paste0('chr',tbaf$Chromosome)
    nbaf$Chromosome=paste0('chr',nbaf$Chromosome)
 
    # smooth data for plot
    write("Smoothing",stdout());
    binned <- NULL
    for (sample in samplename) for (i in 1:nrow(chrsz)) {
      temp <- data.table(
        sample,
        chr=chrsz$chr[i],
        pos=seq(5e5,chrsz$length[i],5e5),
        cumpos=seq(5e5,chrsz$length[i],5e5)+chrsz$starts[i],
        tratio=NA,
        tmaf=NA)
      ctratio <- tratio[sample][Chromosome==chrsz$chr[i]]
      for (j in 1:nrow(temp)) {
        ix <- ctratio$Start>temp$pos[j]-5e5 & ctratio$Start<temp$pos[j]+5e5
        if (sum(ix)>20) {
          d <- density(ctratio$Ratio[ix],na.rm=T)
          temp$tratio[j] <- d$x[which.max(d$y)]
          t=ctratio$BAF[ix]
          t=t[!is.na(t)]
          if (length(t)>10) {
            d <- density(ctratio$BAF[ix],na.rm=T)
            temp$tmaf[j] <- d$x[which.max(d$y)]
          }
        }
      }
      binned <- rbind(binned,temp)
    }
    setkey(binned,'sample')
 
 
    # check which regions are worth reporting
    # cnv file
    write("Filtering CNVs",stdout());
    cnvs=NULL
    for (s in 1:length(samplename)) {
      temp <- as.data.table(read.csv(freec_cnv_file[s],sep='\t',header=F))[,1:8]
      names(temp) <- c('chr','start','end','copies','effect','genotype','value','type')
      temp$sample=samplename[s]
      temp=temp[,c(9,1:8)]
      cnvs=rbind(cnvs,temp)
      setkey(cnvs,'sample')
    }
    freec_loh=cnvs[grep('^A*$',genotype)][type=='somatic'][copies<4]
    freec_loh$chr=paste0('chr',freec_loh$chr)
    cnvs=cnvs[effect!='neutral'] # this removes LOH from CNV table
    cnvs$length=paste0(round((cnvs$end-cnvs$start)/1e6,2),'M')
    cnvs$rank_score <- 0
    cnvs$rank_terms <- ''
    cnvs$cancer_genes <- ''
 
    for (i in cnvs[end-start < 50e6, .I]) {
      genes=tumorgenes[chr==paste0('chr',cnvs$chr[i]) & start<cnvs$end[i] & end>cnvs$start[i],`Gene Symbol`] # genes in segment
      if (length(genes)>0) { # if any of our cancer genes
        cnvs$cancer_genes[i] <- paste(genes,collapse = ' ') # put in table
        if (any(genes %in% alltier1,na.rm = T)) {
          cnvs$rank_score[i]=2
          cnvs$rank_terms[i]='T1_gene'
        } else if (any(genes %in% alltier2,na.rm = T)) {
          cnvs$rank_score[i]=2
          cnvs$rank_terms[i]='T2_gene'
        }
      }
 
      # if the effect is focal
      if (cnvs[i,end-start]<3e6) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+1
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'focal')
      }
 
      # if the effect is loss/gain/amp etc
      if (cnvs$copies[i]==0) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+2
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'hz_loss')
      } else if (cnvs$copies[i] > 5) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+2
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'high_amp')
      } else if (cnvs$copies[i] > 2) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+1
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'dup')
      } else if (cnvs$effect[i] =='loss') {
        cnvs$rank_score[i]=cnvs$rank_score[i]+1
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'del')
      }
    }
		write("Ordering by rank. ",stdout())
    freec_cnv=cnvs[order(rank_score,decreasing = T)] # order by rank
    fwrite(freec_cnv,file=paste0(sampleData$name,'_freec_cnv.csv')) # write to file
 
		# TODO: make it callable if you need a HTML report
    # tableWrapper(freec_cnv) # table in report
		write("Ready. ",stdout())
	}
	# return with
	print(typeof(tratio))
	c(tratio,nratio)
}
# uncomment to enable ControlFREEC
#ratioList <- loadControlFREEC(freec_Tratio_file,freec_Nratio_file, freec_Tbaf_file, freec_Nbaf_file, freec_info_file, freec_cnv_file)
#CFTratio <- ratioList[1]
#CFNratio <- ratioList[2]


############################### Plot Control Freec tumor copy number log ratio, tBAF and nBAF ###################################
#{r plot_control_freec,fig.width=15,fig.height=7}

plotControlFREEC <- function(tratio,nratio) {
	print(typeof(tratio))
	library(ggplot2)
	if (exists('tratio')) for (sample in unique(tratio$sample)) {
		cat('Control Freec: ',sample)
   	par(mai=c(0,4,0,2))
   	temp <- tratio[sample]
   	temp$Ratio[temp$Ratio>3]=3
 
   	g=ggplot()+ylab('B allele ratio and Coverage Ratio')+xlab(sample)+
     	scale_color_gradientn(colours = c('violet','blue','black','red','orange'))+
     	expand_limits(x=c(0,3200e6))+
     	scale_y_continuous(limits=c(-1,3),
                        breaks = 0:2,
                        minor_breaks = seq(0,2.5,0.5),
                        labels = 0:2
     )
   	# normal logR below
   	g=g+geom_point(aes(x=cumstart,y=Ratio-0.3),col='darkgrey',data=nratio,alpha=1/5,shape='.')
   	# tumor logR
   	g=g+geom_point(aes(x=cumstart,y=Ratio,col=CopyNumber),data=temp,alpha=1/5,shape='.')
   	# add smoothed
   	#g=g+geom_line(aes(x=cumpos,y=tratio),stat="identity",colour="green",size=0.05,data=binned[sample])
 
   	g=g+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             axis.title.x=element_blank(),
             axis.text.x=element_blank(),
             axis.ticks.x=element_blank())
   	# add tBAF  <--- downsampled by 10
   	g=g+geom_point(data = tbaf[seq(1,nrow(tbaf),10)][sample],aes(x=cumstart,y=-0.5+0.5*BAF),
                  alpha=1/10,shape='.')
   	# add nBAF  <--- downsampled by 10
   	g=g+geom_point(data = nbaf[seq(1,nrow(nbaf),10)],aes(x=cumstart,y=-1+0.5*BAF),
                  alpha=1/10,shape='.')
   	# add chromosome lines
   	g=g+geom_segment(mapping = aes(x = starts,xend=starts,y= -1,yend=3),data=chrsz,inherit.aes = F,alpha=1/5)
   	# add chromosome labels
   	g=g+geom_text(aes(x=starts+0.5*length,y=0.1,label=label),data = chrsz,inherit.aes = F)
 	
		printList(g)
   	print(g)
	}
}

# TODO: fixit
#plotControlFREEC(CFTratio,CFNratio)


#### Scatter plot
#```{r scatter_plot_control_freec, echo=FALSE,fig.width=15,fig.height=8,warning=F}
#if (exists('binned')) for (sample in unique(binned$sample)) {
#   cat('Control Freec: ',sample)
#   binned$tmaf[binned$tmaf<0]=1
#   g=ggplot()+expand_limits(y=c(.5,1))+expand_limits(x=c(.3,3))+
#     geom_point(aes(x = tratio,y = tmaf,col=chr),data=binned[sample][tratio<1.8])+
#     xlab('Tumor (1Mb mean) DNA ratio')+ylab('Tumor (1Mb mean) major allele ratio')+
#     scale_y_continuous(breaks = seq(0.5,1,0.1)) +
#     scale_x_continuous(breaks = seq(0.5,2.5,0.1))+
#     theme(panel.grid.minor.x = element_blank(),
#           panel.grid.major.x = element_line(colour="grey", size=0.2),
#           panel.grid.major.y = element_line(colour="grey", size=0.2),
#           panel.grid.minor.y = element_blank()
#     )
#   print(g)
#}
#```

write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             ASCAT                                                             x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
ascat_Tratio_file <- files[grep(pattern = "^.*Ascat.*[TR]\\.LogR$",files)]
ascat_Nratio_file <- files[grep(pattern = "^.*Ascat.*[BN]\\.LogR$",files)][1]
ascat_Tbaf_file <- files[grep(pattern = "^.*Ascat.*[TR]\\.BAF$",files)]
ascat_Nbaf_file <- files[grep(pattern = "^.*Ascat.*[BN]\\.BAF$",files)][1]
ascat_segment_file <- files[grep(pattern = "^.*Ascat.*[TR]\\.cnvs\\.txt$",files)]
write("ASCAT files:",stdout())
printList( list(ascat_Tratio_file,ascat_Nratio_file,ascat_Tbaf_file,ascat_Nbaf_file,ascat_segment_file) )

loadASCAT <-function(ascat_Tratio_file,ascat_Nratio_file,ascat_Tbaf_file,ascat_Nbaf_file,ascat_segment_file) {
	ascat_cnv=NULL
	samplename=NULL
	if (length(ascat_Tratio_file)>0) {
		# extract sample names
		for (s in 1:length(ascat_Tratio_file)) samplename[s]=strsplit(basename(ascat_Tratio_file[s]),'[.]')[[1]][1]
		
		ascat_tratio=NULL; ascat_tbaf=NULL; ascat_cnv=NULL
		for (s in 1:length(samplename)) {
			# segmented copy number
			write("Reading segmented copy number",stdout())
			temp=fread(file = ascat_segment_file[s])
			temp$sample=samplename[s]
			temp=temp[,c(6,1:5)]
			temp$chr=paste0('chr',temp$chr)
			ascat_cnv=rbind(ascat_cnv,temp)
			setkey(ascat_cnv,'sample')
			these_cnvs=temp # for use below
			# tumor log ratio
			write("Reading tumor log ratio",stdout())
			temp=fread(file = ascat_Tratio_file[s])[,-1]
			#if (project=='AML') temp$Chr=paste0('chr',temp) 
			colnames(temp)[3]='LogR'
			temp$sample=samplename[s]
			temp=temp[order(Chr,Position),c(4,1:3)]
			temp$Ratio=2^temp$LogR
			temp$smoothed=runmed(x = temp$Ratio,k = 99)
			temp$CopyNumber=2
			temp$MinorCopy=1
			for (i in 1:nrow(these_cnvs)) {
				ix=temp$Chr==these_cnvs$chr[i] & temp$Position > these_cnvs$start[i] & temp$Position < these_cnvs$end[i]
				if (!any(ix)) 
					ix=temp$Chr==substr(these_cnvs$chr[i],4,6) & temp$Position > these_cnvs$start[i] & temp$Position < these_cnvs$end[i]
				temp$CopyNumber[ix]=these_cnvs$nMajor[i]+these_cnvs$nMinor[i]
				temp$MinorCopy[ix]=these_cnvs$nMinor[i]
			}
			temp$CopyNumber[temp$CopyNumber>4]=4
			temp$CopyNumber[1:2]=c(0,4)
			ascat_tratio=rbind(ascat_tratio,temp)
			setkey(ascat_tratio,'sample')
			# tumor BAF
			write("Reading tumor BAF",stdout())
			temp=fread(ascat_Tbaf_file[s])[,-1]
			colnames(temp)[3]='BAF'
			temp=temp[which(BAF!=0 & BAF!=1)]
			temp$sample=samplename[s]
			temp=temp[order(Chr,Position),c(4,1:3)]
			ascat_tbaf=rbind(ascat_tbaf,temp)
			setkey(ascat_tbaf,'sample')
		}
		# normal BAF
		write("Reading normal BAF",stdout())
		ascat_nbaf=fread(ascat_Nbaf_file)[,-1]
		colnames(ascat_nbaf)[3]='BAF'
		ascat_nbaf=ascat_nbaf[which(BAF!=0 & BAF!=1)]

		# join (T) log ratio and (T) BAF
		#ascat_tratio=merge(ascat_tratio,ascat_tbaf,all.x=T)
		#ascat_tratio=merge(ascat_tratio,ascat_nbaf,all.x=T)
		
		# manipulate for plot:
		ascat_tratio$cumstart=ascat_tratio$Position
		ascat_tbaf$cumstart=ascat_tbaf$Position
		ascat_nbaf$cumstart=ascat_nbaf$Position
		#ascat_cnv$LOH=''; ascat_cnv$LOH[ascat_cnv$nMinor==0]='LOH' <---- will use ascat_LOH instead
		for (i in 2:nrow(chrsz)) {
			ix=ascat_tratio$Chr==chrsz$chr[i]
			ascat_tratio$cumstart[ix]=ascat_tratio$Position[ix]+chrsz$starts[i]
			#no n...
			ix=ascat_tbaf$Chr==chrsz$chr[i]
			ascat_tbaf$cumstart[ix]=ascat_tbaf$Position[ix]+chrsz$starts[i]
			ix=ascat_nbaf$Chr==chrsz$chr[i]
			ascat_nbaf$cumstart[ix]=ascat_nbaf$Position[ix]+chrsz$starts[i]
		}
	 

		
		# smooth data for view
		ascat_binned=NULL
		write("Smooth data for view",stdout())
		for (sample in samplename) for (i in 1:nrow(chrsz)) {
			temp=data.table(
				sample,
				chr=chrsz$chr[i],
				pos=seq(5e5,chrsz$length[i],5e5),
				cumpos=seq(5e5,chrsz$length[i],5e5)+chrsz$starts[i],
				tratio=NA,
				tmaf=NA)
			
			tempLogR=ascat_tratio[sample][Chr==chrsz$chr[i]]
			tempBAF=ascat_tbaf[sample][Chr==chrsz$chr[i]]
			tempBAF$BAF=0.5+abs(tempBAF$BAF-0.5)
			for (j in 1:nrow(temp)) {
				ix = tempLogR$Position>temp$pos[j]-5e5 & tempLogR$Position<temp$pos[j]+5e5
				if (sum(ix)>10) {
					d <- density(2^tempLogR$LogR[ix],na.rm=T)
					temp$tratio[j] <- d$x[which.max(d$y)]
				}
				ix = tempBAF$Position>temp$pos[j]-5e5 & tempBAF$Position<temp$pos[j]+5e5
				if (sum(!is.na(tempBAF$BAF[ix]))>20) {
					d=density(tempBAF$BAF[ix],na.rm=T)
					temp$tmaf[j]=d$x[which.max(d$y)]
				}
			}
			ascat_binned=rbind(ascat_binned,temp)
		}
		setkey(ascat_binned,'sample')

		
		ascat_loh=ascat_cnv[nMinor==0 & nMajor<4]
		ascat_loh$chr=paste0('chr',ascat_loh$chr)
		if (!exists('freec_loh')) freec_loh=ascat_loh # to make sure LOH is available if theres no freec
		colnames(ascat_loh)[3:4]=c('start','end')
		cnvs=ascat_cnv[nMajor+nMinor != 2] # this removes cnLOH
		cnvs$length=paste0(round((cnvs$endpos-cnvs$startpos)/1e6,2),'M')
		cnvs$rank_score <- 0
		cnvs$rank_terms <- ''
		cnvs$cancer_genes <- ''
		
		for (i in which(cnvs$endpos-cnvs$startpos < 50e6)) {
			genes=tumorgenes[chr==cnvs$chr[i] & start<cnvs$end[i] & end>cnvs$start[i],`Gene Symbol`] # genes in segment
			if (length(genes)>0) { # if any of our cancer genes
				cnvs$cancer_genes[i] <- paste(genes,collapse = ' ') # put in table
				if (any(genes %in% alltier1,na.rm = T)) {
					cnvs$rank_score[i]=2
					cnvs$rank_terms[i]='T1_gene'
				} else if (any(genes %in% alltier2,na.rm = T)) {
					cnvs$rank_score[i]=2
					cnvs$rank_terms[i]='T2_gene'
				}
			}
			
			# if the effect is focal
			if (cnvs[i,endpos-startpos]<3e6) {
				cnvs$rank_score[i]=cnvs$rank_score[i]+1
				cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'focal')
			}
			
			
			
			# if the effect is loss/gain/amp etc
			if (cnvs[i,nMajor+nMinor]==0) {
				cnvs$rank_score[i]=cnvs$rank_score[i]+2
				cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'hz_loss')
			} else if (cnvs[i,nMajor+nMinor] > 5) {
				cnvs$rank_score[i]=cnvs$rank_score[i]+2
				cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'high_amp')
			} else if (cnvs[i,nMajor+nMinor] > 2) {
				cnvs$rank_score[i]=cnvs$rank_score[i]+1
				cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'dup')
			} else if (cnvs[i,nMajor+nMinor] < 2) {
				cnvs$rank_score[i]=cnvs$rank_score[i]+1
				cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'del')
			}
		}
		
		ascat_cnv=cnvs[order(rank_score,decreasing = T)]
		if (write_tables) fwrite(ascat_cnv,file=paste0(sampleData$name,'_ascat_cnv.csv'))

		#tableWrapper(ascat_cnv)
	}
}
#loadASCAT(ascat_Tratio_file,ascat_Nratio_file,ascat_Tbaf_file,ascat_Nbaf_file,ascat_segment_file)

write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             Manta                                                             x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
manta_tumor_file <- grep(pattern = ".*Manta_.*vs.*somaticSV.vcf.snpEff.ann.vep.ann.vcf$",files,value = T)
manta_normal_file <- grep(pattern = ".*Manta_.*vs.*diploidSV.vcf.snpEff.ann.vep.ann.vcf$",files,value = T)[1]
cosmic_fusions=fread(paste0(ref_data,'cosmic_fusions_table.csv'),key = 'name')

# TODO: sort it out to be local or a general on-demand object
allfusionpairs=NULL
allfusion=tumorgenes[grep('fusion',`Role in Cancer`),`Gene Symbol`]
for (i in 1:length(allfusion)) {
  t=trimws(strsplit(tumorgenes[allfusion[i],`Translocation Partner`],', ')[[1]])
  if (length(t)>0) for (j in 1:length(t))
  allfusionpairs=c(allfusionpairs,paste(sort(c(allfusion[i],t[j])),collapse = ' '))
}

loadManta <- function(manta_tumor_file,manta_normal_file) {
	library(VariantAnnotation)

	# get SweGen AFs: TODO make it optional
	swegen_manta_all=fread(paste0(ref_data,'swegen_sv_counts.csv'),key='name')
	# first collect PASS ids from all samples
	allpass=NULL
	if (length(manta_tumor_file)>0) for (s in 1:length(manta_tumor_file)) {
		vcf=readVcf(file = manta_tumor_file[s],genome = reference_genome)
		pass=rowRanges(vcf)$FILTER=='PASS'
		allpass=c(allpass,names(vcf)[pass])
	}
	# then collect variants...
	manta_tumor_table=NULL
	if (length(manta_tumor_file)>0) for (s in 1:length(manta_tumor_file)) {
		sample=strsplit(basename(manta_tumor_file[s]),'[.]')[[1]][1]
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
		
		
			# Filter out variants seen 2+ (?) times in reference data
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
	
		# Prepare ranking and (some more) filtering
		if (!is.null(manta_tumor_table)) if (nrow(manta_tumor_table)>0) {
			# ## Make table with most important info for ranking, and report
			# selected <- unique(manta_tumor_table[,.(ID,sample,SVTYPE,chr,start,REF,ALT,AFreq,PairedReadSupport_T,SplitReadSupport_T,
			#                        Swegen_count,Rank_score='',Rank_terms='',Gene_Name,Annotation,Annotation_Impact)])
			selection=manta_tumor_table[!is.na(Annotation)]
			
			if (nrow(selection)>0) {
				# Known cancer genes affect ranking by +1
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
				
			# Add high impact
			ix=selection$Annotation_Impact=='HIGH'
			if (any(ix)) {
				selection$rank_score[ix]=selection$rank_score[ix]+2
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high_impact')
			}
			# Add moderate impact
			ix=selection$Annotation_Impact=='MODERATE'
			if (any(ix)) {
				selection$rank_score[ix]=selection$rank_score[ix]+1
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'moderate_impact')
			}
			
			
			# +1 for focal
			ix=selection$chr==selection$altchr & selection$end-selection$start < 3e6 & selection$altpos-selection$start < 3e6
			if (any(ix)) {
				selection$rank_score[ix]=selection$rank_score[ix]+1
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'focal')
			}
			
			# cosmic_>xx, fusion_gene or just fusion
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
			
			# -1 for ablation if long del
			ix=intersect(
				which(selection$SVTYPE=='DEL' & selection$chr==selection$altchr & abs(selection$altpos-selection$start) > 3e6),
				grep('ablation',selection$Annotation)
			)
			if (length(ix)>0) {
				selection$rank_score[ix]=selection$rank_score[ix]-1
				selection$rank_terms[ix]=paste(selection$rank_terms[ix],'long_del')
			}
			
			# -1 for duplication if long dup
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
			if (write_tables) 
				fwrite(manta_tumor_selected[,-c('ANN','Gene_ID')][rank_score>1],file=paste0(sampleData$name,'_manta_tumor.csv'))
			#tableWrapper(manta_tumor_selected[,-c('ANN','Gene_ID')][rank_score>1])
		} 
	}
}

loadManta(manta_tumor_file,manta_normal_file)

