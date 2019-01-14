#!/usr/bin/env Rscript
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

##################################### formatting utility  #####################################
tableWrapper <- function(table) {
  if ('Existing_variation' %in% colnames(table)) {
    write("Existing variation",stderr());
    tic("Existing variation");
    table$Existing_variation=makelinks(table$Existing_variation)
    toc <- toc(); write(paste(toc$msg,toc$toc),stderr());
  }
  if ('Gene_Name' %in% colnames(table)) {
    write("Gene name",stderr());
    tic("Gene name");
    table$Gene_Name[nchar(table$Gene_Name)>200]='many'
    toc <- toc(); write(paste(toc$msg,toc$toc),stderr());
  }
  if ('SYMBOL' %in% colnames(table)) {
    write("SYMBOL",stderr());
    tic("SYMBOL");
    table$SYMBOL=makelinks(table$SYMBOL)
    toc <- toc(); write(paste(toc$msg,toc$toc),stderr());
  }
  if ('cancer_genes' %in% colnames(table)) {
    write("Cancer genes",stderr());
    tic("Cancer genes");
    table$cancer_genes=makelinks(table$cancer_genes,sep=' ')
    toc <- toc(); write(paste(toc$msg,toc$toc),stderr());
  }
  if ('AD_TUMOR' %in% colnames(table)) {
    write("AD_TUMOR",stderr());
    tic("AD_TUMOR");
    table$AD_TUMOR=unlist(lapply(table$AD_TUMOR,paste,collapse=', '))
    toc <- toc(); write(paste(toc$msg,toc$toc),stderr());
  }
  if ('AD_NORMAL' %in% colnames(table)) {
    write("AD_NORMAL",stderr());
    tic("AD_NORMAL");
    table$AD_NORMAL=unlist(lapply(table$AD_NORMAL,paste,collapse=', '))
    toc <- toc(); write(paste(toc$msg,toc$toc),stderr());
  }
  if ('AD' %in% colnames(table)) {
    write("AD",stderr());
    tic("AD");
    if (is.list(table$AD)) table$AD=unlist(lapply(table$AD,paste,collapse=', '))
    toc <- toc(); write(paste(toc$msg,toc$toc),stderr());
  }
  
  write("Sorting table",stderr());
  tic("Sorting")
  tableH = head(table,10)
  samples=sort(unique(tableH$sample))
  toc <- toc(); write(toc$toc,stderr());
   
  write("Creating HTML table",stderr());
  tic("Creating HTML table")
  if (T) htmlTable(tableH, 
            col.rgroup = c("none", "#F7F7F7"),
            tspanner=paste('Rank score:',unique(tableH$rank_score)),
            n.tspanner=rev(table(tableH$rank_score)))
  toc <- toc(); write(toc$toc,stderr());

  # if (length(samples)>1) htmlTable(table, 
  #           col.rgroup = c("none", "#F7F7F7"),
  #           n.rgroup='',
  #           tspanner=paste('Rank score:',unique(table$rank_score)),
  #           n.tspanner=table(table$rank_score))
}


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
  if (!is.na(freec_Tratio_file[1])) {
    # extract sample names
    write("Extract sample names, Tratio and Tbaf",stdout());
    for (s in 1:length(freec_Tratio_file)) samplename[s]=strsplit(basename(freec_Tratio_file[s]),'[.]')[[1]][1]
 
    tratio=NULL; tbaf=NULL
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
    write("Normal log ration and BAF",stdout());
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
}


loadControlFREEC(freec_Tratio_file,freec_Nratio_file, freec_Tbaf_file, freec_Nbaf_file, freec_info_file, freec_cnv_file)
