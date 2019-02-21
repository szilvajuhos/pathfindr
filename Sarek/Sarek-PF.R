#!/usr/bin/env Rscript
# call like:
# Rscript ~/dev/pathfindr/Sarek/Sarek-PF.R -r ~/dev/pathfindr/Sarek/reference_data/ -o bugger -s /data0/btb/P2233/test/
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
    outfile <- paste0(sampleData$name,'_ascat_cnv.csv')
		fwrite(ascat_cnv,file=outfile)
    write(paste0(" *** ASCAT results written to ",outfile),stdout())

		#tableWrapper(ascat_cnv)
	}
}


### Manta (normal)
loadMantaNormal <-function(manta_normal) {
  suppressWarnings(suppressMessages(library(VariantAnnotation)))

  manta_normal_table=NULL
  # first collect PASS ids from all samples
  allpass=NULL
  vcfs=list()
  if (length(manta_normal_file)>0) for (s in 1:length(manta_normal_file)) {
    vcfs[[manta_normal_file[s]]]=readVcf(file = manta_normal_file[s],genome = reference_genome)
    pass=rowRanges(vcfs[[manta_normal_file[s]]])$FILTER=='PASS'
    allpass=c(allpass,names(vcfs[[manta_normal_file[s]]])[pass])
  }
  # then collect variants...
  if (length(manta_normal_file)>0) for (s in 1:length(manta_normal_file)) {
    sample=strsplit(basename(manta_normal_file[s]),'[.]')[[1]][1]
    cat(manta_normal_file[s],'\n')
    
    vcf=vcfs[[manta_normal_file[s]]]
    vcf=vcf[names(vcf) %in% allpass]

    if (length(vcf)>0) {
      
      
      g=geno(vcf)
      inf=info(vcf)
      rr=rowRanges(vcf)
      
      # Read counts
      srtable=as.data.table(g$SR)
      colnames(srtable)=paste0('SplitReadSupport_',c('N')) # Verified on a test sample
      prtable=as.data.table(g$PR)
      colnames(prtable)=paste0('PairedReadSupport_',c('N'))
      
      # Calculate allele ratio
      temp=data.table(
        pr.ref=sapply(prtable$PairedReadSupport_N,"[",1),
        sr.ref=sapply(srtable$SplitReadSupport_N,"[",1),
        pr.alt=sapply(prtable$PairedReadSupport_N,"[",2),
        sr.alt=sapply(srtable$SplitReadSupport_N,"[",2)
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
                   LOH=''),
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
      # loop through all and extract endpoint chr and pos   <------ This one must be remade..
      for (i in 1:nrow(sv_table)) try( {                    # <----  sometimes error here, so try..
        t=strsplit(x = sv_table$ALT[[i]],split = ':')[[1]]
        if (length(t)>1 & t[1]!="<DUP") {
          tchr=str_extract(t[1],'[0-9,X,Y]*$')
          #if (is.na(tchr)) tchr=strsplit(x = t[1],split = 'CHR')[[1]][2]
          #if (project=='BTB') 
          #   sv_table$altchr[i] <- paste0('chr',tchr) 
          #else 
          #   sv_table$altchr[i] <- tchr
          sv_table$altchr[i] <- paste0('chr',tchr) 
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
      snpEff_table=unique(as.data.table(snpEff_table))[Annotation_Impact=='HIGH']  # <--- only HIGH impact kept for the normal
      
      # Filter out annotations of certain types where the ID has >N annotations of that type
      ids=unique(snpEff_table$ID)
      for (id in ids) {
        common=c('protein_protein_contact', 'duplication', 'structural_interaction_variant', 'inversion',
                 'transcript_ablation', 'feature_ablation', 'sequence_feature', 
                 'intergenic_region', 'downstream_gene_variant', 'upstream_gene_variant')
        annotations=snpEff_table$Annotation[snpEff_table$ID==id] # the annotations of this variant
        table=table(annotations[annotations %in% common])
        table=table[table>20] # the common annotations that appear >N times for this variant
        if (length(table)>0) {
          #browser()
          remove=which(snpEff_table$ID==id & snpEff_table$Annotation %in% names(table))
          snpEff_table=snpEff_table[-remove[-1],] # saves one to make sure the variant has some annotation
        }
      }
      
      # Add to data (all samples, one table)
      manta_normal_table=rbind(manta_normal_table,merge(sv_table,snpEff_table,by='ID',all=F))
      setkey(manta_normal_table,'sample')
    }
  }# done parsing each sample
  
  
  # Prepare ranking and (some more) filtering
  if (!is.null(manta_normal_table)) if (nrow(manta_normal_table)>0) {
    
    # ## Make table with most important info for ranking, and report
    # selected <- unique(manta_normal_table[,.(ID,sample,SVTYPE,chr,start,REF,ALT,AFreq,PairedReadSupport_T,SplitReadSupport_T,
    #                        Swegen_count,Rank_score='',Rank_terms='',Gene_Name,Annotation,Annotation_Impact)])
    selection=manta_normal_table[!is.na(Annotation)]
    
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
      
      # Add Control Freec LOH.
      try( {
        if (nrow(selection)>0 & !is.null(freec_loh)) for (i in 1:nrow(selection)) {
          ix=freec_loh$chr==selection$chr[i] &
            freec_loh$start<selection$start[i] &
            freec_loh$end>selection$end[i]
          ix=ix[!is.na(ix)]
          if (sum(ix) >0) selection$LOH[i]='Y'
        }
      },silent=T)
      
      firstcols=c('ID','sample','Gene_Name','rank_score','rank_terms','LOH','AFreq','Annotation','Annotation_Impact','Swegen_count')
      cols=colnames(selection)
      setcolorder(x = selection,neworder = c(firstcols,cols[!cols %in% firstcols]))
      
      
      manta_normal_selected <- selection[order(Feature_ID)][order(rank_score,decreasing = T)]
      outfile <- paste0(sampleData$name,'_manta_normal.csv')
      fwrite(manta_normal_selected[,-c('ANN','Gene_ID')][rank_score>1],file=outfile)
      write(paste0(" *** Manta normal features written to ",outfile),stdout())
      #tableWrapper(manta_normal_selected[,-c('ANN','Gene_ID')][rank_score>1])
    } 
  }
}

loadMantaTumor <- function(manta_tumor_file) {
  suppressWarnings(suppressMessages(library(VariantAnnotation)))

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
      outfile <- paste0(sampleData$name,'_manta_tumor.csv')
      fwrite(manta_tumor_selected[,-c('ANN','Gene_ID')][rank_score>1],file=outfile)
      write(paste0(" *** Manta results written to ",outfile),stdout())
			#tableWrapper(manta_tumor_selected[,-c('ANN','Gene_ID')][rank_score>1])
		} 
	}
}

loadMutect2 <- function(mutect2_file) {
  suppressWarnings(suppressMessages(library(VariantAnnotation)))

  # first collect PASS ids from all samples
  allpass=NULL
  if (exists('mutect2_file')) if (length(mutect2_file)>0) {
    for (s in 1:length(mutect2_file)) {
      vcf=readVcf(file = mutect2_file[s],genome = reference_genome)
      pass=rowRanges(vcf)$FILTER=='PASS'
      allpass=c(allpass,names(vcf)[pass])
    }
    # then collect variants...
    mutect2_table=NULL
    for (s in 1:length(mutect2_file)) {
      sample=strsplit(basename(mutect2_file[s]),'[.]')[[1]][1]
      vcf=readVcf(file = mutect2_file[s],genome = reference_genome)
      vcf=vcf[names(vcf) %in% allpass]

      if (length(vcf)>0) {
        rr=rowRanges(vcf)
        g=geno(vcf)
        inf=info(vcf)	
        # manipulate into a data frame with relevant data
        mutations=as.data.table(ranges(rr))
        mutations$chr <- as.character(seqnames(rr))
        mutations$sample=sample
        mutations=mutations[,.(ID=names,sample,chr,start,end,width)]
    
        mutations$ref=as.character(ref(vcf))
        mutations$alt=as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
        mutations$type=paste0(mutations$ref,'>',mutations$alt)
        mutations$type[nchar(mutations$ref)>nchar(mutations$alt)]='del'
        mutations$type[nchar(mutations$ref)<nchar(mutations$alt)]='ins'
        mutations$type[nchar(mutations$type)>3]='other'
        mutations$type[mutations$type=='T>G']='A>C'
        mutations$type[mutations$type=='T>C']='A>G'
        mutations$type[mutations$type=='T>A']='A>T'
        mutations$type[mutations$type=='G>T']='C>A'
        mutations$type[mutations$type=='G>C']='C>G'
        mutations$type[mutations$type=='G>A']='C>T'
  
        # Headers are sometimes sample names instead of TUMOR-NORMAL
        headers=colnames(g$AD)
        if (!'TUMOR' %in% headers) {
          ix=grep('-02$',headers)
          if (length(ix)==0) 
            ix=grep('[TR]',headers)
          headers[ix]='TUMOR'
          headers[-ix]='NORMAL'
          colnames(g$AD)=headers
          colnames(g$AF)=headers # <-- assumes the same header order in AF
        }
        mutations$AFreq=round(sapply(g$AF[,'TUMOR'], "[", 1),2)
        ad=as.data.table(g$AD)
        colnames(ad)=paste0('AD_',colnames(ad))
        mutations=cbind(mutations,ad)
        # mutations$DP=sapply(g$AD[,'TUMOR'], sum)
        # mutations$AD=sapply(g$AD[,'TUMOR'], "[[", 2)
        # mutations$DPn=sapply(g$AD[,'NORMAL'], sum)
        # mutations$ADn=sapply(g$AD[,'NORMAL'], "[[", 2)
        mutations$reads=''
        temp=sapply(g$AD[,'TUMOR'], "[[", 2)
        mutations$reads[temp<5]='<5'
        mutations$reads[temp>=5]='≥5'
    
        # add counts from the new sources
        suppressWarnings( {
          t=unlist(lapply(X=inf$CAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          mutations$CAF=t
          t=unlist(lapply(X=inf$SWAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          mutations$SWAF=t
          t=unlist(lapply(X=inf$TOPMED,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          mutations$TOPMED=t
        } )
        # Add Swegen counts
        mutations$Swegen_count=snptable[names(vcf),value]
        mutations$Swegen_count[is.na(mutations$Swegen_count)]=0
        # Some have multiple rsIDs (diminishingly few)
        # Also double check Swegen in case pos in reference but rsID or not properly matched in data
        pos=paste0(mutations$chr,':',mutations$start,'_',mutations$ref,'/',mutations$alt)
        snps=snptable[c(pos)][!is.na(value)] # extract from reference[those present]
        if (nrow(snps)>0) 
          mutations$Swegen_count[match(snps$name,pos)]=snps$value # put
      
        ## Annotate by cosmic (for possibly retaining non PASS hotspots, which is not necessary smart)
        key=paste0(substr(mutations$chr,4,6),':',mutations$start,'-',mutations$end)
        key=str_replace(key,'X:','23:')
        key=str_replace(key,'Y:','24:')
      
        counts=cbind(cosmic_coding[key,value],cosmic_noncoding[key,value],0)
        max_=apply(counts,1,max,na.rm=T)
        mutations$Cosmic_count=max_

        # Add Control Freec LOH.
        mutations$rank_score=0
        mutations$rank_terms=''
        mutations$LOH=''
        try( {
          if (nrow(mutations)>0 & !is.null(freec_loh)) 
            for (i in 1:nrow(mutations)) {
              ix = freec_loh$chr==mutations$chr[i] &
                   freec_loh$start<mutations$start[i] &
                   freec_loh$end>mutations$end[i]
              ix = ix[!is.na(ix)]
              if (sum(ix) >0) 
                mutations$LOH[i]='Y'
          }
        },silent=T)
   
        # "info"" annotations also need to be added
        mutations=cbind(mutations,as.data.table(info(vcf))[,.(CSQ)])

        # keep only PASS
        # mutations=mutations[filter=='PASS',]

        # get snpEff headers
        # snpeff_header=strsplit(info(header(vcf))['ANN',][,3],'ions: \'')[[1]][2]
        # snpeff_header=strsplit(snpeff_header,' \\| ')[[1]]
        ## get VEP headers
        vep_header=strsplit(info(header(vcf))['CSQ',][,3],'Format: ')[[1]][2]
        vep_header=strsplit(vep_header,'\\|')[[1]]

        # VEP annotation is put in annotation_table
        annotation_table=matrix(data = NA,nrow = length(unlist(mutations$CSQ)),ncol = length(vep_header)+1)
        colnames(annotation_table)=c('ID',vep_header)
        row=1
        for (i in 1:nrow(mutations)) { # for each variant
          # for each VEP annotation:
          for (j in 1:length(mutations$CSQ[[i]])) {
            line=strsplit(mutations$CSQ[[i]][j],'\\|')[[1]]
            annotation_table[row,1]=mutations$ID[i]
            annotation_table[row,1+(1:length(line))]=line
            row=row+1
          }
        }
        annotation_table=as.data.table(annotation_table)
        # Annotation table then concatenated
        mutect2_table=rbind(mutect2_table,merge(mutations,annotation_table,by='ID'))
      }
    } # done collecting from vcf files

    # filtering by SWAF / Swegen_count (same as with Haplotypecaller)
    mutect2_table=mutect2_table[-which(mutect2_table$SWAF>=0.01)]
    mutect2_table=mutect2_table[-which(mutect2_table$Swegen_count>=10)]
    mutect2_table=mutect2_table[-which(mutect2_table$TOPMED>=0.01)]
  
    mutect2_table$cumstart=NA
    mutect2_table$cumend=NA
    # # for each chr get cumulative pos
    for (i in 1:nrow(chrsz)) {
      ix=mutect2_table$chr==chrsz$chr[i]
      mutect2_table$cumstart[ix]=mutect2_table$start[ix]+chrsz$starts[i]
      mutect2_table$cumend[ix]=mutect2_table$end[ix]+chrsz$starts[i]
    }
    setkey(mutect2_table,'sample')
    selection=mutect2_table[,-c('CSQ')] # not needed after parsing/merging the annotations
    if (nrow(selection)>0) {
      # Cosmic/local Tier2 :
      ix=selection$SYMBOL %in% alltier2
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T2_gene'
      # tier 1 priority:
      ix=selection$SYMBOL %in% alltier1
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T1_gene'

      # Add high impact
      ix=selection$IMPACT=='HIGH'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high_impact')
      }

      # Add moderate impact
      ix=selection$IMPACT=='MODERATE'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'moderate_impact')
      }

      # additional +1 if high impact and TSG
      ix=selection$IMPACT=='HIGH' & selection$SYMBOL %in% alltsg
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high+TSG')
      }

      # Add clinvar pathogenic
      ix=grep('pathogenic',selection$CLIN_SIG)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'clinvar')
      }

      # Add Polyphen/SIFT damaging/deleterious
      ix=union(grep('damaging',selection$PolyPhen),grep('deleterious',selection$SIFT))
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'polyphen/SIFT')
      }

      # Add hotspots and near hotspots (within 2 residues of a hotspot)
      key=paste(selection$SYMBOL,
              str_replace(string = selection$Protein_position,
              pattern = '/.*',replacement = ''))
      ix <-
        key %in% paste(hotspots_snv[,Hugo_Symbol],hotspots_snv[,Amino_Acid_Position]) |
        key %in% paste(hotspots_inframe[,Hugo_Symbol],hotspots_inframe[,Amino_Acid_Position])  ## Warning: Exact match used with inframes
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'hotspot')
      }
      ix = !ix & key %in% near_hotspots # the near_hotspot excludes those that were a hotspot
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'near_hotspot')
      }
    
      # Add cosmic counts
      ix=which(selection$Cosmic_count>50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>50')
      }
      ix=which(selection$Cosmic_count>5 & selection$Cosmic_count<=50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>5')
      }
      
      # Special case for TERT promoter
      ix=which(selection$SYMBOL=='TERT' & selection$Consequence=='5_prime_UTR_variant')
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TERT_5\'UTR')
      }

      # Add TF binding variants near (100kb) cancer genes
      ix=grep('TF',selection$Consequence)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TFBS')
        for (i in 1:length(ix)) {
          near = which(selection$cumend[ix[i]] > (tumorgenes$cumstart-100e3) & 
                 selection$cumstart[ix[i]] < (tumorgenes$cumend + 100e3) &
                 !selection$SYMBOL[ix[i]] %in% alltumorgenes)
          if (length(near)>0) {
            genes=paste(tumorgenes$`Gene Symbol`[near],collapse = ',')
            selection$rank_score[ix[i]]=selection$rank_score[ix[i]]+2
            selection$rank_terms[ix[i]]=paste(selection$rank_terms[ix[i]],paste0('near_',genes))
          }
        }
      }

      ix=which(selection$CANONICAL!='YES')
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]-0
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'not_canonical')
      }
    }

    firstcols=c('ID','sample','SYMBOL','rank_score','rank_terms','LOH','AFreq','Consequence','IMPACT','SWAF','TOPMED','Swegen_count')
    cols=colnames(selection)
    setcolorder(x = selection,neworder = c(firstcols,cols[!cols %in% firstcols]))
    selection <- selection[order(cumstart,Allele)][order(rank_score,decreasing = T)]
  
    # # Take some information to the VCF
    # selected_for_vcf=selection[, .(LOH=LOH,Swegen_count=Swegen_count,rank_terms=paste0(rank_score,':',SYMBOL,' ',trimws(rank_terms)),rank_score=rank_score), by = ID]
    # selected_for_vcf=selected_for_vcf[, .(LOH=LOH,Swegen_count=Swegen_count,rank_terms=paste(unique(rank_terms),collapse=','),rank_score=paste(unique(rank_score),collapse=',')), by=ID]
    # selected_for_vcf$rank_terms[selected_for_vcf$rank_score=='0']=''
  
    mutect2_selected <- selection
    outfile <- paste0(sampleData$name,'_mutect2_tumor.csv')
    fwrite(mutect2_selected,file=outfile)
    write(paste0(" *** MuTect2 GATK v3.8 results written to ",outfile),stdout())
  }
}

loadStrelka <- function() {
  suppressWarnings(suppressMessages(library(VariantAnnotation)))

  allpass=NULL
  vcfs=list()
  if (length(strelka_snv_file)>0  & length(strelka_indel_file)>0) for (s in 1:length(strelka_snv_file)) {
    # read snvs
    vcfs[[strelka_snv_file[s]]]=readVcf(file = strelka_snv_file[s],genome=reference_genome)
    pass=rowRanges(vcfs[[strelka_snv_file[s]]])$FILTER=='PASS'
    allpass=c(allpass,names(vcfs[[strelka_snv_file[s]]])[pass])
    # read indels
    vcfs[[strelka_indel_file[s]]]=readVcf(file = strelka_indel_file[s],genome=reference_genome)
    pass=rowRanges(vcfs[[strelka_indel_file[s]]])$FILTER=='PASS'
    allpass=c(allpass,names(vcfs[[strelka_indel_file[s]]])[pass])
  }
  strelka_table=NULL
  write("Creating Strelka table",stderr())
  if (length(strelka_snv_file)>0  & length(strelka_indel_file)>0) {
    write("processing Strelka SNVs",stderr());
    for (s in 1:length(strelka_snv_file)) {
      sample=strsplit(basename(strelka_snv_file[s]),'_somatic')[[1]][1]
      # first snvs
      vcf=vcfs[[strelka_snv_file[s]]]
      vcf=vcf[names(vcf) %in% allpass] # only those with PASS in either sample
      write(paste("Changing TUMOR and NORMAL", colnames(vcf)),stderr());
      if(!'TUMOR' %in% colnames(vcf)) 
        colnames(vcf)=c('NORMAL','TUMOR') # assumption that normal comes first
      # Collect sample-unspecific data for the [PASS in any] IDs into a table
      write(paste("colnames",colnames(vcf)), stderr());
      g=geno(vcf)
      rr=rowRanges(vcf)
      inf=info(vcf)
      table_snvs=NULL
      csq=NULL
      if (length(vcf)>0) {
        table_snvs = data.table(ID=as.character(names(vcf)),
                            sample=sample,
                            chr=as.character(seqnames(rr)),
                            start=start(rr),
                            end=end(rr),
                            ref=as.data.table(rr$REF)$x,
                            alt=as.data.table(rr$ALT)$value,
                            AFreq=NA,AD=NA,DP=as.data.table(g$DP)$TUMOR,
                            AD_normal=NA,DP_normal=as.data.table(g$DP)$NORMAL
                     )
        # Collect variant allele depths for SNVs:
        for (this_ref in unique(table_snvs$ref)) {
          for (this_alt in unique(table_snvs[ref==this_ref,alt])) {
            which_ones = table_snvs$ref == this_ref & table_snvs$alt==this_alt
            refcounts=as.data.table(data.frame(g[[paste0(this_ref,'U')]]))
            altcounts=as.data.table(data.frame(g[[paste0(this_alt,'U')]]))
            table_snvs$AFreq[which_ones]=round(altcounts[which_ones,TUMOR.1] /
                                             (altcounts[which_ones,TUMOR.1]+refcounts[which_ones,TUMOR.1]),2)
            table_snvs$AD[which_ones]=altcounts[which_ones,TUMOR.1]
            table_snvs$AD_normal[which_ones]=altcounts[which_ones,NORMAL.1]
          }
        }
        table_snvs$reads='≥5'
        table_snvs$reads[table_snvs$AD<5]='<5'
        # CAF,SWAF,TOPMED
        #if (project=='BTB') 
        suppressWarnings({
          t=unlist(lapply(X=inf$CAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          table_snvs$CAF=t
          t=unlist(lapply(X=inf$SWAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          table_snvs$SWAF=t
          t=unlist(lapply(X=inf$TOPMED,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          table_snvs$TOPMED=t
        } )
      
        # The snv annotations will be used later:
        csq=as.data.table(info(vcf))[,CSQ]
      }
    
      # Same procedure for indels
      vcf=vcfs[[strelka_indel_file[s]]]
      vcf=vcf[names(vcf) %in% allpass] # only those with PASS in either sample
      if (!'TUMOR' %in% colnames(vcf)) 
        colnames(vcf)=c('NORMAL','TUMOR') # assumption that normal comes first
      # Collect sample-unspecific data for the [PASS in any] IDs into a table
      g=geno(vcf)
      rr=rowRanges(vcf)
      inf=info(vcf)
      table_indels=NULL
      if (length(vcf)>0) {
        table_indels=data.table(ID=as.character(names(vcf)),
                              sample=sample,
                              chr=as.character(seqnames(rr)),
                              start=start(rr),
                              end=end(rr),
                              ref=as.data.table(rr$REF)$x,
                              alt=as.data.table(rr$ALT)$value,
                              AFreq=NA,AD=NA,DP=as.data.table(g$DP)$TUMOR,
                              AD_normal=NA,DP_normal=as.data.table(g$DP)$NORMAL
                     )
        refcounts=as.data.table(data.frame(g$TAR))
        altcounts=as.data.table(data.frame(g$TIR))
        table_indels$AFreq=round(altcounts$TUMOR.1 / (altcounts$TUMOR.1 + refcounts$TUMOR.1),2)
        table_indels$AD_normal=altcounts$NORMAL.1
        table_indels$AD=altcounts$TUMOR.1
        table_indels$reads='≥5'
        table_indels$reads[table_indels$AD<5]='<5'
        # CAF,SWAF,TOPMED
        #if (project=='BTB') 
        suppressWarnings({
          t=unlist(lapply(X=inf$CAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          table_indels$CAF=t
          t=unlist(lapply(X=inf$SWAF,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          table_indels$SWAF=t
          t=unlist(lapply(X=inf$TOPMED,FUN=max,na.rm=T))
          t[is.infinite(t)]=0
          table_indels$TOPMED=t
        })
      
        # indel annotations added to snv annotations
        csq=c(csq,as.data.table(info(vcf))[,CSQ])
      }

      # Concatenate the snvs and indels tables for this patient and sample
      table_both=rbind(table_snvs,table_indels)
      table_both$Swegen_count=NA
      table_both$Cosmic_count=NA
      table_both$rank_score=0
      table_both$rank_terms=''
      table_both$LOH=''
 
      # Variant type
      table_both$type=paste0(table_both$ref,'>',table_both$alt)
      table_both$type[nchar(table_both$ref)>nchar(table_both$alt)]='del'
      table_both$type[nchar(table_both$ref)<nchar(table_both$alt)]='ins'
      table_both$type[nchar(table_both$type)>3]='other'
      table_both$type[table_both$type=='T>G']='A>C'
      table_both$type[table_both$type=='T>C']='A>G'
      table_both$type[table_both$type=='T>A']='A>T'
      table_both$type[table_both$type=='G>T']='C>A'
      table_both$type[table_both$type=='G>C']='C>G'
      table_both$type[table_both$type=='G>A']='C>T'
 
      table_both$CSQ=csq
 
      # get VEP headers from a vcf object (assumed never to differ between snvs and indels)
      vep_header=strsplit(info(header(vcf))['CSQ',][,3],'Format: ')[[1]][2]
      vep_header=strsplit(vep_header,'\\|')[[1]]
 
      # VEP annotation is put in annotation_table
      annotation_table=matrix(data = NA,nrow = length(unlist(csq)),ncol = length(vep_header)+1)
      colnames(annotation_table)=c('ID',vep_header)
      row=1
      for (i in 1:length(csq)) { # for each variant
      # for each VEP annotation:
        for (j in 1:length(csq[[i]])) {
          line=strsplit(csq[[i]][j],'\\|')[[1]]
          annotation_table[row,1]=table_both$ID[i]
          annotation_table[row,1+(1:length(line))]=line
          row=row+1
        }
      }

      # Annotation table then merged with the "table_both"
      strelka_table=rbind(strelka_table,merge(table_both,as.data.table(annotation_table),by='ID'))
      setkey(strelka_table,'sample')
    } # done parsing each sample
    rm(vcfs)

    # Add Control Freec LOH.  <-- slow.
    try( {
      if (nrow(strelka_table)>0 & !is.null(freec_loh)) 
        for (i in 1:nrow(strelka_table)) {
          ix = freec_loh$chr==strelka_table$chr[i] &
               freec_loh$start<strelka_table$start[i] &
               freec_loh$end>strelka_table$end[i]
          ix=ix[!is.na(ix)]
          if (sum(ix) >0) 
            strelka_table$LOH[i]='Y'
        }
     },silent=T)

    # Add Swegen counts
    by_pos=snptable[strelka_table$ID,value] # if the variant is as 10:10001801_C/T in database
    by_name1=snptable[str_extract(strelka_table$Existing_variation,"^rs[0-9]+"),value] # if first rsid
    by_name2=snptable[str_extract(strelka_table$Existing_variation,"rs[0-9]+$"),value] # if last rsid (often both)
    d=data.table(by_pos,by_name1,by_name2,default=0)
    d$max=apply(X = d,MARGIN = 1,FUN = max,na.rm=T)
    strelka_table$Swegen_count=d$max
 
    # filtering by SWAF / Swegen_count (same as with Haplotypecaller)
    strelka_table=strelka_table[-which(strelka_table$SWAF>=0.01)]
    strelka_table=strelka_table[-which(strelka_table$Swegen_count>=10)]
    strelka_table=strelka_table[-which(strelka_table$TOPMED>=0.01)]
    
    # # Add Cosmic counts
    key=paste0(substr(strelka_table$chr,4,6),':',strelka_table$start,'-',strelka_table$end)
    #if (project!='BTB') 
    key=paste0(strelka_table$chr,':',strelka_table$start,'-',strelka_table$end)
    key=str_replace(key,'X:','23:')
    key=str_replace(key,'Y:','24:')
    counts=cbind(cosmic_coding[key,value],cosmic_noncoding[key,value],0)
    max_=apply(counts,1,max,na.rm=T)
    strelka_table$Cosmic_count=max_
 
    strelka_table$cumstart=strelka_table$start
    strelka_table$cumend=strelka_table$end
    # # for each chr get cumulative pos
    for (i in 1:nrow(chrsz)) {
      ix=strelka_table$chr==chrsz$chr[i]
      strelka_table$cumstart[ix]=strelka_table$start[ix]+chrsz$starts[i]
      strelka_table$cumend[ix]=strelka_table$end[ix]+chrsz$starts[i]
    }
    selection=strelka_table[,-c('CSQ')] # not needed after parsing/merging the annotations

    if (nrow(selection)>0) {
      # Cosmic/local Tier2 :
      ix=selection$SYMBOL %in% alltier2
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T2_gene'
      # tier 1 priority:
      ix=selection$SYMBOL %in% alltier1
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T1_gene'
 
      # Add high impact
      ix=selection$IMPACT=='HIGH'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high_impact')
      }
 
      # Add moderate impact
      ix=selection$IMPACT=='MODERATE'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'moderate_impact')
      }
 
      # additional +1 if high impact and TSG
      ix=selection$IMPACT=='HIGH' & selection$SYMBOL %in% alltsg
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high+TSG')
      }
 
      # Add clinvar pathogenic
      ix=grep('pathogenic',selection$CLIN_SIG)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'clinvar')
      }
 
      # Add Polyphen/SIFT damaging/deleterious
      ix=union(grep('damaging',selection$PolyPhen),grep('deleterious',selection$SIFT))
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'polyphen/SIFT')
      }
 
      # Add hotspots and near hotspots (within 2 residues of a hotspot)
      key=paste(selection$SYMBOL,
                str_replace(string = selection$Protein_position,
                pattern = '/.*',replacement = ''))    
      ix <-
        key %in% paste(hotspots_snv[,Hugo_Symbol],hotspots_snv[,Amino_Acid_Position]) |
        key %in% paste(hotspots_inframe[,Hugo_Symbol],hotspots_inframe[,Amino_Acid_Position])  ## Warning: Exact match used with inframes
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'hotspot')
      }
      ix = !ix & key %in% near_hotspots # the near_hotspot excludes those that were a hotspot
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'near_hotspot')
      }

      # Add cosmic counts
      ix=which(selection$Cosmic_count>50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>50')
      }
      ix=which(selection$Cosmic_count>5 & selection$Cosmic_count<=50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>5')
      }
     
      # Special case for TERT promoter
      ix=which(selection$SYMBOL=='TERT' & selection$Consequence=='5_prime_UTR_variant')
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TERT_5\'UTR')
      }

      # Add TF binding variants near (100kb) cancer genes
      ix=grep('TF',selection$Consequence)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TFBS')
        for (i in 1:length(ix)) {
          near = which(selection$cumend[ix[i]] > (tumorgenes$cumstart-100e3) & 
                       selection$cumstart[ix[i]] < (tumorgenes$cumend + 100e3) &
                       !selection$SYMBOL[ix[i]] %in% alltumorgenes)
          if (length(near)>0) {
            genes=paste(tumorgenes$`Gene Symbol`[near],collapse = ',')
            selection$rank_score[ix[i]]=selection$rank_score[ix[i]]+2
            selection$rank_terms[ix[i]]=paste(selection$rank_terms[ix[i]],paste0('near_',genes))
          }
        }
      }
      
      ix=which(selection$CANONICAL!='YES')
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]-0
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'not_canonical')
      }
    }

    firstcols=c('ID','sample','SYMBOL','rank_score','rank_terms','LOH','AFreq','Consequence','IMPACT','SWAF','TOPMED','Swegen_count')
    cols=colnames(selection)
    setcolorder(x = selection,neworder = c(firstcols,cols[!cols %in% firstcols]))
  
    strelka_selected <- selection[order(cumstart,Allele)][order(rank_score,decreasing = T)]
    outfile <- paste0(sampleData$name,'_strelka_tumor.csv')
    fwrite(strelka_selected,file=outfile)
    write(paste0(" *** Strelka results written to ",outfile),stdout())
  }
}

loadHaplotypecaller <- function() {
  suppressWarnings(suppressMessages(library(VariantAnnotation)))

  if (exists('haplotypecaller_N_file')) {
    haplotypecaller_files=haplotypecaller_N_file
    if (exists('haplotypecaller_T_file')) haplotypecaller_files=c(haplotypecaller_N_file,haplotypecaller_T_file)
    haplotypecaller_table=NULL
    haplotypecaller_ids=NULL
    for (s in 1:length(haplotypecaller_files)) {
      
      vcf=readVcf(file = haplotypecaller_files[s],genome=reference_genome)
      sample=strsplit(basename(haplotypecaller_files[s]),'[.]')[[1]][1]
      haplotypecaller_ids[[sample]]=names(vcf)
      
      if (s>1) # if not the normal, keep only IDs already present (the normal is first)
        vcf=vcf[names(vcf) %in% haplotypecaller_table$ID]
      
      # pre-filtering by SWAF / TOPMED to avoid extreme number of variants
      swaf=as.data.table(info(vcf)$SWAF)
      keep_swaf=unique(swaf$group[which(swaf$value<0.01)])
      topmed=as.data.table(info(vcf)$TOPMED)
      keep_topmed=unique(topmed$group[which(topmed$value<0.01)])
      vcf=vcf[union(keep_swaf,keep_topmed)]
      
      g=geno(vcf)
      inf=(info(vcf))
      rr=rowRanges(vcf)
      
      mutations=as.data.table(ranges(rr))
      mutations$chr <- as.character(seqnames(rr))
      mutations$sample=sample
      mutations=mutations[,.(ID=names,sample,chr,start,end,width)]
      
      mutations$ref=as.character(ref(vcf))
      mutations$alt=as.data.table(alt(vcf))[, .(values = list(value)), by = group][,values]
      mutations$type=paste0(mutations$ref,'>',mutations$alt)
      mutations$type[nchar(mutations$ref)>nchar(mutations$alt)]='del'
      mutations$type[nchar(mutations$ref)<nchar(mutations$alt)]='ins'
      mutations$type[nchar(mutations$type)>3]='other'
      mutations$type[mutations$type=='T>G']='A>C'
      mutations$type[mutations$type=='T>C']='A>G'
      mutations$type[mutations$type=='T>A']='A>T'
      mutations$type[mutations$type=='G>T']='C>A'
      mutations$type[mutations$type=='G>C']='C>G'
      mutations$type[mutations$type=='G>A']='C>T'
      
      mutations$AFreq=round(sapply(g$AD, "[[", 2) / unlist(g$DP),2) # not a perfect estimate...
      ad=as.data.table(g$AD)
      colnames(ad)='AD'
      mutations=cbind(mutations,ad)
      # mutations$reads=''
      # temp=sapply(g$AD, "[[", 2)
      # mutations$reads[temp<5]='<5'
      # mutations$reads[temp>=5]='≥5'
      #
      # add counts from the new sources
      #if (project=='BTB') 
      suppressWarnings( {
        t=unlist(lapply(X=inf$CAF,FUN=max,na.rm=T))
        t[is.infinite(t)]=0
        mutations$CAF=t
        t=unlist(lapply(X=inf$SWAF,FUN=max,na.rm=T))
        t[is.infinite(t)]=0
        mutations$SWAF=t
        t=unlist(lapply(X=inf$TOPMED,FUN=max,na.rm=T))
        t[is.infinite(t)]=0
        mutations$TOPMED=t
      } )
      
      # "info"" annotations also need to be added
      mutations=cbind(mutations,as.data.table(info(vcf))[,.(CSQ)])
      
      # Add Swegen counts
      mutations$Swegen_count=snptable[names(vcf),value]
      mutations$Swegen_count[is.na(mutations$Swegen_count)]=0
      # Some have multiple rsIDs
      ix=grep(';',mutations$ID)
      ids=data.table(ix,id=strsplit(mutations$ID[ix],';'))
      temp_snptable=snptable[name %in% unlist(ids$id)]
      ids$counts=unlist(lapply(ids$id,function(x) return(max(temp_snptable[x,value],na.rm=T))))
      ids$counts[is.infinite(ids$counts)]=0
      mutations$Swegen_count[ix]=ids$counts
      
      # Also double check Swegen in case pos in reference but rsID or not properly matched in data
      pos=paste0(mutations$chr,':',mutations$start,'_',mutations$ref,'/',mutations$alt)
      snps=snptable[c(pos)][!is.na(value)] # extract from reference[those present]
      if (nrow(snps)>0) mutations$Swegen_count[match(snps$name,pos)]=snps$value # put
      # # Some have multiple alt alleles
      # ix=grep(',',mutations$alt)
      # alts=data.table(ix,alt=mutations$alt[ix])
      # alts$counts[ix]=unlist(lapply(alts$pos,function(x) return(max(snptable[x,value],na.rm=T))))
      
      # More thorough filtering
      if (any(mutations$Swegen_count>=10)) mutations=mutations[-which(Swegen_count>=10)]
      
      
      ## Annotate by cosmic (for possibly retaining non PASS hotspots, which is not necessarily smart)
      key=paste0(substr(mutations$chr,4,6),':',mutations$start,'-',mutations$end)
      key=str_replace(key,'X:','23:')
      key=str_replace(key,'Y:','24:')
      
      #if (project!='BTB') 
      key=paste0(mutations$chr,':',mutations$start,'-',mutations$end)
      counts=cbind(cosmic_coding[key,value],cosmic_noncoding[key,value],0)
      max_=apply(counts,1,max,na.rm=T)
      mutations$Cosmic_count=max_
      
      mutations$rank_score=0
      mutations$rank_terms=''
      mutations$LOH=''
      
      # Add Control Freec LOH.
      try( {
        if (nrow(mutations)>0 & !is.null(freec_loh)) for (i in 1:nrow(freec_loh)) {
          ix = freec_loh$chr[i]==mutations$chr &
               freec_loh$start[i]<mutations$start &
               freec_loh$end[i]>mutations$end
          ix=ix[!is.na(ix)]
          if (sum(ix) >0) 
            mutations$LOH[ix]='Y'
        }
      },silent=T)
      
      vep_header=strsplit(info(header(vcf))['CSQ',][,3],'Format: ')[[1]][2]
      vep_header=strsplit(vep_header,'\\|')[[1]]
      
      # VEP annotation is put in annotation_table
      annotation_table=matrix(data = NA,nrow = length(unlist(mutations$CSQ)),ncol = length(vep_header)+1)
      colnames(annotation_table)=c('ID',vep_header)
      row=1
      for (i in 1:nrow(mutations)) { # for each variant
        # for each VEP annotation:
        for (j in 1:length(mutations$CSQ[[i]])) {
          line=strsplit(mutations$CSQ[[i]][j],'\\|')[[1]]
          annotation_table[row,1]=mutations$ID[i]
          annotation_table[row,1+(1:length(line))]=line
          row=row+1
        }
      }
      annotation_table=as.data.table(annotation_table)
      
      # Annotation table then concatenated
      haplotypecaller_table=rbind(haplotypecaller_table,merge(mutations,annotation_table,by='ID'))
    } # done collecting from vcf files
    
    haplotypecaller_table$cumstart=haplotypecaller_table$start
    haplotypecaller_table$cumend=haplotypecaller_table$end
    # # for each chr get cumulative pos
    for (i in 1:nrow(chrsz)) {
      ix=haplotypecaller_table$chr==chrsz$chr[i]
      haplotypecaller_table$cumstart[ix]=haplotypecaller_table$start[ix]+chrsz$starts[i]
      haplotypecaller_table$cumend[ix]=haplotypecaller_table$end[ix]+chrsz$starts[i]
    }
    setkey(haplotypecaller_table,'sample')
    
    # Due to size: remove intron/intergenic variants
    selection=haplotypecaller_table[!Consequence %in% c('intron_variant','intergenic_variant'),-c('CSQ')] # not needed after parsing/merging the annotations
    
    if (nrow(selection)>0) {
      # Cosmic/local Tier2 :
      ix=selection$SYMBOL %in% alltier2
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T2_gene'
      # tier 1 priority:
      ix=selection$SYMBOL %in% alltier1
      selection$rank_score[ix]=2
      selection$rank_terms[ix]='T1_gene'
      
      # Add high impact
      ix=selection$IMPACT=='HIGH'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high_impact')
      }
      
      # Add moderate impact
      ix=selection$IMPACT=='MODERATE'
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'moderate_impact')
      }
      
      # additional +1 if high impact and TSG
      ix=selection$IMPACT=='HIGH' & selection$SYMBOL %in% alltsg
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'high+TSG')
      }
      
      # Add clinvar pathogenic
      ix=grep('pathogenic',selection$CLIN_SIG)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'clinvar')
      }
      
      # Add Polyphen/SIFT damaging/deleterious
      ix=union(grep('damaging',selection$PolyPhen),grep('deleterious',selection$SIFT))
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'polyphen/SIFT')
      }
      
      # Add hotspots and near hotspots (within 2 residues of a hotspot)
      key=paste(selection$SYMBOL,
                str_replace(string = selection$Protein_position,
                            pattern = '/.*',replacement = ''))
      ix <-
        key %in% paste(hotspots_snv[,Hugo_Symbol],hotspots_snv[,Amino_Acid_Position]) |
        key %in% paste(hotspots_inframe[,Hugo_Symbol],hotspots_inframe[,Amino_Acid_Position])  ## Warning: Exact match used with inframes
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'hotspot')
      }
      ix = !ix & key %in% near_hotspots # the near_hotspot excludes those that were a hotspot
      if (any(ix)) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'near_hotspot')
      }
      # Add cosmic counts
      ix=which(selection$Cosmic_count>50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>50')
      }
      ix=which(selection$Cosmic_count>5 & selection$Cosmic_count<=50)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+1
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'cosmic_>5')
      }
      # Special case for TERT promoter
      ix=which(selection$SYMBOL=='TERT' & selection$Consequence=='5_prime_UTR_variant')
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TERT_5\'UTR')
      }
      # Add TF binding variants near (100kb) cancer genes
      ix=grep('TF',selection$Consequence)
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]+2
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'TFBS')
        for (i in 1:length(ix)) {
          near=which(selection$cumend[ix[i]] > (tumorgenes$cumstart-100e3) & selection$cumstart[ix[i]] < (tumorgenes$cumend + 100e3) &
                       !selection$SYMBOL[ix[i]] %in% alltumorgenes)
          if (length(near)>0) {
            genes=paste(tumorgenes$`Gene Symbol`[near],collapse = ',')
            selection$rank_score[ix[i]]=selection$rank_score[ix[i]]+2
            selection$rank_terms[ix[i]]=paste(selection$rank_terms[ix[i]],paste0('near_',genes))
          }
        }
      }
      ix=which(selection$CANONICAL!='YES')
      if (length(ix)>0) {
        selection$rank_score[ix]=selection$rank_score[ix]-0
        selection$rank_terms[ix]=paste(selection$rank_terms[ix],'not_canonical')
      }
    }
    firstcols=c('ID','sample','SYMBOL','rank_score','rank_terms','LOH','AFreq','Consequence','IMPACT','SWAF','TOPMED','Swegen_count')
    cols=colnames(selection)
    setcolorder(x = selection,neworder = c(firstcols,cols[!cols %in% firstcols]))
    
    haplotypecaller_selected <- selection[order(cumstart,Allele)][order(rank_score,decreasing = T)]
    outfile <- paste0(sampleData$name,'_haplotypecaller.csv')
		fwrite(haplotypecaller_selected,file=outfile)
    write(paste0(" *** HaploTyper results written to ",outfile),stdout())
    #ix=haplotypecaller_selected$sample==strsplit(basename(haplotypecaller_files[1]),'[.]')[[1]][1] # first file is the normal
    #tableWrapper(haplotypecaller_selected[ix][,-c('cumstart','cumend','DOMAINS')][rank_score>1])
  }
}

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
loadControlFREEC(freec_Tratio_file,freec_Nratio_file, freec_Tbaf_file, freec_Nbaf_file, freec_info_file, freec_cnv_file)
toc()

tic("ASCAT")
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
loadASCAT(ascat_Tratio_file,ascat_Nratio_file,ascat_Tbaf_file,ascat_Nbaf_file,ascat_segment_file)
toc()

tic("Manta")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             Manta                                                             x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
manta_tumor_file <- grep(pattern = ".*Manta_.*vs.*somaticSV.vcf.snpEff.ann.vep.ann.vcf$",files,value = T)
manta_normal_file <- grep(pattern = ".*Manta_.*vs.*diploidSV.vcf.snpEff.ann.vep.ann.vcf$",files,value = T)[1]
cosmic_fusions = fread(paste0(ref_data,'cosmic_fusions_table.csv'),key = 'name')
swegen_manta_all=fread(paste0(ref_data,'swegen_sv_counts.csv'),key='name')
write("Files for Manta structural variants: ", stdout())
printList( list(manta_tumor_file,manta_normal_file) )
# TODO: sort it out to be local or a general on-demand object
allfusionpairs=NULL
allfusion=tumorgenes[grep('fusion',`Role in Cancer`),`Gene Symbol`]

for (i in 1:length(allfusion)) {
  t=trimws(strsplit(tumorgenes[allfusion[i],`Translocation Partner`],', ')[[1]])
  if (length(t)>0) 
		for (j in 1:length(t))
 		allfusionpairs=c(allfusionpairs,paste(sort(c(allfusion[i],t[j])),collapse = ' '))
}
loadMantaTumor(manta_tumor_file)
loadMantaNormal(manta_normal_file)
toc()

tic("Database read")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                     Reading databases for SNV and small indel prioritisation                  x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
# These are also needed for Strelka
write("Reading SweGen SNP counts ...",stdout())
snptable         = fread(paste0(ref_data,'swegen_snp_counts.csv'),key='name')
write("Reading COSMIC tables ...",stdout())
cosmic_coding    = fread(paste0(ref_data,'cosmic_coding_table.csv'),key = 'name')
cosmic_noncoding = fread(paste0(ref_data,'cosmic_noncoding_table.csv'),key = 'name')
cosmic_fusions   = fread(paste0(ref_data,'cosmic_fusions_table.csv'),key = 'name')
alltsg = tumorgenes[grep('TSG',`Role in Cancer`),`Gene Symbol`]
write("Reading hotspots ...",stdout())
hotspots_snv = unique( fread(paste0(ref_data,'hotspots_v2_snv.csv'))[,.(Hugo_Symbol,Amino_Acid_Position)])[-grep('splice',Amino_Acid_Position)]
hotspots_snv$pos <- as.numeric(hotspots_snv$Amino_Acid_Position)
hotspots_inframe = unique(fread(paste0(ref_data,'hotspots_v2_inframe.csv'))[,.(Hugo_Symbol,Amino_Acid_Position)])
hotspots_inframe$start=as.numeric(str_replace(string = hotspots_inframe$Amino_Acid_Position,pattern = '-[0-9]*',replacement = ''))
hotspots_inframe$end=as.numeric(str_replace(string = hotspots_inframe$Amino_Acid_Position,pattern = '[0-9]*-',replacement = ''))
near_hotspots=NULL
for (i in 1:nrow(hotspots_inframe)) 
  near_hotspots = c(near_hotspots,
                    paste(hotspots_inframe$Hugo_Symbol[i],
                    c(seq(hotspots_inframe$start[i]-2,hotspots_inframe$start[i]+2),
                    hotspots_inframe$end[i]-2,hotspots_inframe$end[i]+2))
                  )
for (i in 1:nrow(hotspots_snv))
  near_hotspots = c(near_hotspots,
                    paste(hotspots_snv$Hugo_Symbol[i],
                    seq(hotspots_snv$pos[i]-2,hotspots_snv$pos[i]+2))
                  )

toc()
tic("Mutect2")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             Mutect2 with GATK 3.8                                             x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("Reading MuTect2 file ...",stdout())
mutect2_file <- grep(pattern = ".*mutect2_.*AF.*vep.ann.vcf$",files,value = T)

loadMutect2(mutect2_file)
toc()

tic("Strelka")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             Strelka                                                           x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
# snpEff + VEP Strelka snv file
strelka_snv_file <- grep(pattern = ".*Strelka_.*_somatic_snvs.*AF.*vep.ann.vcf$",files,value = T)
# snpEff + VEP Strelka indel file
strelka_indel_file <- grep(pattern = ".*Strelka_.*_somatic_indels.*AF.*vep.ann.vcf$",files,value = T)
write(paste("Strelka SNV",strelka_snv_file),stderr());
write(paste("Strelka indel",strelka_indel_file),stderr());

loadStrelka()
toc()

tic("HaplotypeCaller")
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
write("x                                             HaplotypeCaller                                                   x",stderr());
write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx",stderr());
# Haplotype caller variant files
haplotypecaller_T_file <- grep(pattern = ".*haplotypecaller.*[TR].*AF.*vep.ann.vcf$",files,value = T)
haplotypecaller_N_file <- grep(pattern = ".*haplotypecaller.*[NB].*AF.*vep.ann.vcf$",files,value = T)
write("Tumour and normal germline calls:",stdout())
printList(list(haplotypecaller_T_file,haplotypecaller_N_file))
loadHaplotypecaller()
toc()

toc()	# measuring the whole stuff
write("\nReady.",stdout())
