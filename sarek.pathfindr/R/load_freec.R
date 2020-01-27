load_freec <-function(result_files) {
  cnvs=NULL
  freec_cnv=NULL
  samplename=NULL
  if (!is.na(result_files["freec_Tratio_file"][1])) {
    # extract sample names
    for (s in 1:length(result_files["freec_Tratio_file"])) 
      samplename[s]=strsplit(basename(result_files["freec_Tratio_file"][s]),'[.]')[[1]][1]
    
    tratio=NULL; tbaf=NULL
    for (s in 1:length(result_files["freec_Tratio_file"])) {
      # tumor log ratio
      cat("Loading tumor log ratio\n")
      temp <- data.table::fread(file = result_files["freec_Tratio_file"][s])
      temp$sample=samplename[s]
      temp=temp[,c(10,1:6)]
      temp$CopyNumber[1:2] <- c(0,4) # to secure range of colors
      tratio=rbind(tratio,temp)
      setkey(tratio,'sample')
      
      # tumor BAF
      cat("Loading tumor BAF\n")
      temp <- data.table::fread(result_files["freec_Tbaf_file"][s])[,1:3]
      temp$sample=samplename[s]
      temp=temp[,c(4,1:3)]
      tbaf=rbind(tbaf,temp)
      setkey(tbaf,'sample')
    }
    # normal log ratio
    cat("Loading normal log ratio\n")
    nratio <- data.table::fread(file = result_files["freec_Nratio_file"])[,1:6]
    # normal BAF
    cat("Loading normal BAF\n")
    nbaf <- data.table::fread(result_files["freec_Nbaf_file"])[,1:3]
    
    # manipulate for plot:
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
    cat("Smooth data for plot\n")
    binned <- NULL
    cat("Writing chromosomes: ")
    for (sample in samplename) for (i in 1:nrow(chrsz)) {
      cat(chrsz$chr[i]," ")
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
    cat("\nCheck what to report\n")
    cnvs=NULL
    for (s in 1:length(samplename)) {
      temp <- as.data.table(read.csv(result_files["freec_cnv_file"][s],sep='\t',header=F))[,1:8]
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
    
    cat("Sorting effects (focal/loss/gain etc.)\n")
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
    cat("Writing results\n")
    freec_cnv=cnvs[order(rank_score,decreasing = T)] # order by rank
    csv_filename=paste0(csv_dir,'/',samplename,'_freec_cnv.csv')
    data.table::fwrite(freec_cnv,file=csv_filename) # write to file
    cat(paste0("Results written to CSV table ",getwd(),"/",csv_filename))
    #tableWrapper(freec_cnv) # table in report
  }
}