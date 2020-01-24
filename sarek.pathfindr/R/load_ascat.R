load_ascat <- function(result_files) {
  ascat_cnv=NULL
  samplename=NULL
  
  if (length(result_files["ascat_Tratio_file"])>0) {
    # extract sample names
    for (s in 1:length(result_files["ascat_Tratio_file"])) samplename[s] = strsplit(basename(result_files["ascat_Tratio_file"][s]),'[.]')[[1]][1]
    
    ascat_tratio=NULL; 
    ascat_tbaf=NULL; 
    ascat_cnv=NULL
    for (s in 1:length(samplename)) {
      # segmented copy number
      cat("Loading segmented copy number\n")
      temp=data.table::fread(file = result_files["ascat_segment_file"][s])
      temp$sample=samplename[s]
      temp=temp[,c(6,1:5)]
      temp$chr=paste0('chr',temp$chr)
      ascat_cnv=rbind(ascat_cnv,temp)
      setkey(ascat_cnv,'sample')
      these_cnvs=temp # for use below
      # tumor log ratio
      cat("Loading tumor log ratio\n")
      temp=fread(file = result_files["ascat_Tratio_file"][s])[,-1]
      colnames(temp)[3]='LogR'
      temp$sample=samplename[s]
      temp=temp[order(temp$Chr,temp$Position),c(4,1:3)]
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
      cat("Loading tumor BAF\n")
      temp=fread(result_files["ascat_Tbaf_file"][s])[,-1]
      colnames(temp)[3]='BAF'
      #browser()
      #temp=temp[which(temp$BAF!=0 & temp$BAF!=1)]
      temp = subset.data.frame(temp,temp$BAF!=0 & temp$BAF !=1)
      temp$sample=samplename[s]
      temp=temp[order(temp$Chr,temp$Position),c(4,1:3)]
      ascat_tbaf=rbind(ascat_tbaf,temp)
      setkey(ascat_tbaf,'sample')
    }
    # normal BAF
    cat("Loading normal BAF\n")
    ascat_nbaf=fread(result_files["ascat_Nbaf_file"])[,-1]
    colnames(ascat_nbaf)[3]='BAF'
    #ascat_nbaf=ascat_nbaf[which(BAF!=0 & BAF!=1)]
    ascat_nbaf = subset.data.frame(ascat_nbaf, ascat_nbaf$BAF!=0 & ascat_nbaf$BAF!=1)
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
    cat("Smooth data for view\n")
    ascat_binned=NULL
    for (sample in samplename) for (i in 1:nrow(chrsz)) {
      temp=data.table(
        sample,
        chr=chrsz$chr[i],
        pos=seq(5e5,chrsz$length[i],5e5),
        cumpos=seq(5e5,chrsz$length[i],5e5)+chrsz$starts[i],
        tratio=NA,
        tmaf=NA)
      
      #browser()
      
      tempLogR = subset.data.frame(ascat_tratio,ascat_tratio$sample == sample & ascat_tratio$Chr == chrsz$chr[i])
      # tempLogR=ascat_tratio[sample][Chr==chrsz$chr[i]]
      tempBAF = subset.data.frame(ascat_tbaf,ascat_tbaf$sample == sample & ascat_tbaf$Chr == chrsz$chr[i])
      #tempBAF=ascat_tbaf[sample][Chr==chrsz$chr[i]]
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
    cat("Adding LOH (if no info from ControlFREEC)\n")
    
    #browser()
    ascat_loh = subset.data.frame(ascat_cnv,ascat_cnv$nMinor == 0 & ascat_cnv$nMajor < 4)
    #ascat_loh=ascat_cnv[nMinor==0 & nMajor<4]
    ascat_loh$chr=paste0('chr',ascat_loh$chr)
    if (!exists('freec_loh')) freec_loh=ascat_loh # to make sure LOH is available if theres no freec
    colnames(ascat_loh)[3:4]=c('start','end')
    #cnvs=ascat_cnv[nMajor+nMinor != 2] # this removes cnLOH
    cnvs=subset.data.frame(ascat_cnv,ascat_cnv$nMajor+ascat_cnv$nMinor != 2)
    cnvs$length=paste0(round((cnvs$endpos-cnvs$startpos)/1e6,2),'M')
    cnvs$rank_score <- 0
    cnvs$rank_terms <- ''
    cnvs$cancer_genes <- ''
    
    cat("Building table view\n")
    for (i in which(cnvs$endpos-cnvs$startpos < 50e6)) {
      #genes=tumorgenes[chr==cnvs$chr[i] & start<cnvs$end[i] & end>cnvs$start[i],`Gene Symbol`] # genes in segment
      tg_chr=subset.data.frame(tumorgenes,tumorgenes$chr == cnvs$chr[i])
      focus_starts=subset.data.frame(tg_chr,tg_chr$start < cnvs$endpos[i])
      focus_ends=subset.data.frame(focus_starts,focus_starts$end > cnvs$start[i])
      genes=focus_ends[,"Gene Symbol"]
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
      if ( cnvs$endpos[i]-cnvs$startpos[i]<3e6 ) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+1
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'focal')
      }
      effect_scale = cnvs$nMajor[i]+cnvs$nMinor[i]
      # if the effect is loss/gain/amp etc
      if (effect_scale == 0) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+2
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'hz_loss')
      } else if (effect_scale > 5) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+2
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'high_amp')
      } else if (effect_scale > 2) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+1
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'dup')
      } else if (effect_scale < 2) {
        cnvs$rank_score[i]=cnvs$rank_score[i]+1
        cnvs$rank_terms[i]=paste(cnvs$rank_terms[i],'del')
      }
    }
    cat("Writing CSV table\n")
    browser()
    ascat_cnv=cnvs[order(cnvs$rank_score,decreasing = T),]
    fwrite(ascat_cnv,file=paste0(csv_dir,'/',samplename,'_ascat_cnv.csv'))
  }
}
