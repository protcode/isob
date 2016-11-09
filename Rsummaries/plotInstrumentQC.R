library(ggplot2)
library(plyr)
library(tools)
library(scales)
require(cowplot)
library(rhdf5)
library(gridExtra)
library(grid)
library(dplyr)
library(reshape2)


mytheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.3)), colhead = list(fg_params=list(cex = 0.3)))
my_minimal_theme<-ttheme_minimal(core = list(fg_params=list(hjust=0, x=0)), rowhead = list(fg_params = list(fontface = 2L, hjust=0, x=0)))



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
  legend <- tmp$grobs[[leg]]
  return(legend)
}
g_legend<-function(a.gplot){  tmp <- ggplot_gtable(ggplot_build(a.gplot));
leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box');
legend <- tmp$grobs[[leg]];
return(legend)};



getInterpolNoise<-function(x, mz){
  h<-approx(x$mz,x$noise, mz );
  ino = data.frame('interp_noise' = h$y);
  return(ino);
}

getNoiseData<-function(noise, spectra, msmsheader){
  n_128<-ddply(noise, .(spec_id), function(noise) getInterpolNoise(noise, 128))
  n_500<-ddply(noise, .(spec_id), function(noise) getInterpolNoise(noise, 500))
  n_800<-ddply(noise, .(spec_id), function(noise) getInterpolNoise(noise, 800))


  noisedata<-merge(n_128,n_500, by='spec_id')
  names(noisedata) <- c('spec_id', 'noise_128', 'noise_500')
  noisedata<-merge(noisedata,n_800, by='spec_id')
  names(noisedata) <- c('spec_id', 'noise_128', 'noise_500', 'noise_800')

  spectra<-merge(spectra, noisedata, by='spec_id');
  msms_spectra<-merge(spectra, msmsheader[c('spec_id', 'scanevent', 'survey_spec')], by='spec_id');
  msms_spectra$scanevent = as.factor(msms_spectra$scanevent -1);
  noise_ms2msms<-merge(msms_spectra, spectra[c('spec_id','noise_500','noise_128', 'type')],by.y='spec_id',
                       by.x='survey_spec')
  noise_ms2msms$noise_ratio = noise_ms2msms$noise_500.x / noise_ms2msms$noise_500.y
  noise_data_summary<-ddply(spectra, .(round(rt), type), summarise, sdv = sd(noise_500),
                            med500=median(noise_500), mn500=mean(noise_500), med128=median(noise_128),mn128=mean(noise_128),
                            med800=median(noise_800),mn800=mean(noise_800))
  noise_data_summary[is.na(noise_data_summary)]<-0
  return(list(noise_ms2msms=noise_ms2msms, noise_data_summary=noise_data_summary))
}


getMaxDis<-function(x){
  h<-density(x)
  maxdist <- max(h$y, na.rm=TRUE)
  return(data.frame(maxdist=maxdist))
  }

justify <- function(x, hjust=center, vjust=center, draw=TRUE){
  w <- sum(x$widths)
  h <- sum(x$heights)
  xj <- switch(hjust,
               center = 0.5,
               left = 0.6*w,
               right=unit(1,'npc') - 0.5*w)
  yj <- switch(vjust,
               center = 0.5,
               bottom = 0.5*h,
               top=unit(1,'npc') - 0.5*h)
  x$vp <- viewport(x=xj, y=yj)
  if(draw) grid.draw(x)
  return(x)
}

getlabel<-function(min, max, step){
  v<-vector();
  for (r in seq(min,max)){
    if(r%%step ==  0){
      v<-c(v,r);
    }
    else{
      v<-c(v,' ');
    }  }
  return(v);
}


getRunParameters<-function(myfile){

  required_data<-c('Creation_date','Instrument Name','MS Run Time (min)')
  mylist = list()
  params<-h5read(file=myfile, name='/rawdata/parameters')
  for (i in required_data){
    v<-params[params$parameter==i,'value']
    if (length(v)==0){
      v<-'Not Set'
    }
    mylist[i]<-v
  }

  aqsum<-tableGrob(t(data.frame(mylist)),  theme=my_minimal_theme, colnames<-c('Acquisition date','Instrument','Analysis time'))
  return(aqsum)
}


getPeptides_cz<-function(myfile){
  allpeptides = data.frame('queryno'=numeric(), 'pepno' = numeric(), 'sequence'=factor(),'mass'=numeric(), 'score'=numeric(), 'spec_id'=factor(), name=factor(), delta=numeric(), peplen=numeric())
  h = try(h5read(myfile, '/imports' ))
  if (class(h) != 'try-error'){
    for (name in h$name){
      if (name !='calibration' & name!=0){
        peptides = try(h5read(myfile, paste('/',name,'/','peptides',  sep='') ))
        queries = try(h5read(myfile, paste('/',name,'/','queries',  sep='') ))
        if (class(peptides) !='try-error'){
          if (nrow(peptides) >0){
            peptides$name = name
            peptides$peplen<-nchar(peptides$sequence)
            peptides<-merge(peptides, queries[,c('queryno', 'spec_id')], by='queryno')
            print(name)
            allpeptides<-rbind(allpeptides, peptides[,c('queryno', 'ishook','pepno', 'sequence', 'mass', 'score', 'name', 'peplen', 'delta', 'spec_id')])
          }
        }
      }
    }
  }else{return(allpeptides)}
  return(allpeptides[allpeptides$pepno==1,])
}

getExpSummary<-function(myfile){
  mytheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.7)), colhead = list(fg_params=list(cex = 0.7)))
  matched_peptides=0
  numhooks=0
  numspectra=0
  uniqueseqs=0
  database='n/a'
  name = 'n/a'
  peptol='n/a'
  itol='n/a'
  ppmerror2pre='n/a'
  #params<-getParameters(myfile)
  h = try(h5read(myfile, '/imports' ))
  if (class(h) != 'try-error'){
    for (name in h$name){
      if (name !='calibration' & name!=0){
        peptides = try(h5read(myfile, paste('/',name,'/','peptides',  sep='') ))
        statistics = try(h5read(myfile, paste('/',name,'/','statistics',  sep='') ))
        matched_peptides = nrow(peptides[peptides$pepno == 1,])
        numhooks<-statistics[statistics$statistic == 'numhooks', c('value')]
        numspectra<-statistics[statistics$statistic == 'numspectra', c('value')]
        ppmerror2pre<-round(statistics[statistics$statistic == 'reldelta_stdev', c('value')],2)
        parameters = h5read(myfile, paste('/',name,'/','parameters',  sep='') )
        itol<-paste(parameters[parameters$parameter=='ITOL',c('value')], parameters[parameters$parameter=='ITOLU',c('value')])
        peptol<-paste(parameters[parameters$parameter=='TOL',c('value')], parameters[parameters$parameter=='TOLU',c('value')])
        database<-paste(basename(parameters[parameters$parameter=='fastafile',c('value')]))
        seq2acc<-h5read(myfile, paste('/',name,'/','seq2acc',  sep='') )
        seq2unique<-ddply(seq2acc[substring(seq2acc$accession, 1,2) !='DD',], .(sequence, name), summarise, m=mean(pepscore), ac_count=length(accession))
        uniqueseqs<-nrow(seq2unique[seq2unique$ac_count==1,])
      }
    }
  }
  if (numspectra == 0){
    #need to look up data from msmsheader table, as no data yet imported
    numspectra  = nrow(h5read(file=myfile, name='/rawdata/msmsheader'))
  }
 summary_frame<-data.frame(datfile=name,uniqueseqs=uniqueseqs,  numspectra=numspectra,matched_peptides=matched_peptides,numhooks=numhooks, ppmerror2pre=ppmerror2pre )
 colnames(summary_frame)<-c('Dat File name', 'Unique\nPeptides', 'MS2', 'Assigned\nMS2','number of hook peptides', 'ppmerror')
 search_para_frame<-data.frame(database=database, peptold=peptol, itol=itol)

summrygrob<-tableGrob(summary_frame, rows=NULL, theme=mytheme)
searchgrob<-tableGrob(search_para_frame, rows=NULL, theme=mytheme,  colnames<-c('Search Database',
                                                                                'Precusor Tolerance', 'Fragment Ion Tolerance'))



return(list(summrygrob, searchgrob))

}

getPeptides<-function(myfile){
  allpeptides = data.frame('query'=numeric(), 'pepno' = numeric(), 'sequence'=factor(),'mass'=numeric(), 'score'=numeric(), 'spec_id'=factor(), name=factor(), da_delta=numeric(), peplen=numeric())
  h = try(h5read(myfile, '/imports' ))
  if (class(h) != 'try-error'){
    for (name in h$name){
      if (name !='calibration' & name!=0){
        peptides = try(h5read(myfile, paste('/',name,'/','peptides',  sep='') ))
        queries = try(h5read(myfile, paste('/',name,'/','queries',  sep='') ))
        if (class(peptides) !='try-error'){
          if (nrow(peptides) >0){
            peptides$name = name
            peptides$peplen<-nchar(peptides$sequence)
            peptides<-merge(peptides, queries[,c('query', 'spec_id')], by='query')
            print(name)
            allpeptides<-rbind(allpeptides, peptides[,c('query', 'is_hook','pepno', 'sequence', 'mass', 'score', 'name', 'peplen', 'da_delta', 'spec_id')])
          }
        }
      }
    }
  }else{return(allpeptides)}
  return(allpeptides[allpeptides$pepno==1,])
}


getPepIDtype<-function(myfile){
  allseq2acc = data.frame('sequence'=factor(),'accession'=factor(), name=factor())
  h = try(h5read(myfile, '/imports' ))
  if (class(h) != 'try-error'){
  for (name in h$name){
    if ((name !='calibration') & (name != '0')){
      seq2acc = try(h5read(myfile, paste('/',name,'/','seq2acc',  sep='') ))
      if (class(seq2acc) !='try-error'){
        if (nrow(seq2acc) >0){
          seq2acc$name = name
          allseq2acc<-rbind(allseq2acc,seq2acc[,c('sequence', 'hittype','name')])
        }}
    }
  }
  #allseq2acc$idtype <- 'FWD'
  #allseq2acc[substring(allseq2acc$accession, 1,2)=='DD',c('idtype')]<-'REV'
  #allseq2acc[grepl('###R\\w{2}###',allseq2acc$accession),c('idtype')]<-'REV'

  # take the FWD hit if there is FWD and reverse:
  seq2acc_tru<-ddply(allseq2acc, c('sequence', 'name'), summarize, hittype=max(hittype))
  return(seq2acc_tru)
  }else {return(data.frame(sequence=c(),name=c(),hittype=c()))}
}

getParameters<-function(myresultfile){
  allparameters = data.frame('section'=factor(),'parameter'=factor(),'value'=factor(), name=factor())
  searchsummary = data.frame(name=factor(), itol=factor(), peptol=factor(), database=factor())
  h = h5read(myresultfile, '/imports' )
  for (name in h$name){
    if ((name !='calibration') & (name != '0')){
      parameters = h5read(myresultfile, paste('/',name,'/','parameters',  sep='') )
      if (1){
        parameters$name = name

        itol<-paste(parameters[parameters$parameter=='ITOL',c('value')], parameters[parameters$parameter=='ITOLU',c('value')])
        peptol<-paste(parameters[parameters$parameter=='TOL',c('value')], parameters[parameters$parameter=='TOLU',c('value')])
        database<-paste(parameters[parameters$parameter=='fastafile',c('value')])
        searchsummary<-rbind(searchsummary,data.frame(name=name,itol=itol, peptol=peptol, database=database))
        allparameters<-rbind(allparameters, parameters)
      }
    }
  }
  return(list(allparameters,searchsummary))
}

getFDRdata<-function(myhdf5file){
  seq2acc_tru<-getPepIDtype(myhdf5file)
  peptides<-getPeptides(myhdf5file)
  if (nrow(peptides) > 0 & nrow(seq2acc_tru)> 0){

    peptides$score<-round(peptides$score)
    peptides<-merge(seq2acc_tru, peptides, by=c('sequence','name'))
    score_data<-ddply(peptides, c('hittype', 'score', 'name'), summarize, total=length(score))
    fwd_score_data<-score_data[score_data$hittype=='FWD',]
    rev_score_data<-score_data[score_data$hittype=='REV',]
    h<-merge(rev_score_data, fwd_score_data, by=c('score', 'name'),all=T)
    h[is.na(h$total.x),c('total.x')] <- 0
    h[is.na(h$hittype.x),c('hittype.x')] <- 'REV'
    h[is.na(h$total.y),c('total.y')] <- 0
    h[is.na(h$hittype.y),c('hittype.y')] <- 'FWD'
    h<-h[rev(order(h$score)),]
    h<-ddply(h,.(name), transform, cumtotal.y = cumsum(total.y))
    h<-ddply(h,.(name), transform, cumtotal.x = cumsum(total.x))
    h$fdr<-h$cumtotal.x/h$cumtotal.y
    h$truhits<-h$cumtotal.y-h$cumtotal.x


  return(h)}else {return(data.frame(score=c(), name=c(), hittype.x=c(), total.x=c(), hittype.y=c(), total.y=c(),
                                   cumtotal.y=c(), cumtotal.x=c(),fdr=c(), truhits=c()))}

}

processMS1summaryData<-function(myfile){
  print('processMS1summaryData')
  ms1summary<-h5read(file=myfile, name='/rawdata/ms1summary')
  msmsheader<-h5read(file=myfile, name='/rawdata/msmsheader')
  spectra = h5read(file=myfile, name='/rawdata/spectra')
  peptides<-getPeptides(myfile)
  j<-NULL
  h<-NULL
  p<-NULL
  i<-NULL
  if (nrow(peptides)>0){
    transmission <- merge(peptides, msmsheader, all.y=T)
    transmission$resolution<-transmission$rt/transmission$fwhm
    transmission[is.infinite(transmission$resolution), 'resolution']  <- 0
    transmission[is.infinite(transmission$theoplates), 'theoplates']  <- 0
    transmission[(!is.na(transmission$is_hook)) & transmission$is_hook==0,'is_hook']<-NA
    assigned_spectra<-nrow(transmission[!is.na(transmission$score) & (transmission$score>-1),])
    triggered<-nrow(transmission)

    hookspectra<-nrow(transmission[!is.na(transmission$is_hook) & (transmission$is_hook==1),])

    m_transmission<-melt(transmission[,c('spec_id', 'rt','is_hook','score')], id='rt', na.rm = T)
    m_transmission$name = 'expt'

    labs<-getlabel(round(min(m_transmission$rt)),round(max(m_transmission$rt)),10)
    j<-ggplot(data=m_transmission, aes(x=rt, ..count.., fill=factor(variable, levels = c('spec_id', 'score','is_hook'), labels=c('MS/MS triggered', 'MS/MS matched', 'MS/MS hook'))))+
     geom_histogram(binwidth = 1, position='identity') +
    guides(fill = guide_legend(title = 'MS/MS success level')) +
    scale_fill_manual(values=c('red', 'blue', 'green')) +
    scale_x_discrete(breaks= labs, labels=labs, name='Retention time (mins)') +
    ylab('frequency')+ggtitle('Frequency of MS/MS Scans per minute') +  theme(legend.position = 'bottom', title = element_text(size=12, face='bold'), axis.title=element_text(size=12, face='plain'),  legend.title=element_text(size=12, face='plain') )

    transmission[(!is.na(transmission$is_hook)) & transmission$is_hook==1, c('spec_id', 'score')]<-NA
    transmission[(!is.na(transmission$score)) & transmission$score>-1, c('spec_id')]<-NA
    theop_transmission<-melt(transmission[transmission$resolution>0, c('spec_id', 'resolution','is_hook','score')], id='resolution', na.rm = T)

    #p<-ggplot(data=theop_transmission, aes(x=resolution, ..count.., col=factor(variable, levels = c('spec_id', 'score','is_hook'), labels=c('MS/MS triggered', 'MS/MS matched', 'MS/MS hook')))) + geom_freqpoly(binwidth=10) + xlim(0,1000)  +
    #guides(col = guide_legend(title = 'MS/MS success level')) + ggtitle('Chromatographic Resolution') + theme(legend.position = 'bottom', title=element_text(size=18,face='bold'), axis.title=element_text(size=18, face='plain') ) +
    #  xlab('Resolution (retention time (seconds) / FWHM (seconds) )')
  }
  atleast1ms2<-ms1summary[ms1summary$numms2>0,]
  if (nrow(atleast1ms2) >0){
    ggplot(data=atleast1ms2, aes(x=factor(x=numms2))) + geom_bar()
  h<-ggplot(data=atleast1ms2, aes(x=factor(x=numms2))) + geom_bar() + xlab('Number MS/MS events') +
  ggtitle('Number of MS/MS Events from\nsingle MS event (all ms2 spectra)') + theme(legend.position = 'bottom', title = element_text(size=8, face='bold'), axis.title=element_text(size=8, face='plain'))
  }
  basepeak_data<-merge(ms1summary, spectra, by='spec_id')
  i<-ggplot(data=basepeak_data, aes(y=basepeak_inten, x=rt)) + geom_line() +
    xlab('Retention Time (mins)') + ylab('base peak')+ ggtitle('Basepeak Chromatogram') +
    theme(legend.position = 'bottom', title = element_text(size=12, face='bold'), axis.title=element_text(size=12, face='plain') )
  H5close()
  return(list(j,i,h))
}

processSpecParams<-function(myfile){
  print('processSpecParams')
  specparams<-h5read(file=myfile,'/rawdata/specparams')

  ms2specids = specparams[specparams$value=='ms2','spec_id']
  k<-specparams[specparams$parameter =='Ion Injection Time (ms)', ]
  injectiontime<-k[k$spec_id %in% ms2specids,'value']


  l<-matrix(quantile(as.numeric(injectiontime), 0.99))[1]
  max_ijt_for_plot<-ifelse(l + 20 > 150, l + 20, 150)

  filltimedata<-data.frame(injectiontime=as.numeric(injectiontime), name='expt')
  sijts<-ddply(filltimedata, .(name), summarise, min_injt=round(min(injectiontime),2),
                       max_injt=round(max(injectiontime)), med_injt=round(median(injectiontime)))
  maxDs<-ddply(filltimedata, .(name), function(filltimedata) getMaxDis(round(filltimedata$injectiontime)))
  maxY<-max(maxDs$maxdist)
  mytheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.4)),
                                               colhead = list(fg_params=list(cex = 0.5)))
  mid_graph = max(filltimedata$injectiontime) / 2.0
  k<-ggplot(data=filltimedata) +
  geom_density(aes(x=(injectiontime), ..density.., fill=name), col=NA, alpha=0.55) +
  xlab('fill time (ms, bin width=1)') + theme(legend.position = 'None', title=element_text(size=18,face='bold'), axis.title=element_text(size=18, face='plain'))  +
  ggtitle('MS2 fill time Distribution') +
  xlim(0,max_ijt_for_plot) + annotation_custom(tableGrob(sijts[2:4] , rows=NULL, colnames<-c('min\ninj.time',
                                                               'max\ninj.time', 'med\ninj.time'), theme=mytheme),
                                               xmin=mid_graph - 8, xmax=mid_graph + 8, ymin= maxY / 2, ymax=1 )
  H5close()
  return(k)
}


processNoiseData<-function(myfile){
  print('processNoiseData')
  spectra = h5read(file=myfile, name='/rawdata/spectra')
  msmsheader<-h5read(file=myfile, name='/rawdata/msmsheader')
  noise<-h5read(file=myfile, name='/rawdata/noise')
  noisedata<-getNoiseData(noise, spectra, msmsheader)
  noise_summary<-noisedata$noise_data_summary
  v1<-ggplot(data=noise_summary, aes(`round(rt)`, med500))+geom_line() +
               facet_grid(`type`~., scale='free') + ggtitle('Noise Cut off at m/z 500' ) +
  ylab('Median noise (Intensity) at m/z 500') +
  xlab('Retention time (mins)')+  theme(legend.position = 'bottom', title = element_text(size=12, face='bold'), axis.title=element_text(size=12, face='plain') )

  v2<-ggplot(data=noise_summary, aes(`round(rt)`, med128))+geom_line() +
    facet_grid(`type`~., scale='free') + ggtitle('Noise Cut off at m/z 128' ) +
    ylab('Median noise (Intensity) at m/z 128') +
    xlab('Retention time (mins)')+  theme(legend.position = 'bottom', title = element_text(size=12, face='bold'), axis.title=element_text(size=12, face='plain') )

  v3<-ggplot(data=noise_summary, aes(`round(rt)`, med800))+geom_line() +
    facet_grid(`type`~., scale='free') + ggtitle('Noise Cut off at m/z 800' ) +
    ylab('Median noise (Intensity) at m/z 800') +
    xlab('Retention time (mins)')+  theme(legend.position = 'bottom', title = element_text(size=12, face='bold'), axis.title=element_text(size=12, face='plain') )
 # noise_ms2msms<-noisedata$noise_ms2msms
  #w<-ggplot(data=ddply(noise_ms2msms, .(scanevent), summarise, me =mean(noise_ratio)), aes(x=scanevent, y=me))+ geom_bar(stat='identity', position='dodge') +
  #theme(legend.position = 'bottom', title = element_text(size=12, face='bold'), axis.title=element_text(size=12, face='plain')) + ylab('Mean ratio MS : MS/MS noise cut off') + xlab('MS/MS Event') + ggtitle('Plot of mean ratio of MS Noise to MS/MS spectral noise (at m/z 500)')
  H5close()
  return(list(v1,v2,v3))

}

getCycletimeData<-function(myfile){
  L_spec = h5read(myfile,"rawdata/spectra")

  # filter to contain only MS1

  L_spec.1 = L_spec[L_spec$type == "ms", ]



  # calculate cycle time

  for(i in 1 : length (L_spec.1$type)-1){

    L_spec.1$cycletime[i] = (L_spec.1$rt[i+1]-L_spec.1$rt[i])*60

  }

  maxcycletime = round(max(L_spec.1$cycletime))
  mincycletime = round(min(L_spec.1$cycletime))
  return(ggplot(data=L_spec.1[L_spec.1$rt>20 & L_spec.1$rt<90,]) + geom_density(aes(x=cycletime))+xlim(mincycletime-0.5, maxcycletime+0.5 )+xlab('cycle time (s)')+
           ggtitle('Distribution of cycle times\nbetween consecutive MS1 Scans') +
           theme(legend.position = 'None', title=element_text(size=18,face='bold'), axis.title=element_text(size=18, face='plain')))
}



processMs1Ms2Data<-function(myfile){
  print('processMs1Ms2Data')
  # all to do with ms1 / ms2 raw data
  msmsheader<-h5read(file=myfile, name='/rawdata/msmsheader')
  ms2summary<-h5read(file=myfile, name='/rawdata/ms2summary')
  msmsheader$name <-'expt'
  msmsheader$fwhm_seconds = msmsheader$fwhm * 60
  msmsheader$log_precursor<-log10(msmsheader$precinten)
  ms2summary$log_tic_mz<-log10(ms2summary$tic)
  ms2summary$name <-'expt'

  msmsheader$p2t<-msmsheader$c12 / msmsheader$thresh
  msmsheader$rtdelta = (msmsheader$rt - msmsheader$rtapex) * 60

  a<-ggplot(msmsheader, aes(log_precursor, ..count..)) + geom_freqpoly(binwidth = 0.2) +
               xlab('Log10 Precursor Intensity (bin width=0.2)') + xlab('Log10 Precursor') +
  ggtitle('Log10 Precursor\nIntensity Distribution') + theme(legend.position = 'none',
                                                              title=element_text(size=8,face='bold'))
  precursor_summaries<-ddply(msmsheader, .(name), summarise, avg=round(mean(msmsheader[!is.infinite(msmsheader$log_precursor),'log_precursor'],  na.rm = TRUE),1),
               median=round(median(log_precursor),1), quant5=round(quantile(log_precursor,0.05),1),
                             quant95=round(quantile(log_precursor, 0.95),1))
  tic_summaries<-ddply(ms2summary, .(name), summarise, avg=round(mean(log_tic_mz),1), median= round(median(log_tic_mz),1), quant5= round(quantile(log_tic_mz,0.05),1),
                       quant95=round(quantile(log_tic_mz, 0.95),1))
  b<-ggplot(ms2summary, aes(log_tic_mz, ..count..)) + geom_freqpoly(binwidth = 0.2) +
               xlab('Log10 TIC (MS2) Intensity (bin width=0.2)') +   ggtitle('Log10 TIC (ms2)\nIntensity Distribution') +
  xlab('Log10 TIC')+theme(legend.position = 'none', title=element_text(size=8,face='bold'))

  precursor_data<-merge(ms2summary, msmsheader[!is.infinite(msmsheader$log_precursor),], by='spec_id')
  #write.table(precursor_data, file='pcd.txt', sep='\t', row.names = F)
  precursor_data$name = 'expt'
  coefs <- ddply(precursor_data, .(name), function(df) {m <- lm(log_tic_mz ~ log_precursor, data=df);
  data.frame(precursor_intercept = coef(m)[1], precursor_slope = coef(m)[2]); })


  c<-ggplot(data=precursor_data, aes(x=log_precursor, y=log_tic_mz)) + geom_point(alpha=0.2, size=1) +
               geom_abline(data=coefs, aes(intercept=precursor_intercept, slope=precursor_slope)) +
  xlab('Log10 Precursor Intensity') + ylab('Log10 TIC (ms2) Intensity')  + xlim(4,9) + ylim(4,9) +
  ggtitle('Log10 Precursor vs\nLog10 TIC Intensities') + theme(legend.position = 'bottom',
                                                                  title=element_text(size=8,face='bold'))   +
  guides(col= guide_legend(title = 'Source', override.aes = list(alpha=1, size=6)))


  d<-ggplot(data=precursor_data, aes(y=log_precursor, x=precmz))  +
  geom_point( size=1, alpha = 0.2) + geom_smooth(color='green',                                            size=0.5)  +
  ylab('Log10 Precursor Intensity') + xlab('Precursor m/z') + ylab('Log10 Precursor') +
  ggtitle('Log10 Precursor Intensity vs m/z')

  e<-ggplot(data=precursor_data, aes(y=log_tic_mz, x=precmz))  +
               geom_point(alpha = 0.2, size=1, show.legend=FALSE)  + geom_smooth(show.legend=FALSE, size=0.5,
                                                                                 col='green') +
  ylab('Log10 TIC (ms2) Intensity') + xlab('Precursor m/z') +
  ggtitle('Log10 TIC (ms2) Intensity vs m/z') + theme(legend.position = 'none')

  precursor_summaries$name<-'Precursor'
  tic_summaries$name<-'TIC (ms2)'
  tableSummaries<-rbind(precursor_summaries, tic_summaries)


  mytheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.6)),
         colhead = list(fg_params=list(cex = 0.6)))
  tabA<-tableGrob(tableSummaries, rows=NULL, theme=mytheme, colnames<-c('name', 'average', 'median',  '5%\nQuantile','95%\nQuantile'))

  coefs[2:3]<-round(coefs[2:3], 1)
  tabB<-justify(tableGrob(coefs[2:3], rows=NULL, theme=mytheme, colnames<-c( 'intercept', 'slope')), 'left', 'top', FALSE)

  tab1<-arrangeGrob(tabA, tabB)


  fwhm_spread_data <- ddply(msmsheader[msmsheader$fwhm_seconds >0,], .(name), summarise, median = median(round(fwhm_seconds,1)), quantile05=as.numeric(quantile(round(fwhm_seconds,1), 0.05)), quantile95=as.numeric(quantile(round(fwhm_seconds,1), 0.95)))

  mytheme <- gridExtra::ttheme_default(core = list(fg_params=list(cex = 0.5)),
                                       colhead = list(fg_params=list(cex = 0.4)))


  tg<-tableGrob(fwhm_spread_data[2:4], rows=NULL, colnames<-c('median', '5%\nquantile', '95%\nquantile'), theme=mytheme)
  #tg<-tableGrob(fwhm_spread_data, rows=NULL,  theme=mytheme)
  #print(tg)



  maxcount<-max(hist(msmsheader$fwhm_seconds,breaks=seq(min(msmsheader$fwhm_seconds), max(msmsheader$fwhm_seconds)+1), plot=FALSE)$counts)
  #s<-ggplot(data=msmsheader) + geom_freqpoly(aes(x=fwhm_seconds, ..count..), binwidth=1) +
  #ggplot(data=msmsheader) + geom_freqpoly(aes(x=fwhm_seconds, ..sum..), binwidth=1)
  s<-ggplot(data=msmsheader[msmsheader$fwhm_seconds >0,]) + geom_freqpoly(aes(fwhm_seconds), binwidth=1) +
         xlim(0,60) + ggtitle('Distribution of FWHM of\nChromatographic Peaks') +
          geom_vline(data=fwhm_spread_data, aes(xintercept=as.numeric(quantile95), col=name), linetype=2) +
          geom_vline(data=fwhm_spread_data, aes(xintercept=as.numeric(quantile05), col=name), linetype=2) +
          theme(legend.position = '', title=element_text(size=18, face='bold'), axis.title=element_text(size=18, face='plain') ) + xlab('FWHM (seconds, bin width = 1s)') +
     annotation_custom(tg, xmin = 25, ymin = 4000, ymax=5000, xmax=60)
  #print(s)

  maxcount<-max(hist(msmsheader$rtdelta,breaks=seq(min(msmsheader$rtdelta), max(msmsheader$rtdelta)+1), plot=FALSE)$counts)

  apex_spread_data<- ddply(msmsheader, .(name), summarise, median = round(median(rtdelta),1),
                           quantile05=round(quantile(rtdelta, 0.05),1), quantile95=round(quantile(rtdelta, 0.95),1))
  if (apex_spread_data$quantile95 > 30) {
    apex_spread_data$quantile95_p <- 31
  } else { apex_spread_data$quantile95_p <- apex_spread_data$quantile95 }
  if (apex_spread_data$quantile05 <  -30) {
    apex_spread_data$quantile05_p =  -31
  } else { apex_spread_data$quantile05_p <- apex_spread_data$quantile05 }



  t<-ggplot(data=msmsheader) + geom_freqpoly(aes(x=rtdelta, ..count..),binwidth=1)+
         xlim(-31,31)+xlab('Delta to Apex (seconds, bin width = 1s)')+ geom_vline(data=apex_spread_data, aes(xintercept=quantile95_p),linetype=2)+ geom_vline(data=apex_spread_data, aes(xintercept=quantile05_p),linetype=2)+
         ggtitle('Distribution of distance MS/MS to Apex') +
         theme(legend.position = 'None', title=element_text(size=18, face='bold'), axis.title=element_text(size=18, face='plain') ) +
         annotation_custom(tableGrob(apex_spread_data[2:4], rows=NULL, colnames<-c( 'median', '5%\nquantile', '95%\nquantile'), theme=mytheme) ,
                          xmin = 5,  xmax=20, ymin=maxcount * 0.3, ymax=maxcount * 0.6 )



  g<-ggplot(msmsheader, aes(factor(charge), ..count..))  + geom_bar( stat='count',
         position='dodge') + xlab('Charge State') + ggtitle('Distribution of Precursor Charge States\n(for all ms2 spectra)') +
         theme(legend.position = 'bottom',title=element_text(size=8, face='bold') )
  p2t_meds <- ddply(msmsheader, .(name), summarise, n = length(p2t), med = median(p2t), avg = mean(p2t) )
  n<-ggplot(data=msmsheader, aes(y=p2t, x=name)) + geom_violin( show.legend=FALSE)  +
  theme(legend.position = 'bottom', title=element_text(size=18, face='bold'), axis.text=element_text(size=18, face='plain'), axis.text.x=element_blank(), axis.ticks.x = element_blank() ) +
  geom_boxplot(width=.1, show.legend=FALSE) + ggtitle('Distribution of P2T values') +
  xlab('') + ylab('P2T value') + ylim(-10,100) + geom_text(data = p2t_meds, aes(x = name, y = -10, label = paste0('n=',n))) + geom_text(data = p2t_meds, aes( x=name, y = -5, label = paste0('median=',round(med,1))))
  s2i_meds <- ddply(msmsheader, .(name), summarise, n = length(s2i), med = median(s2i),
               avg = mean(s2i) )
  o<-ggplot(data=msmsheader, aes(y=s2i, x=name)) + geom_violin(show.legend = FALSE) +
  geom_boxplot(width=.1,show.legend = FALSE)+ggtitle('Distribution of S2I values') +
  theme(legend.position = 'bottom',title=element_text(size=18, face='bold'), axis.text = element_text(size=18, face='plain'), axis.text.x=element_blank(), axis.ticks.x = element_blank() )  +
  xlab('')+ ylab('S2I value') + geom_text(data = s2i_meds, aes(x = name, y = -0.07, label = paste0('n=',n))) + geom_text(data = s2i_meds, aes(x = name, y = -0.025, label = paste0('median=',round(med,1))))
  H5close()
  return(list(a,b,c,d,e,g,n,o,s,t,tab1))
}


processDatfileData<-function(myfile){
  print('processDatfileData')
  msmsheader<-h5read(file=myfile, name='/rawdata/msmsheader')
  f<-NULL
  m1<-NULL
  m2<-NULL
  fdr1<-NULL
  fdr2<-NULL
  peptides<-getPeptides(myfile)

  seq2acc_tru<-getPepIDtype(myfile)
  if (nrow(peptides)>0){

    peptides$score<-round(peptides$score)
    peptides<-merge(seq2acc_tru, peptides, by=c('sequence','name'))
    dodge <- position_dodge(width = 0.8)

    peptides$ppmdelta<-(peptides$da_delta / peptides$mass ) * 1E6
    peptides<-peptides[peptides$name != FALSE,]
    mz2delta<-merge(msmsheader[,c('precmz','spec_id')], peptides, by='spec_id')
    f<-ggplot(data=mz2delta, aes(x=precmz, y=ppmdelta))+ geom_point(alpha=0.2, size=1) +
                 geom_smooth(show.legend=FALSE, size=0.5, colour='green', linetype=1) +
    ylab('ppm delta')+ xlab('Precursor m/z')+ggtitle('Precursor PPM Delta vs m/z') +
    theme(legend.position = 'none')
    score_meds <- ddply(peptides, .(name), summarise, n = length(score), med = median(score), avg = mean(score) )
    m1<-ggplot(peptides, aes(y=score, x = name)) + geom_violin() + geom_boxplot(width=.5) +
      geom_text(size=4,data = score_meds, aes(x = name, y = -15, label = paste0('n=', n), size = 18)) +
      geom_text(size=4,data = score_meds, aes(x = name, y = -5, label = paste0('median=', round(med,1)), size = 10)) +
      theme(legend.position='none', legend.position = 'bottom') + theme(title=element_text(size=18,face='bold'), axis.text=element_text(size=18, face='plain'), legend.title = element_text(size=18, face='plain') ) +
    ggtitle('Distribution of SSM Scores') + xlab('') + ylab('Mascot Ion Score') +
    coord_cartesian(ylim = c(-20,120))
    score_bytype_meds <- ddply(peptides, .(name, hittype), summarise, n = length(score), med = median(score), avg = mean(score) )
    m2<-ggplot(data=peptides, aes(y=score, x=name, fill=hittype)) + geom_violin(position = dodge ) +
      geom_boxplot(width=.1, position = dodge, show.legend=FALSE) + xlab('') +
      guides(fill= guide_legend(title = 'Hit Type')) + theme(legend.position = 'bottom',
                                                             title=element_text(size=18,face='bold'), axis.text=element_text(size=18, face='plain') ) +
    ggtitle('Distribution of SSM scores') + coord_cartesian(ylim = c(-20,120)) +
    geom_text(size=4, data = score_bytype_meds, aes(y = -5, label = paste0('median=',round(med,1)),
                                        ymax = 120),position = dodge) + geom_text(size=4,data = score_bytype_meds, aes( y = -15,  label = paste0('n=',n), ymax = 120),position = dodge)

    fdrd<-getFDRdata(myfile)
#    pp<-getParameters(myfile)
#    searchsummary = pp[[2]]
#    searchsummary$summary<-paste(searchsummary$peptol,searchsummary$itol, basename(as.character(searchsummary$database)))
#    fdrd<-merge(fdrd, searchsummary[,c('name', 'summary','itol','peptol')], by='name')


    # this will not currently work for more than 1 dat file (name)
    print(nrow(fdrd))
    if (sum(fdrd$fdr)>0) {
    d<-approx(fdrd$fdr,fdrd$score, 0.01);
    interpolfdr<-data.frame(y=d$x, x=round(d$y))

    # this will not currently work for more than 1 dat file (name)
    d<-approx(fdrd$fdr,fdrd$truhits, 0.01);
    interpoltrucount<-data.frame(x=d$x, y=round(d$y))

    fdr1<-ggplot(data=fdrd)+geom_line(aes(y=truhits, x=fdr, col=name))+ geom_vline(aes(xintercept=0.01), linetype=2) + xlab('FDR') + ylab('# true positives (true spectra minus decoy)') +
      geom_text(data=interpoltrucount, aes(x=x,y=y,label=y), nudge_y = 1,nudge_x = 0.1) + ggtitle('True spectra against FDR') + theme(title=element_text(size=18,face='bold'), axis.title=element_text(size=18, face='plain'), legend.title=element_text(size=18, face='plain') )  +
      geom_segment(data=interpoltrucount, aes(x = x+0.08, y = y+1, xend = x, yend = y), arrow = arrow(length = unit(0.3, "cm")))
    fdr2<-ggplot(data=subset(fdrd, score<51))+geom_line(aes(y=fdr, x=score, col=name)) + geom_hline(aes(yintercept=0.01), linetype=2) + ylab('FDR') + xlab('Mascot Ion Score (upper limit 50)')+
      geom_text(data=interpolfdr, aes(x=x,y=y,label=x), nudge_x = 1, nudge_y=0.025) + ggtitle('FDR against Mascot Score') + theme(title=element_text(size=18,face='bold'), axis.title=element_text(size=18, face='plain') ) +
       geom_segment(data=interpolfdr, aes(x = x+1, y = y+0.022, xend = x, yend = y), arrow = arrow(length = unit(0.3, "cm")))
    }
  }

  H5close()
  return(list(f, m1, m2, fdr1, fdr2))
}

runall<-function(myfile, outfile){
  print(paste('running for', myfile))
  pdf(outfile)
  aqsum<-getRunParameters(myfile)
  sumaries<-getExpSummary(myfile)
  allplots<-arrangeGrob(grobs=list(aqsum,sumaries[[1]], sumaries[[2]]), ncol=1, nrow=3)
  grid.arrange(allplots)
  out<-processMs1Ms2Data(myfile)

  a<-out[[1]]
  b<-out[[2]]
  c<-out[[3]]
  d<-out[[4]]
  e<-out[[5]]
  g<-out[[6]]
  n<-out[[7]]
  o<-out[[8]]
  s<-out[[9]]
  t<-out[[10]]
  tab1<-out[[11]]
  out<-processDatfileData(myfile)
  f<-out[[1]]
  m1<-out[[2]]
  m2<-out[[3]]
  fdr1<-out[[4]]
  fdr2<-out[[5]]
  out<-processMS1summaryData(myfile)
  j<-out[[1]]
  #p<-out[[2]]
  i<-out[[2]]
  h<-out[[3]]
  out<-processNoiseData(myfile)
  v1<-out[[1]]
  v2<-out[[2]]
  v3<-out[[3]]
  #w<-out[[4]]
  grid.arrange(grobs=list(a, b, c, tab1), nrow = 2, ncol=2)
  if (!is.null(f)){
    grid.arrange(grobs=list(d + theme(legend.position = 'none'), e, f), nrow = 2, ncol=2)
    if (!is.null(h)){
    grid.arrange(grobs=list(g, h),  nrow=2)
    }else{print(g)}
    grid.arrange(grobs=list(i, j),  nrow=2)
  }else{
  grid.arrange(grobs=list(d + theme(legend.position = 'none'), e), nrow = 1, ncol=2)
  if (!is.null(h)){
    grid.arrange(grobs=list(g, h),  nrow=2)
  }else{print(g)}
  }

  k<-processSpecParams(myfile)
  ct<-getCycletimeData(myfile)
  print(ct)
  print(k)
  print(m1)
  print(m2)
  print(fdr1)
  print(fdr2)
  print(n)
  print(o)
  #print(p)
  print(s)
  print(t)

  grid.arrange(grobs=list(v2, v1  + theme(legend.position = 'none'), v3), ncol=2, nrow=2)

  dev.off()
}
