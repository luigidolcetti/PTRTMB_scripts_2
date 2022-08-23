areas<- read.csv("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/bigTables/areaTable.csv",
                 )

bloodyTable <- read.csv("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/bigTables/IMCSamples_data_PFS_BOR.csv")
bloodyTable <- bloodyTable[-1,]

colnames(bloodyTable)[1]<-'sample'

CD8_1and2 <- areas[areas$super.cluster=='CD8' & areas$sub.cluster %in% c(1,2),]

CD8_1and2 <- aggregate(sub.area~sample+uid+sample.area+super.area+pfs,CD8_1and2,sum)
PercParental <- CD8_1and2$sub.area/CD8_1and2$super.area*100
PercParental[is.na(PercParental)]<-0
PercParental <- cbind.data.frame(CD8_1and2,percParental=PercParental)
# percTotal <- CD8_1and2$sub.area/CD8_1and2$sample.area*100
# percTotal[is.na(percTotal)]<-0
# percParentalagg <- aggregate(x=list(perc=PercParental),by=list(sample=CD8_1and2$sample,pfs=CD8_1and2$pfs),FUN=mean)

# cor.test(PercParental, CD8_1and2$pfs)
# cor.test(percTotal, CD8_1and2$pfs)

bloodyMix <- dplyr::inner_join(bloodyTable,PercParental,by='sample')
# plot(bloodyMix$CD33.Monocytes.c1,bloodyMix$percParental)
# cor.test(bloodyMix$CD33.Monocytes.c1, bloodyMix$percParental)

mn <- aggregate(percParental~pfs+sample, bloodyMix,mean)
sdmn <- aggregate(percParental~pfs+sample, bloodyMix,sd)
sts <- cbind.data.frame(pfs = mn$pfs,
                        mean = mn$percParental,
                        meanP = mn$percParental+sdmn$percParental,
                        meanM = mn$percParental-sdmn$percParental)

# bloodyMix <- dplyr::inner_join(bloodyMix,mn,by='sample')
# bloodyMix <- dplyr::inner_join(bloodyMix,sdmn,by='sample')


jPalette<-setNames( RColorBrewer::brewer.pal(9,'Paired'),sort(unique(bloodyMix$pfs)))

  postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA",
                              "CD8_1and2_vs_pfs.eps"),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  par(xpd=TRUE)
  plot(data.frame(jitter(bloodyMix$pfs,3),
                  bloodyMix$percParental),
       bg=jPalette[as.character(bloodyMix$pfs)],
       col='gray40',
       pch=21,
       ylab='% of parental area',
       xlab='pfs',
       
       bty='l',
       las=1,
       ylim=c(min(c(sts$meanM,bloodyMix$percParental)),max(sts$meanP,bloodyMix$percParental)))
  # mtext(side = 3,adj=1,prettyLabels[i])
  for (ii in 1:nrow(sts)){
    lines(matrix(c(sts$pfs[ii]-0.5,sts$pfs[ii]+0.5,sts$mean[ii],sts$mean[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(sts$pfs[ii],sts$pfs[ii],sts$meanP[ii],sts$meanM[ii]),
                 ncol = 2,
                 byrow = F),
          col=jPalette[as.character(sts$pfs[ii])],
          lwd=2)
    lines(matrix(c(sts$pfs[ii]-0.2,sts$pfs[ii]+0.2,sts$meanP[ii],sts$meanP[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(sts$pfs[ii]-0.2,sts$pfs[ii]+0.2,sts$meanM[ii],sts$meanM[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
  }
  
  
  cr <- cor.test(bloodyMix$pfs,bloodyMix$percParental)
  lr <- lm(percParental~pfs,bloodyMix)
  xrng<-data.frame(pfs=par('usr')[c(1,2)])
  yrng<-predict(lr,xrng)
  
  lines(x=unlist(xrng),y=yrng,col='red',lty=3,lwd=3)
  legend('topleft',
         legend = c(paste0('Correlation = ',
                           round(cr$estimate,3)),
                    paste0('Correlation p val = ',
                           formatC(cr$p.value,format = 'e', digits = 2),
                           ' (',gtools::stars.pval(cr$p.value),')')),
                  
         bty = 'n',
         cex=0.8,
         inset=c(-0.05,-0.1))
  dev.off()

##########################
##########################

  mn <- aggregate(percParental~pfs+sample+CD33.Monocytes.c1, bloodyMix,mean)
  sdmn <- aggregate(percParental~pfs+sample+CD33.Monocytes.c1, bloodyMix,sd)
  sts <- cbind.data.frame(pfs = mn$pfs,
                          cd33 = mn$CD33.Monocytes.c1,
                          mean = mn$percParental,
                          meanP = mn$percParental+sdmn$percParental,
                          meanM = mn$percParental-sdmn$percParental)
  
  postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA",
                              "CD8_1and2_vs_CD33MONO.eps"),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  par(xpd=TRUE)
  plot(data.frame(jitter(bloodyMix$CD33.Monocytes.c1,3),
                  bloodyMix$percParental),
       bg=jPalette[as.character(bloodyMix$pfs)],
       col='gray40',
       pch=21,
       ylab='% of parental area',
       xlab='% of CD33+ monocytes ',
       
       bty='l',
       las=1,
       ylim=c(min(c(sts$meanM,bloodyMix$percParental)),max(sts$meanP,bloodyMix$percParental)))
  # mtext(side = 3,adj=1,prettyLabels[i])
  for (ii in 1:nrow(sts)){
    lines(matrix(c(sts$cd33[ii]-0.5,sts$cd33[ii]+0.5,sts$mean[ii],sts$mean[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(sts$cd33[ii],sts$cd33[ii],sts$meanP[ii],sts$meanM[ii]),
                 ncol = 2,
                 byrow = F),
          col=jPalette[as.character(sts$pfs[ii])],
          lwd=2)
    lines(matrix(c(sts$cd33[ii]-0.2,sts$cd33[ii]+0.2,sts$meanP[ii],sts$meanP[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(sts$cd33[ii]-0.2,sts$cd33[ii]+0.2,sts$meanM[ii],sts$meanM[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
  }
  
  
  cr <- cor.test(bloodyMix$CD33.Monocytes.c1,bloodyMix$percParental)
  lr <- lm(percParental~CD33.Monocytes.c1,bloodyMix)
  xrng<-data.frame(CD33.Monocytes.c1=par('usr')[c(1,2)])
  yrng<-predict(lr,xrng)
  
  lines(x=unlist(xrng),y=yrng,col='red',lty=3,lwd=3)
  legend('topleft',
         legend = c(paste0('Correlation = ',
                           round(cr$estimate,3)),
                    paste0('Correlation p val = ',
                           formatC(cr$p.value,format = 'e', digits = 2),
                           ' (',gtools::stars.pval(cr$p.value),')')),
         
         bty = 'n',
         cex=0.8,
         inset=c(-0.05,-0.1))
  dev.off()
  
  
  
  
  ##########################
  ##########################
  
  mn <- aggregate(percParental~pfs+sample+CD8.Cent.Mem.c1, bloodyMix,mean)
  sdmn <- aggregate(percParental~pfs+sample+CD8.Cent.Mem.c1, bloodyMix,sd)
  sts <- cbind.data.frame(pfs = mn$pfs,
                          cd33 = mn$CD8.Cent.Mem.c1,
                          mean = mn$percParental,
                          meanP = mn$percParental+sdmn$percParental,
                          meanM = mn$percParental-sdmn$percParental)
  
  postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA",
                              "CD8_1and2_vs_CD8CM.eps"),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  par(xpd=TRUE)
  plot(data.frame(jitter(bloodyMix$CD8.Cent.Mem.c1,3),
                  bloodyMix$percParental),
       bg=jPalette[as.character(bloodyMix$pfs)],
       col='gray40',
       pch=21,
       ylab='% of parental area',
       xlab='% of CD8 central memory',
       
       bty='l',
       las=1,
       ylim=c(min(c(sts$meanM,bloodyMix$percParental)),max(sts$meanP,bloodyMix$percParental)))
  # mtext(side = 3,adj=1,prettyLabels[i])
  for (ii in 1:nrow(sts)){
    lines(matrix(c(sts$cd33[ii]-0.5,sts$cd33[ii]+0.5,sts$mean[ii],sts$mean[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(sts$cd33[ii],sts$cd33[ii],sts$meanP[ii],sts$meanM[ii]),
                 ncol = 2,
                 byrow = F),
          col=jPalette[as.character(sts$pfs[ii])],
          lwd=2)
    lines(matrix(c(sts$cd33[ii]-0.2,sts$cd33[ii]+0.2,sts$meanP[ii],sts$meanP[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(sts$cd33[ii]-0.2,sts$cd33[ii]+0.2,sts$meanM[ii],sts$meanM[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
  }
  
  
  cr <- cor.test(bloodyMix$CD8.Cent.Mem.c1,bloodyMix$percParental)
  lr <- lm(percParental~CD8.Cent.Mem.c1,bloodyMix)
  xrng<-data.frame(CD8.Cent.Mem.c1=par('usr')[c(1,2)])
  yrng<-predict(lr,xrng)
  
  lines(x=unlist(xrng),y=yrng,col='red',lty=3,lwd=3)
  legend('topleft',
         legend = c(paste0('Correlation = ',
                           round(cr$estimate,3)),
                    paste0('Correlation p val = ',
                           formatC(cr$p.value,format = 'e', digits = 2),
                           ' (',gtools::stars.pval(cr$p.value),')')),
         
         bty = 'n',
         cex=0.8,
         inset=c(-0.05,-0.1))
  dev.off()
  
  
  #############################
  #############################
  #############################
  
  
  CD8_1and2 <- areas[areas$super.cluster=='CD8' & areas$sub.cluster %in% c(1,2),]
  
  CD8_1and2 <- aggregate(sub.area~sample+uid+sample.area+super.area+pfs,CD8_1and2,sum)
  PercParental_CD8 <- CD8_1and2$sub.area/CD8_1and2$super.area*100
  PercParental_CD8[is.na(PercParental_CD8)]<-0
  PercParental_CD8 <- cbind.data.frame(CD8_1and2,percParental_CD8=PercParental_CD8)
  
  UNKNW_8 <- areas[areas$super.cluster=='UKNW' & areas$sub.cluster==8,]
  UNKNW_8 <- aggregate(sub.area~sample+uid+sample.area+super.area+pfs,UNKNW_8,sum)
  PercParental_UNKNW <- UNKNW_8$sub.area/UNKNW_8$super.area*100
  PercParental_UNKNW[is.na(PercParental_UNKNW)]<-0
  PercParental_UNKNW <- cbind.data.frame(UNKNW_8,percParental_U8=PercParental_UNKNW)
  
  postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA",
                              "CD8_1and2_vs_UKNWN_8.eps"),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  par(xpd=TRUE)
  plot(PercParental_CD8$percParental_CD8,PercParental_UNKNW$percParental_U8,
       bg=jPalette[as.character(PercParental_CD8$pfs)],
       col='gray40',
       pch=21,
       ylab='CD8+ (1&2) % of parental area',
       xlab='UNKN (8) % of parental area',
       
       bty='l',
       las=1)

  cr <- cor.test(PercParental_CD8$percParental_CD8,PercParental_UNKNW$percParental_U8)
  test <- cbind.data.frame(y=PercParental_UNKNW$percParental_U8,x=PercParental_CD8$percParental_CD8)
  lr <- lm(y~x,test)
  xrng<-data.frame(x=par('usr')[c(1,2)])
  yrng<-predict(lr,xrng)
  
  lines(x=unlist(xrng),y=yrng,col='red',lty=3,lwd=3)
  legend('topleft',
         legend = c(paste0('Correlation = ',
                           round(cr$estimate,3)),
                    paste0('Correlation p val = ',
                           formatC(cr$p.value,format = 'e', digits = 2),
                           ' (',gtools::stars.pval(cr$p.value),')')),
         
         bty = 'n',
         cex=0.8,
         inset=c(-0.05,-0.1))
  dev.off()
  
  
  #############################
  #############################
  #############################
  jPalette<-setNames( RColorBrewer::brewer.pal(9,'Paired'),bloodyTable$sample[order(unique(bloodyTable$t.PFS))])
  
  CD8_1and2 <- areas[areas$super.cluster=='CD8' & areas$sub.cluster %in% c(1,2),]
  
  CD8_1and2 <- aggregate(sub.area~sample+uid+sample.area+super.area+pfs,CD8_1and2,sum)
  PercParental_CD8 <- CD8_1and2$sub.area/CD8_1and2$sample.area*100
  PercParental_CD8[is.na(PercParental_CD8)]<-0
  PercParental_CD8 <- cbind.data.frame(CD8_1and2,percParental_CD8=PercParental_CD8)
  
  UNKNW_8 <- areas[areas$super.cluster=='UKNW' & areas$sub.cluster==8,]
  UNKNW_8 <- aggregate(sub.area~sample+uid+sample.area+super.area+pfs,UNKNW_8,sum)
  PercParental_UNKNW <- UNKNW_8$sub.area/UNKNW_8$sample.area*100
  PercParental_UNKNW[is.na(PercParental_UNKNW)]<-0
  PercParental_UNKNW <- cbind.data.frame(UNKNW_8,percParental_U8=PercParental_UNKNW)
  
  
  postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA",
                              "CD8_1and2 TOTAL_vs_UKNWN_8_TOTAL.eps"),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  par(xpd=TRUE)
  plot(PercParental_CD8$percParental_CD8,PercParental_UNKNW$percParental_U8,
       bg=jPalette[as.character(PercParental_CD8$sample)],
       col='gray40',
       pch=21,
       ylab='CD8+ (1&2) % of total area',
       xlab='UNKN (8) % of total area',
       
       bty='l',
       las=1)
  
  cr <- cor.test(PercParental_CD8$percParental_CD8,PercParental_UNKNW$percParental_U8)
  test <- cbind.data.frame(y=PercParental_UNKNW$percParental_U8,x=PercParental_CD8$percParental_CD8)
  lr <- lm(y~x,test)
  xrng<-data.frame(x=par('usr')[c(1,2)])
  yrng<-predict(lr,xrng)
  
  lines(x=unlist(xrng),y=yrng,col='red',lty=3,lwd=3)
  legend('topleft',
         legend = c(paste0('Correlation = ',
                           round(cr$estimate,3)),
                    paste0('Correlation p val = ',
                           formatC(cr$p.value,format = 'e', digits = 2),
                           ' (',gtools::stars.pval(cr$p.value),')')),
         
         bty = 'n',
         cex=0.8,
         inset=c(-0.05,-0.1))
  dev.off()
  
  
  
  
  