primary_id<-c('CD4','CD8','CD20','EPTL','MYLD','SMA','UKNW')
prettyLabels<-c('CD4+ cells','CD8+ cells','CD20+ cells','E-cadherin+ cells','Myeloid cells',
                'aSMA+ cells','Unclassified cells','Dump cells')
names(prettyLabels)<-c(primary_id,'DUMP')
names(primary_id)<-prettyLabels[-8]
grPalette<- setNames(RColorBrewer::brewer.pal(length(primary_id),name = 'Set2'),primary_id)

studyTable<-read.csv("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003/bigTables/studyTable.csv")
expGate <- read.csv("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/bigTables/expGate.csv",
                    row.names = NULL)
areaS <- sf::read_sf(
  "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003/bigTables/fullSet_biexp_trans.sqlite"
)
colnames(areaS)[1] <- 'pIDs'
expGate <- expGate[!expGate$dump,]
expGate <- dplyr::inner_join(expGate,areaS[,c('pIDs','area'),drop=TRUE],by='pIDs')

samples <- unique(expGate$sample)
pfSample <- unique(studyTable[,c('sample','bioGroup')])
pfSample <- setNames(pfSample$bioGroup,pfSample$sample)
pfSample <- pfSample[order(pfSample)]

# postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/XsampleUmap.eps"),
#            onefile = F,
#            width = 7.5,
#            height = 7.5,
#            horizontal = F,
#            paper = 'special',
#            bg='white')

tiff(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/XsampleUmap_sideBar.tif"),
           
           width = 7.5,
           height = 7.5,
           units = 'in',
     compression = 'lzw',
     res = 600,
           bg='white')

par(mfrow=c(3,3),mar=c(1,1,1,1),xpd=TRUE)

for (i in 1:length(pfSample)){ 
  
  subPop <- expGate[expGate$sample==names(pfSample)[i],]
  popTot <- aggregate(area~sample,subPop,sum)
  popProportion <- aggregate(area~sample+primary_id,subPop,sum)
  popProportion$area <- popProportion$area/popTot$area
  popAbsent <- primary_id[!(primary_id %in% popProportion$primary_id)]
  if (length(popAbsent)!=0){
  popAbsent <- data.frame(sample = unique(popProportion$sample),
                          primary_id = popAbsent,
                          area = 0)
  popProportion<- rbind.data.frame(popProportion,popAbsent)}
  plot(expGate$uMAP1,expGate$uMAP2,
       asp=1,
       cex=0.3,
       pch=20,
       xaxt='n',
       yaxt='n')
  limS <- par('usr')
  boxStart <-limS[4]
  for (ii in 1:length(primary_id)){
    
    boxHeight <-(limS[4]-limS[3])*popProportion$area[popProportion$primary_id==primary_id[ii]]
    rect(limS[2]+(limS[2]-limS[1])/15,boxStart,limS[2],boxStart-boxHeight, col = grPalette[primary_id[ii]])
    boxStart <- boxStart -boxHeight
  }
  axis(1,at = c(-10,-5,0,5,10,15),labels = FALSE)
  axis(2,at = c(-10,-5,0,5,10),labels = FALSE)
  mtext(paste0(names(pfSample)[i],' - pfs: ',pfSample[i]),side = 3,line = 0,adj = 0.5)
  XYarray<-expGate[expGate$sample==names(pfSample)[i],c('uMAP1','uMAP2')]
  TEMPD<-ks::kde(XYarray)
  
  newZ<-predict(object = TEMPD,x=XYarray)
  
  newZ<-cut(newZ,
            include.lowest = T,
            breaks = c(-Inf,TEMPD$cont,+Inf),
            labels=F)
  colZ<-setNames(1:max(newZ),colorRampPalette(c('cadetblue','aquamarine','brown2'))(max(newZ)))
  newZ<-names(colZ)[newZ]
  points(XYarray,
         cex=0.1,pch=20,col=newZ)
}
dev.off()
########
########
tiff(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/XsampleUmap_sideBar_points_only.tif"),
     
     width = 7.5,
     height = 7.5,
     units = 'in',
     compression = 'lzw',
     res = 600,
     bg='white')

par(mfrow=c(3,3),mar=c(1,1,1,1),xpd=TRUE)

for (i in 1:length(pfSample)){ 
  
  subPop <- expGate[expGate$sample==names(pfSample)[i],]
  popTot <- aggregate(area~sample,subPop,sum)
  popProportion <- aggregate(area~sample+primary_id,subPop,sum)
  popProportion$area <- popProportion$area/popTot$area
  popAbsent <- primary_id[!(primary_id %in% popProportion$primary_id)]
  if (length(popAbsent)!=0){
    popAbsent <- data.frame(sample = unique(popProportion$sample),
                            primary_id = popAbsent,
                            area = 0)
    popProportion<- rbind.data.frame(popProportion,popAbsent)}
  plot(expGate$uMAP1,expGate$uMAP2,
       asp=1,
       cex=0.3,
       pch=20,
       xaxt='n',
       yaxt='n',
       bty = 'n')
  limS <- par('usr')
  boxStart <-limS[4]
  # for (ii in 1:length(primary_id)){
  #   
  #   boxHeight <-(limS[4]-limS[3])*popProportion$area[popProportion$primary_id==primary_id[ii]]
  #   rect(limS[2]+(limS[2]-limS[1])/15,boxStart,limS[2],boxStart-boxHeight, col = grPalette[primary_id[ii]])
  #   boxStart <- boxStart -boxHeight
  # }
  # axis(1,at = c(-10,-5,0,5,10,15),labels = FALSE)
  # axis(2,at = c(-10,-5,0,5,10),labels = FALSE)
  # mtext(paste0(names(pfSample)[i],' - pfs: ',pfSample[i]),side = 3,line = 0,adj = 0.5)
  XYarray<-expGate[expGate$sample==names(pfSample)[i],c('uMAP1','uMAP2')]
  TEMPD<-ks::kde(XYarray)
  
  newZ<-predict(object = TEMPD,x=XYarray)
  
  newZ<-cut(newZ,
            include.lowest = T,
            breaks = c(-Inf,TEMPD$cont,+Inf),
            labels=F)
  colZ<-setNames(1:max(newZ),colorRampPalette(c('cadetblue','aquamarine','brown2'))(max(newZ)))
  newZ<-names(colZ)[newZ]
  points(XYarray,
         cex=0.1,pch=20,col=newZ)
}
dev.off()

postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/XsampleUmap_only_frame.eps"),
           onefile = F,
           width = 7.5,
           height = 7.5,
           horizontal = F,
           paper = 'special',
           bg='white')

par(mfrow=c(3,3),mar=c(1,1,1,1),xpd=TRUE)

for (i in 1:length(pfSample)){ 
  
  subPop <- expGate[expGate$sample==names(pfSample)[i],]
  popTot <- aggregate(area~sample,subPop,sum)
  popProportion <- aggregate(area~sample+primary_id,subPop,sum)
  popProportion$area <- popProportion$area/popTot$area
  popAbsent <- primary_id[!(primary_id %in% popProportion$primary_id)]
  if (length(popAbsent)!=0){
    popAbsent <- data.frame(sample = unique(popProportion$sample),
                            primary_id = popAbsent,
                            area = 0)
    popProportion<- rbind.data.frame(popProportion,popAbsent)}
  plot(NA,
       xlim = limS[1:2],
       ylim = limS[3:4],
       asp=1,
       cex=0.3,
       pch=20,
       xaxt='n',
       yaxt='n')
  limS <- par('usr')
  boxStart <-limS[4]
  for (ii in 1:length(primary_id)){
    
    boxHeight <-(limS[4]-limS[3])*popProportion$area[popProportion$primary_id==primary_id[ii]]
    rect(limS[2]+(limS[2]-limS[1])/15,boxStart,limS[2],boxStart-boxHeight, col = grPalette[primary_id[ii]])
    boxStart <- boxStart -boxHeight
  }
  axis(1,at = c(-10,-5,0,5,10,15),labels = FALSE)
  axis(2,at = c(-10,-5,0,5,10),labels = FALSE)
  mtext(paste0(names(pfSample)[i],' - pfs: ',pfSample[i]),side = 3,line = 0,adj = 0.5)
  # XYarray<-expGate[expGate$sample==names(pfSample)[i],c('uMAP1','uMAP2')]
  # TEMPD<-ks::kde(XYarray)
  # 
  # newZ<-predict(object = TEMPD,x=XYarray)
  # 
  # newZ<-cut(newZ,
  #           include.lowest = T,
  #           breaks = c(-Inf,TEMPD$cont,+Inf),
  #           labels=F)
  # colZ<-setNames(1:max(newZ),colorRampPalette(c('cadetblue','aquamarine','brown2'))(max(newZ)))
  # newZ<-names(colZ)[newZ]
  # points(XYarray,
  #        cex=0.1,pch=20,col=newZ)
}
dev.off()









postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/XsampleUmap_legend.eps"),
           onefile = F,
           width = 7.5,
           height = 7.5,
           horizontal = F,
           paper = 'special',
           bg='white')
b=5
plot(NA,xlim=c(0,1),ylim=c(0,ceiling(length(TEMPD$cont)/b)+1),asp=1,bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
ii=0
for(i in seq(0,length(colZ),by=b)){
  rect(0,ii,1,ii+1,col = c('black',names(colZ))[i+1],border=NA)
  text(-0.2,ii+0.5,c('absent',colZ)[i+1],adj=c(1,0.5))
ii=ii+1
  }
rect(0,0,1,ceiling(length(TEMPD$cont)/b)+1)
dev.off()
