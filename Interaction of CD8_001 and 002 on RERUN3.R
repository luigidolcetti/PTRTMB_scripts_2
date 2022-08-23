library(sf)
primary_id<-c('CD4','CD8','CD20','EPTL','MYLD','SMA','UKNW')
prettyLabels<-c('CD4+ cells','CD8+ cells','CD20+ cells','E-cadherin+ cells','Myeloid cells',
                'aSMA+ cells','Unclassified cells','Dump cells')
grPalette<- setNames(RColorBrewer::brewer.pal(length(prettyLabels),name = 'Set2'),prettyLabels)

names(prettyLabels)<-c(primary_id,'DUMP')
names(primary_id)<-prettyLabels[-8]
geom<- st_read("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003/bigTables/fullSet_biexp_trans.sqlite")
expGate <- read.csv("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/bigTables/expGate.csv",
                    row.names = NULL)
expGate <- expGate[!expGate$dump,]
geom <- geom[geom$pids %in% expGate$pIDs,]

expGate$uMap_som <- formatC(expGate$uMap_som,width = 3,flag='0')
extracolumn <- paste(expGate$primary_id,expGate$uMap_som,sep='_')
expGate <- cbind.data.frame(expGate,data.frame(class = extracolumn))
expGate$class <- factor(expGate$class)
uids <- unique(expGate$uid)

intr_CD8_001 <- lapply(setNames(uids,uids),function(u, origin = c('CD8_001','CD8_002')){
  subG <- expGate[expGate$uid == u,]
  subSF <- geom[geom$uid == u,]
  
  subG_focal <- subG[subG$class %in% origin,]
  if (nrow(subG_focal)==0) return(NULL)
  rlt <- st_relate(subSF[subSF$pids %in% subG_focal$pIDs,],subSF,pattern = "F***T****")
  #intr <- st_intersects(subSF[subSF$pids %in% subG_focal$pIDs,],subSF)
  out<- lapply(rlt,function(x)table(subG$class[x]))
  out<- do.call(rbind,out)
  out<- colMeans(out)
  return(out)
  })

whichOut <- unlist(lapply(intr_CD8_001,is.null))

postscript(file = file.path("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/interBar.eps"),
           onefile = F,
           width = 8,
           height = 15,
           horizontal = F,
           paper = 'special',
           bg='white')
meanInter <- do.call(rbind,intr_CD8_001[!whichOut])
meanI <- colMeans(meanInter)
sdI <- apply(meanInter,2,sd)/sqrt(nrow(meanInter))

meanI <- meanI[order(meanI)]
sdI <- sdI[order(sdI)]
newNames <- strsplit(names(meanI),'_')
newColors <- unlist(lapply(newNames,'[',1))
newColors <- grPalette[prettyLabels[newColors]]
newCodes <- unlist(lapply(newNames, '[',2))
names(meanI)<-newCodes
x <- barplot(meanI,horiz = TRUE,names.arg = names(meanI),las=1,cex.names = 0.8,
             xlim=c(0,max(meanI+sdI)),col=newColors,xlab='mean number of interactions')
legend('bottomright',
       legend = names(grPalette)[-8],
       pch=21,
       pt.bg = grPalette[-8],
       col ='black',
       cex=1.5,
       inset = c(0,0),
       bty = 'n')
for (i in 1:length(x)){
  lines(matrix(c(meanI[i]+sdI[i],
                 meanI[i]-sdI[i],
                 x[i],x[i]),
               ncol = 2,
               byrow = FALSE))
}
dev.off()
