expGate <- read.csv(
  "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/bigTables/expGate.csv"
)


primary_id<-c('CD4','CD8','CD20','EPTL','MYLD','SMA','UKNW')
prettyLabels<-c('CD4+ cells','CD8+ cells','CD20+ cells','E-cadherin+ cells','Myeloid cells',
                'aSMA+ cells','Mixed phenotype','Dump cells')
names(prettyLabels)<-c(primary_id,'DUMP')
names(primary_id)<-prettyLabels[-8]
grPalette<- setNames(RColorBrewer::brewer.pal(length(prettyLabels),name = 'Set2'),prettyLabels)

tiff(file = file.path(
  "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/def",
  'uMap_points.tif'),
  res = 600,
  units = 'in',
  width = 7.5,
  height = 7.5)

plot(expGate[,c('uMAP1','uMAP2')],
     pch=16,
     cex=0.5,
     col='black',main='',
     xlab='',ylab='',xaxt='n',yaxt='n',bty='n',asp=1)
lims <- par('usr')

for (i in primary_id){
  points(expGate[expGate$primary_id==i & !expGate$dump,c('uMAP1','uMAP2')],
         pch=16,
         cex=0.1,
         col=grPalette[names(primary_id)[primary_id==i]])
  
}
dev.off()



postscript(file = file.path(
  "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/EXTRA/def",
  'uMap_frame.eps'),
  onefile = F,
  width = 7.5,
  height = 7.5,
  horizontal = F,
  paper = 'special',
  bg='white')

plot(NA,
     pch=16,
     cex=0.5,
     col='black',
     main = 'Total cells',
     xlim = lims[1:2],ylim=lims[3:4],
     xlab='uMAP1',ylab='uMAP2',asp=1)

legend('topright',
       legend = names(grPalette)[-8],
       pch=21,
       pt.bg = grPalette[-8],
       col ='black',
       cex=1,
       inset = c(0,0),
       bty = 'n')

dev.off()
