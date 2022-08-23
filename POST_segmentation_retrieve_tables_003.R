library(sf)
library(FlowSOM)
library(RUNIMCTEMP)
library(pbapply)

##### declare vars ####

rootFolder<-"C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4"
analFolder<- "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003"

dirList<-list.dirs(rootFolder,recursive = F)
whichDir<-grepl('RUN_s',dirList,fixed = T)
dirList<-dirList[whichDir]
fileList<-paste0(dirList,"/analysis/Anal_001/archive/Expression.RDS")

trnsf_mod<- flowCore::biexponentialTransform(a=0.3,b=0.01,c=0.1,d=5,f=-0.3,w=10)

primary_id<-c('CD4','CD8','CD20','EPTL','MYLD','SMA','UKNW')
prettyLabels<-c('CD4+ cells','CD8+ cells','CD20+ cells','E-cadherin+ cells','Myeloid cells',
                'aSMA+ cells','Unclassified cells','Dump cells')
names(prettyLabels)<-c(primary_id,'DUMP')
names(primary_id)<-prettyLabels[-8]

##### load single files ####

expTable<- pblapply(fileList,function(fl){
  
  objectTotal<-readRDS(fl)
  tableTotal<-objectTotal@exprs$toplayers
  valueTotal<-objectTotal@exprs$toplayers_mean[,-1]
  valueTotal<-apply(valueTotal,2,trnsf_mod)
  idsTotal<-objectTotal@exprs$toplayers_mean[,1,drop=F]
  valueTotal<-cbind.data.frame(idsTotal,valueTotal)
  attr(valueTotal[,1],'class')<-NULL
  
  uids<-unique(tableTotal$uid)
  
  cleanTable<-lapply(uids,function(uids){
    tableSub<-tableTotal[tableTotal$uid==uids,]
    
    
    intrct<-st_intersects(tableSub)
    notBanale<-unlist(lapply(intrct,function(x)length(x)>1))
    intrct<-intrct[notBanale]
    
    intrct<-lapply(intrct,function(x){
      origin<-x[1]
      target<-x[-1]
      intrct<-st_relate(tableSub[origin,],tableSub[target,])
      out<-matrix(c(rep(origin,length(target)),target,intrct),ncol=3,byrow = F)
      out<-out[regexpr("^2",out[,3])==1,1:2]
      return(out)
    })
    intrct<-do.call(rbind,intrct)
    intrct<-unique(intrct)
    
    dup1<-intrct[duplicated(intrct[,1]),1]
    dup2<-intrct[duplicated(intrct[,2]),2]
    dup3<-intrct[!(intrct[,1] %in% dup1) & !(intrct[,2] %in% dup2),,drop=F]
    
    dup3<-apply(dup3,1,function(x){
      a1<-st_area(tableSub[x[1],])
      a2<-st_area(tableSub[x[2],])
      x[which.min(c(a1,a2))]
    })
    
    outPIDs<-tableSub$pIDs[as.numeric(c(dup1,dup2,dup3)),drop=T]
    tableSub<-tableSub[!(tableSub$pIDs %in% outPIDs),]
    
    out<-dplyr::left_join(tableSub,valueTotal,by='pIDs')
    out$area<-st_area(out)
    
    return(out)
  })
  cleanTable<-do.call(dplyr::bind_rows,cleanTable)
  return(cleanTable)
})

expData<-do.call(dplyr::bind_rows,expTable)

dir.create(file.path(analFolder,'bigTables'))

st_write(expData,
         file.path(analFolder,'bigTables','fullSet_biexp_trans.sqlite'),
         append=F)

##### in case of restart from here
expData<-st_read(file.path(analFolder,'bigTables','fullSet_biexp_trans.sqlite'))
colnames(expData)[1]<-'pIDs'
colnames(expData)[40]<-'geom'
st_geometry(expData)<-'geom'
#######

rm(expTable)

gc()

##### retrieve study table #####

fileList<-paste0(dirList,"/archive/studyTable.RDS")

studyTable<-lapply(fileList,function(fl){
  out<-readRDS(fl)
  return(out)
})

studyTable<-do.call(rbind.data.frame,studyTable)

write.csv(studyTable,file.path(analFolder,'bigTables','studyTable.csv'),row.names = F)

expData<- dplyr::full_join(expData,studyTable,by="uid")

expGate<-data.frame(pIDs = expData$pIDs,
                    uid = expData$uid,
                    sample = expData$sample,
                    primary_id = as.factor(expData$primary_id),
                    primary_som = factor(NA),
                    uMap = NA,
                    uMap_som = factor(NA),
                    uMAP1 = NA,
                    uMAP2 = NA,
                    dump = NA)

##### initial som on total ####

expMatrix<-expData[,10:39,drop=T]
expMatrix<-as.matrix(expMatrix)
oldCnames<-colnames(expMatrix)
chTable<-readRDS("C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4/RUN_test/archive/ChannelTable.RDS")
newCnames<-chTable$marker[chTable$RcolumnNames %in% oldCnames]
names(newCnames)<-oldCnames
colnames(expMatrix)<-newCnames[colnames(expMatrix)]
colToUse<-colnames(expMatrix)[-c(10,20,29,30)]

som<-ReadInput(input = expMatrix,
               compensate = F,
               transform = F,
               scale = T,
               silent = F)

som<-BuildSOM(fsom = som,
              colsToUse = colToUse,
              xdim = 10,
              ydim = 10)

dir.create(file.path(analFolder,'primary_som'))

saveRDS(som,
        file.path(analFolder,'primary_som',"som.rds"))

somClass<-GetClusters(som)

somClass<-formatC(somClass,width=3,flag='0')

somClass<-factor(somClass)
somClass<-addNA(somClass)

expGate$primary_som<-somClass
rm(somClass)
rm(som)

TEMP<-dplyr::left_join(expData[,c('pIDs','area'),drop=T],
                       expGate[,c('pIDs','primary_id','primary_som')],
                       by ='pIDs')

agg_mixClass<-aggregate(area~primary_id+primary_som,TEMP,sum)

rm(TEMP)

agg_mixClass$primary_id<-factor(agg_mixClass$primary_id)
agg_mixClass$primary_som<-factor(agg_mixClass$primary_som)

agg_mixClass<-lapply(setNames(levels(agg_mixClass$primary_id),levels(agg_mixClass$primary_id)),function(x){
  out<-lapply(setNames(levels(agg_mixClass$primary_som),levels(agg_mixClass$primary_som)),function(y){
    out<-agg_mixClass[agg_mixClass$primary_id==x & agg_mixClass$primary_som==y,'area']
    if (length(out)==0) out<-0
    return(out)
  })
  out<-matrix(unlist(out),nrow = 1,dimnames = list(Pixel_level=x,Cell_level=levels(agg_mixClass$primary_som)))
})
agg_mixClass<-do.call(rbind,agg_mixClass)

agg_mixClass_scale<-scale(agg_mixClass)
dst<-dist(t(agg_mixClass_scale),method = 'maximum')
hc<-hclust(dst,method = 'complete')

TEMP<-lapply(1:ncol(dst),function(kats){
  gr<-cutree(hc,k=kats)
  gr<-gr[hc$order]
  crossNames<-unlist(lapply(unique(gr),function(x){
    somClass<-names(gr[gr==x])
    subMat<-agg_mixClass_scale[,somClass,drop=F]
    winner<-apply(subMat,1,sum)
    if (length(which(winner>0))>1) {
      out<-setNames('DUMP',x)
      return(out)
    } else {
      winner<-names(which.max(winner))
      out<-setNames(winner,x)
      return(out)}
  }))
  out<-table(crossNames)
  out<-c(sum(out[names(out)!="DUMP"])/length(unique(expGate$primary_id)),out['DUMP'])
  out[is.na(out)]<-0
  names(out)<-c('fragmentation','DUMP')
  return(out)
})
TEMP<-do.call(rbind,TEMP)

candidateDUMP<-which(abs(1-TEMP[,1])==min(abs(1-TEMP[,1])))
toDUMP<-max(TEMP[candidateDUMP,2])
rangeDUMP<-which(TEMP[,2]==toDUMP)
DumpIt<-rangeDUMP[rangeDUMP %in% candidateDUMP]

Kats<-DumpIt

gr<-cutree(hc,k=Kats)
gr<-gr[hc$order]

crossNames<-unlist(lapply(unique(gr),function(x){
  somClass<-names(gr[gr==x])
  subMat<-agg_mixClass_scale[,somClass,drop=F]
  winner<-apply(subMat,1,sum)
  if (length(which(winner>0))>1) {
    out<-setNames('DUMP',x)
    return(out)
  } else {
    winner<-names(which.max(winner))
    out<-setNames(winner,x)
    return(out)}
}))

grNames<-prettyLabels[crossNames[as.character(gr)]]
names(grNames)<-names(gr)

grPalette<- setNames(RColorBrewer::brewer.pal(length(unique(grNames)),name = 'Set2'),prettyLabels)
grAnnotation<-data.frame(Composite=as.factor(grNames))
grAnnotation_color<-list(Composite=grPalette[unique(grNames)])

##### heat map som ####

pheatmap::pheatmap(t(agg_mixClass_scale),
                   cluster_rows = T,
                   cluster_cols = T,
                   cutree_rows = Kats,
                   clustering_distance_rows = 'maximum',
                   clustering_method = 'complete',
                   annotation_row = grAnnotation,
                   annotation_colors = grAnnotation_color,
                   labels_col = prettyLabels[colnames(t(agg_mixClass_scale))],
                   filename = file.path(analFolder,'primary_som','somMFIxforest.pdf'),
                   cellwidth = 13,
                   cellheight =13,
                   height = 22,
                   width = 6,
                   main = 'Coposite clusters',
                   scale = 'none')


##### calculate compound class  and uMAP#####

clusterToDump<-crossNames[as.character(gr)]
names(clusterToDump)<-names(gr)

expGate$dump<-!(apply(apply(cbind(clusterToDump,names(clusterToDump)),1,function(x)expGate$primary_id==x[1] & expGate$primary_som==x[2] & x!='DUMP'),1,any))

ump<-umap::umap(expData[expData$pIDs %in% expGate$pIDs[!expGate$dump],9:39,drop=T])

dir.create(file.path(analFolder,'UMAP'))

saveRDS(ump,file.path(analFolder,'UMAP','umap_model.R'))

ump_total<- predict(ump,expData[,9:39,drop=T])

expGate$uMAP1<-ump_total[,1]
expGate$uMAP2<-ump_total[,2]

rm(ump)
rm(ump_total)

h=20
m1min<-min(expGate[,'uMAP1'])
m1max<-max(expGate[,'uMAP1'])
mm1<-m1max-m1min
m2min<-min(expGate[,'uMAP2'])
m2max<-max(expGate[,'uMAP2'])
mm2<-m2max-m2min

dumpOnU<-lapply(seq(m1min,m1max-mm1/h,mm1/h),function(x){
  out<-lapply(seq(m2min,m2max-mm2/h,mm2/h),function(y){
    pID_nghbrs<-expGate$pIDs[expGate$uMAP1>=x &
                               expGate$uMAP1<=x+mm1/h&
                               expGate$uMAP2>=y &
                               expGate$uMAP2<=y+mm2/h&
                               !expGate$dump
    ]
    if(length(pID_nghbrs)==0) {
      out<-data.frame(pIDs=character(0),dump=logical(0))} else {
        primary_id<-expGate$primary_id[expGate$pIDs %in% pID_nghbrs]
        winner_id<-table(primary_id)
        winner_id<-names(winner_id)[which.max(winner_id)]
        out<-data.frame(pIDs=pID_nghbrs,dump=(primary_id!=winner_id))
      }
    
    return(out)
  })
  out<-do.call(rbind.data.frame,out)
})

dumpOnU<-do.call(rbind.data.frame,dumpOnU)
expGate$dump[match(dumpOnU$pIDs,expGate$pIDs)]<-dumpOnU$dump
expGate$uMap[match(dumpOnU$pIDs,expGate$pIDs)]<-dumpOnU$dump

postscript(file = file.path(analFolder,'Umap','uMap_Reteined.eps'),
           onefile = F,
           width = 7.5,
           height = 7.5,
           horizontal = F,
           paper = 'special',
           bg='white')

plot(expGate[,c('uMAP1','uMAP2')],
     pch=16,
     cex=0.5,
     col='black',
     main = 'Valid cells',
     xlab='uMAP1',ylab='uMAP2')

legend('bottomright',
       legend = names(grPalette)[-8],
       pch=21,
       pt.bg = grPalette[-8],
       col ='black',
       cex=1,
       inset = c(0,0),
       bty = 'n')

for (i in primary_id){
  points(expGate[expGate$primary_id==i & !expGate$dump,c('uMAP1','uMAP2')],
         pch=16,
         cex=0.1,
         col=grPalette[names(primary_id)[primary_id==i]])
  
}
dev.off()

postscript(file = file.path(analFolder,'Umap','uMap_ReJected.eps'),
           onefile = F,
           width = 7.5,
           height = 7.5,
           horizontal = F,
           paper = 'special',
           bg='white')

plot(expGate[,c('uMAP1','uMAP2')],
     pch=16,
     cex=0.5,
     col='black',
     main = 'Rejected cells',
     xlab='uMAP1',ylab='uMAP2')

legend('bottomright',
       legend = names(grPalette)[-8],
       pch=21,
       pt.bg = grPalette[-8],
       col ='black',
       cex=1,
       inset = c(0,0),
       bty = 'n')

for (i in primary_id){
  points(expGate[expGate$primary_id==i & expGate$dump,c('uMAP1','uMAP2')],
         pch=16,
         cex=0.1,
         col=grPalette[names(primary_id)[primary_id==i]])
  
}
dev.off()

###### detailed heatmaps ######

#### uncomment next 2 in case 
# importanceTable<-matrix(0,
#                         ncol = ncol(expMatrix),
#                         nrow = length(unique(expData$primary_id)),
#                         dimnames = list(cluster=sort(unique( expData$primary_id )),
#                                         marker=colnames(expMatrix)))

# write.csv(importanceTable,
#           file.path(analFolder,'bigTables','ImportanceTable.csv'))

importanceTable<-read.csv(file.path(analFolder,'bigTables','importanceTable.csv'),
                          header = T,
                          row.names = 1,
                          check.names = F)
importanceTable<-cbind.data.frame(importanceTable,
                                  uMAP1 = rep(1,nrow(importanceTable)),
                                  uMAP2 = rep(1,nrow(importanceTable)))

mrkr<-unique(expData$primary_id)
colToUse<-colnames(importanceTable)
dir.create(file.path(analFolder,'subClusters'))
for (i in 1:length(mrkr)){
  
  subMatrix<-expMatrix[!expGate$dump & expGate$primary_id==mrkr[i],]
  subMatrix<-cbind(subMatrix,expGate[!expGate$dump & expGate$primary_id==mrkr[i],c('uMAP1','uMAP2')])
  subMatrix<-as.matrix(subMatrix)
  
  som<-ReadInput(input = subMatrix,
                 compensate = F,
                 transform = F,
                 scale = T,
                 silent = F)
  
  som<-BuildSOM(fsom = som,
                colsToUse = colnames(subMatrix),
                importance = unlist(importanceTable[mrkr[i],,drop=T]),
                xdim = 4,
                ydim = 2)
  
  saveRDS(som,
          file.path(analFolder,'subClusters',paste0("som_",mrkr[i],".rds")))
  
  somClass<-GetClusters(som)
  
  somClass<-formatC(somClass,width=3,flag='0')
  
  somClass<-factor(somClass)
  somClass<-addNA(somClass)
  levels(expGate$uMap_som)<-c(unique(levels(somClass),levels(expGate$uMap_som)))
  expGate$uMap_som[!expGate$dump & expGate$primary_id==mrkr[i]]<-somClass
  
  somMFI<-dplyr::inner_join(expData[,c(1,2,4,9:39),drop=T],expGate[!expGate$dump & expGate$primary_id==mrkr[i],c('pIDs','uMap_som')],by='pIDs')
  
  colnames(somMFI)[match(names(newCnames),colnames(somMFI))]<-newCnames
  
  colToUse<-c(setNames("area","area"),newCnames)
  
  somMFI_median<-aggregate(somMFI[,colToUse],list(somMFI$uMap_som),median)
  TEMP<-somMFI_median[,1]
  somMFI_median<-as.matrix(somMFI_median[,-1])
  rownames(somMFI_median)<-TEMP
  
  colmin<-apply(expData[!expGate$dump,c(9:39),drop=T],2,quantile,0.01)
  names(colmin)<-colToUse
  colmax<-apply(expData[!expGate$dump,c(9:39),drop=T],2,quantile,0.99)
  names(colmax)<-colToUse
  somMFI_median<-do.call(cbind,lapply(setNames(colnames(somMFI_median),colnames(somMFI_median))
                                      ,function(x){
                                        do.call(c,lapply(setNames(rownames(somMFI_median),rownames(somMFI_median)),function(y){
                                          if (somMFI_median[y,x]>=colmax[x]) return(colmax[x])
                                          if (somMFI_median[y,x]<=colmin[x]) return(colmin[x])
                                          return(somMFI_median[y,x])
                                        }))
                                      }))
  
  somMFI_median<-do.call(cbind,lapply(setNames(colnames(somMFI_median),colnames(somMFI_median))
                                      ,function(x)(somMFI_median[,x]-colmin[x])/(colmax[x]-colmin[x])))
  
  pheatmap::pheatmap(t(somMFI_median),
                     cluster_rows = F,
                     cluster_cols = T,
                     filename = file.path(analFolder,'subClusters',paste0("som_",mrkr[i],".pdf")),
                     cellwidth = 13,
                     cellheight =13,
                     height = 8,
                     width = 5,
                     scale = 'none',
                     main = mrkr[i])
  
}

# saveRDS(somList,file.path(analFolder,'subClusters','Sub_cluster_Som_list.R'))

write.csv(expGate,
          file.path(analFolder,'bigTables','expGate.csv'),
          row.names = F)


for (i in 1:length(mrkr)){
  
  
  uPalette<- setNames(c(RColorBrewer::brewer.pal(length(levels(expGate$uMap_som)),name = 'Set2')),
                      levels(expGate$uMap_som))
  
  postscript(file = file.path(analFolder,'subClusters',paste0('Umap_',mrkr[i],'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  
  plot(expGate[expGate$primary_id==mrkr[i] & !expGate$dump,c('uMAP1','uMAP2')],
       pch=16,
       cex=0.6,
       col='black',
       main = prettyLabels[mrkr[i]],
       xlab='uMAP1',ylab='uMAP2')
  
  legend('bottomleft',
         legend = names(uPalette),
         pch=21,
         pt.bg = uPalette,
         col = 'black',
         cex=1,
         inset = c(0,0),
         bty = 'n')
  
  for (ii in unique(expGate$uMap_som)){
    
    points(expGate[expGate$primary_id==mrkr[i] & expGate$uMap_som==ii&
                     !expGate$dump ,c('uMAP1','uMAP2')] ,
           pch=16,
           cex=0.1,
           col=uPalette[ii],
           xlab='uMAP1',ylab='uMAP2')
  }
  dev.off()
}


##### SUPER POP percentages ####
##### VS sample ####
colToKeep_data<-c('uid','sample','bioGroup','pIDs','area')
colToKeep_gate<-c('pIDs','primary_id','uMap_som')

aggArea<-lapply(setNames(unique(expGate$uid),unique(expGate$uid)),function(u){
  expSubSet<-dplyr::left_join(expGate[expGate$uid == u & !expGate$dump,colToKeep_gate],
                              expData[,colToKeep_data,drop=T])
  
  sampleArea<-sum(expSubSet$area)
  primary_idArea<-by(expSubSet$area,list(expSubSet$primary_id),sum)
  uMap_somArea<-by(expSubSet$area,list(expSubSet$uMap_som,expSubSet$primary_id),sum)
  
  names_primary_idArea<-dimnames(primary_idArea)[[1]]
  names_uMap_somArea<-dimnames(uMap_somArea)[[1]]
  
  out<-data.frame(uid=u,
                  sample = unique(expSubSet$sample),
                  sample.area = sampleArea,
                  super.cluster = rep(names_primary_idArea,each=length(names_uMap_somArea)),
                  super.area = rep(as.vector(primary_idArea),each=length(names_uMap_somArea)),
                  sub.cluster = names_uMap_somArea,
                  sub.area = as.vector(uMap_somArea),
                  pfs = unique(expSubSet$bioGroup))
  out[is.na(out)]<-0
  
  return(out)})

aggArea<-do.call(rbind.data.frame,aggArea)
rownames(aggArea)<-NULL

write.csv(aggArea,
          file.path(analFolder,'bigTables','areaTable.csv'),
          row.names = F)

percentage_super<-aggArea[,-c(6,7)]
percentage_super<-unique(percentage_super)

av_DF<-data.frame(pfs=percentage_super$pfs,
                  smp=percentage_super$sample,
                  class=percentage_super$super.cluster,
                  area=percentage_super$super.area/percentage_super$sample.area*100)

av_list_mp<-lapply(setNames(unique(av_DF$class),unique(av_DF$class)),function(x){
  av<-aov(area~smp,av_DF[av_DF$class==x,])
  
  av_HSD<-TukeyHSD(av)
  
  out<-data.frame(av_first_smp = unlist(lapply(strsplit(rownames(av_HSD$smp),'-'),'[',1)),
                  av_first_som = x,
                  av_second_smp = unlist(lapply(strsplit(rownames(av_HSD$smp),'-'),'[',2)),
                  av_second_som = x,
                  p.adj= av_HSD$smp[,"p adj"],
                  s.adj= gtools::stars.pval(av_HSD$smp[,"p adj"]),
                  row.names = NULL)
  return(out)
})

av_list_mp<-do.call(rbind.data.frame,av_list_mp)

smpOrder<-unique(percentage_super[,c('sample','pfs')])
smpOrder<-factor(smpOrder$sample[order(as.numeric(smpOrder$pfs),decreasing = F)],
                 levels=smpOrder$sample[order(as.numeric(smpOrder$pfs),decreasing = F)],
                 ordered = T)

dir.create(file.path(analFolder,"SUPER.POP"))
for (i in unique(av_DF$class)){
  subSet<-av_DF$area[av_DF$class==i]
  names(subSet)<-av_DF$smp[av_DF$class==i]
  smp<-factor(names(subSet),levels = levels(smpOrder),ordered = T)
  df_box<-data.frame(Sample=smp,
                     `% of total area` = subSet,
                     check.names = F)
  postscript(file = file.path(analFolder,"SUPER.POP",paste0(i,'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  
  par(mfrow=c(2,1),mar=c(0,5,2,1),las=1,bty='l')
  
  boxplot(`% of total area`~Sample,
          data = df_box,
          main=prettyLabels[i],
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n')
  points(data.frame(jitter(as.integer(df_box$Sample)),df_box$`% of total area`),bg='gray60',col='gray40',pch=21)
  boxplot(`% of total area`~Sample,
          data = df_box,
          main=i,
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n',
          add=T)
  ll<-par('usr')
  plot(NA,
       xlim=c(ll[1],ll[2]),
       ylim=c(ll[2],ll[1]),
       bty='n',
       xlab=NA,
       ylab=NA,
       xaxt='n',
       yaxt='n',
       xaxs='i')
  axis(2,
       at=as.integer(smpOrder),
       labels = paste0("(",smpOrder,") ",LETTERS[1:length(smpOrder)][as.integer(smpOrder)]),
       cex.axis=0.6,ylab=NA,
       col=NA,
       col.ticks = 'black')
  axis(3,
       at=as.integer(smpOrder),
       labels = LETTERS[1:length(smpOrder)],
       cex.axis=0.6,
       col=NA,
       col.ticks = 'black')
  subSet<-av_list_mp[av_list_mp$av_first_som==i,]
  idx<-as.numeric(smpOrder)
  names(idx)<-levels(smpOrder)
  subSet[,"av_first_smp"]<-idx[subSet[,"av_first_smp"]]
  subSet[,"av_second_smp"]<-idx[subSet[,"av_second_smp"]]
  text(subSet[,"av_first_smp"],subSet[,"av_second_smp"],labels=subSet[,"s.adj"])
  for (iii in 1:nrow(subSet)){
    if (subSet$s.adj[iii]!=" " & subSet$s.adj[iii]!="."){
      ff<-subSet$av_first_smp[iii]
      ss<-subSet$av_second_smp[iii]
      
      jtr<-sample(seq(-0.3,0.3,0.01),1)
      sm<-c(ff+jtr,ll[3])
      pm<-c(ff+jtr,ss+jtr)
      em<-c(ll[1],ss+jtr)
      mm<-rbind(sm,pm,em)
      lines(mm,lty=3,col='gray60')
    }
  }
  dev.off()
}


##### VS pfs ####

av_DF<-data.frame(pfs=as.numeric(percentage_super$pfs),
                  smp=percentage_super$sample,
                  class=percentage_super$super.cluster,
                  area=as.numeric(percentage_super$super.area)/as.numeric(percentage_super$sample.area)*100)

av_list_mp<-lapply(setNames(unique(av_DF$class),unique(av_DF$class)),function(x){
  av_lm<-lm(area~pfs,av_DF[av_DF$class==x,])
  av_av<-anova(av_lm)
  av_cor<-cor.test(av_DF$pfs[av_DF$class==x],av_DF$area[av_DF$class==x])
  out<-list(av_som = x,
            p.anova= av_av[1,"Pr(>F)"],
            s.anova= gtools::stars.pval(av_av[1,"Pr(>F)"]),
            lm = av_lm,
            c.cor = av_cor$estimate,
            p.cor = av_cor$p.value,
            s.cor = gtools::stars.pval(av_cor$p.value))
  return(out)
})


dir.create(file.path(analFolder,"SUPER.POP_pfs"))
for (i in unique(av_DF$class)){
  
  df_box<-data.frame(pfs=av_DF$pfs[av_DF$class==i],
                     smp=av_DF$smp[av_DF$class==i],
                     `% of total area` = av_DF$area[av_DF$class==i],
                     check.names = F)
  avr<-aggregate(`% of total area`~smp+pfs,df_box,mean)
  avr<-data.frame(smp = avr$smp,
                  pfs = as.numeric(avr$pfs),
                  mean = as.numeric(avr$`% of total area`))
  sdr<-aggregate(`% of total area`~smp+pfs,df_box,sd)
  sdr<-data.frame(smp = sdr$smp,
                  pfs = as.numeric(sdr$pfs),
                  sd = as.numeric(sdr$`% of total area`))
  avrsdr<- data.frame(smp=avr$smp,
                      pfs=avr$pfs,
                      mean=avr$mean,
                      sd=sdr$sd,
                      `mean+sd`=avr$mean+sdr$sd,
                      `mean-sd`=avr$mean-sdr$sd,
                      check.names = F)
  
  postscript(file = file.path(analFolder,"SUPER.POP_pfs",paste0(i,'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  par(xpd=T)
  plot(data.frame(jitter(as.integer(df_box$pfs)),
                  df_box$`% of total area`),
       bg='gray60',
       col='gray40',
       pch=21,
       ylab='% of parental area',
       xlab='pfs',
       
       bty='l',
       las=1,
       ylim=c(min(c(avrsdr$`mean-sd`,df_box$`% of total area`)),max(avrsdr$`mean+sd`,df_box$`% of total area`)))
  mtext(side = 3,adj=1,prettyLabels[i])
  for (ii in 1:nrow(avrsdr)){
    lines(matrix(c(avrsdr$pfs[ii]-0.5,avrsdr$pfs[ii]+0.5,avrsdr$mean[ii],avrsdr$mean[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(avrsdr$pfs[ii],avrsdr$pfs[ii],avrsdr$`mean+sd`[ii],avrsdr$`mean-sd`[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(avrsdr$pfs[ii]-0.2,avrsdr$pfs[ii]+0.2,avrsdr$`mean+sd`[ii],avrsdr$`mean+sd`[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(avrsdr$pfs[ii]-0.2,avrsdr$pfs[ii]+0.2,avrsdr$`mean-sd`[ii],avrsdr$`mean-sd`[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
  }
  
  xrng<-data.frame(pfs=par('usr')[c(1,2)])
  yrng<-predict(av_list_mp[[i]]$lm,xrng)
  
  lines(x=unlist(xrng),y=yrng,col='red',lty=3,lwd=3)
  legend('topleft',
         legend = c(paste0('Correlation = ',round(av_list_mp[[i]]$c.cor,3)),
                    paste0('Correlation p val = ',formatC(av_list_mp[[i]]$p.cor,format = 'e', digits = 2),' (',av_list_mp[[i]]$s.cor,')'),
                    paste0('Adjusted p val = ',
                           formatC(p.adjust(av_list_mp[[i]]$p.cor,n = length(unique(av_DF$class))),format = 'e', digits = 2),
                           ' (',gtools::stars.pval(round(p.adjust(av_list_mp[[i]]$p.cor,n = length(unique(av_DF$class))),3)),')')),
         bty = 'n',
         cex=0.8,
         inset=c(-0.05,-0.1))
  dev.off()
}


##### SUB POP percentages ####
##### VS sample ####

av_DF<-lapply(setNames(unique(aggArea$uid),unique(aggArea$uid)),function(u){
  out<-lapply(setNames(unique(aggArea$super.cluster),unique(aggArea$super.cluster)),function(sc){
    subSet<- aggArea[aggArea$uid==u & aggArea$super.cluster==sc,]
    totArea<-sum(subSet$sub.area)
    out<-data.frame(pfs=subSet$pfs,
                    smp=subSet$sample,
                    super.cluster = subSet$super.cluster,
                    class= paste0(subSet$super.cluster,'_',subSet$sub.cluster),
                    area=subSet$sub.area/totArea*100)
  })
  out<-do.call(rbind,out)
})

av_DF<-do.call(rbind,av_DF)
av_DF[is.na(av_DF)]<-0

av_list_mp<-lapply(setNames(unique(av_DF$class),unique(av_DF$class)),function(x){
  av<-aov(area~smp,av_DF[av_DF$class==x,])
  
  av_HSD<-TukeyHSD(av)
  
  out<-data.frame(av_first_smp = unlist(lapply(strsplit(rownames(av_HSD$smp),'-'),'[',1)),
                  av_first_som = x,
                  av_second_smp = unlist(lapply(strsplit(rownames(av_HSD$smp),'-'),'[',2)),
                  av_second_som = x,
                  p.adj= av_HSD$smp[,"p adj"],
                  s.adj= gtools::stars.pval(av_HSD$smp[,"p adj"]),
                  row.names = NULL)
  out
})

av_list_mp<-do.call(rbind.data.frame,av_list_mp)


smpOrder<-unique(percentage_super[,c('sample','pfs')])
smpOrder<-factor(smpOrder$sample[order(as.numeric(smpOrder$pfs),decreasing = F)],
                 levels=smpOrder$sample[order(as.numeric(smpOrder$pfs),decreasing = F)],
                 ordered = T)
dir.create(file.path(analFolder,"SUB.POP"))

for (i in unique(av_DF$class)){
  subSet<-av_DF$area[av_DF$class==i]
  names(subSet)<-av_DF$smp[av_DF$class==i]
  smp<-factor(names(subSet),levels = levels(smpOrder),ordered = T)
  df_box<-data.frame(Sample=smp,
                     `% of parental area` = subSet,
                     check.names = F)
  postscript(file = file.path(analFolder,"SUB.POP",paste0(i,'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  
  par(mfrow=c(2,1),mar=c(0,5,2,1),las=1,bty='l')
  
  sbc<-strsplit(i,'_')[[1]]
  boxplot(`% of parental area`~Sample,
          data = df_box,
          main=paste0('Cluster ',sbc[2],' within ',prettyLabels[sbc[1]]),
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n')
  points(data.frame(jitter(as.integer(df_box$Sample)),df_box$`% of parental area`),bg='gray60',col='gray40',pch=21)
  boxplot(`% of parental area`~Sample,
          data = df_box,
          main=i,
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n',
          add=T)
  ll<-par('usr')
  plot(NA,
       xlim=c(ll[1],ll[2]),
       ylim=c(ll[2],ll[1]),
       bty='n',
       xlab=NA,
       ylab=NA,
       xaxt='n',
       yaxt='n',
       xaxs='i')
  axis(2,
       at=as.integer(smpOrder),
       labels = paste0("(",smpOrder,") ",LETTERS[1:length(smpOrder)][as.integer(smpOrder)]),
       cex.axis=0.6,ylab=NA,
       col=NA,
       col.ticks = 'black')
  axis(3,
       at=as.integer(smpOrder),
       labels = LETTERS[1:length(smpOrder)],
       cex.axis=0.6,
       col=NA,
       col.ticks = 'black')
  subSet<-av_list_mp[av_list_mp$av_first_som==i,]
  idx<-as.numeric(smpOrder)
  names(idx)<-levels(smpOrder)
  subSet[,"av_first_smp"]<-idx[subSet[,"av_first_smp"]]
  subSet[,"av_second_smp"]<-idx[subSet[,"av_second_smp"]]
  text(subSet[,"av_first_smp"],subSet[,"av_second_smp"],labels=subSet[,"s.adj"])
  for (iii in 1:nrow(subSet)){
    if (subSet$s.adj[iii]!=" " & subSet$s.adj[iii]!="."){
      ff<-subSet$av_first_smp[iii]
      ss<-subSet$av_second_smp[iii]
      
      jtr<-sample(seq(-0.3,0.3,0.01),1)
      sm<-c(ff+jtr,ll[3])
      pm<-c(ff+jtr,ss+jtr)
      em<-c(ll[1],ss+jtr)
      mm<-rbind(sm,pm,em)
      lines(mm,lty=3,col='gray60')
    }
  }
  dev.off()
}

##### VS pfs #####

av_DF<-lapply(setNames(unique(aggArea$uid),unique(aggArea$uid)),function(u){
  out<-lapply(setNames(unique(aggArea$super.cluster),unique(aggArea$super.cluster)),function(sc){
    subSet<- aggArea[aggArea$uid==u & aggArea$super.cluster==sc,]
    totArea<-as.numeric(sum(subSet$sub.area))
    out<-data.frame(pfs=as.numeric(subSet$pfs),
                    smp=subSet$sample,
                    super.cluster = subSet$super.cluster,
                    class= paste0(subSet$super.cluster,'_',subSet$sub.cluster),
                    area=as.numeric(subSet$sub.area)/totArea*100)
  })
  out<-do.call(rbind,out)
})

av_DF<-do.call(rbind,av_DF)
av_DF[is.na(av_DF)]<-0

av_list_mp<-lapply(setNames(unique(av_DF$class),unique(av_DF$class)),function(x){
  av_lm<-lm(area~pfs,av_DF[av_DF$class==x,])
  av_av<-anova(av_lm)
  av_cor<-cor.test(av_DF$pfs[av_DF$class==x],av_DF$area[av_DF$class==x])
  out<-list(av_som = x,
            p.anova= av_av[1,"Pr(>F)"],
            s.anova= gtools::stars.pval(av_av[1,"Pr(>F)"]),
            lm = av_lm,
            c.cor = av_cor$estimate,
            p.cor = av_cor$p.value,
            s.cor = gtools::stars.pval(av_cor$p.value))
  return(out)
})


dir.create(file.path(analFolder,"SUB.POP_pfs"))
for (i in unique(av_DF$class)){
  
  df_box<-data.frame(pfs=av_DF$pfs[av_DF$class==i],
                     smp=av_DF$smp[av_DF$class==i],
                     `% of total area` = av_DF$area[av_DF$class==i],
                     check.names = F)
  avr<-aggregate(`% of total area`~smp+pfs,df_box,mean)
  avr<-data.frame(smp = avr$smp,
                  pfs = as.numeric(avr$pfs),
                  mean = as.numeric(avr$`% of total area`))
  sdr<-aggregate(`% of total area`~smp+pfs,df_box,sd)
  sdr<-data.frame(smp = sdr$smp,
                  pfs = as.numeric(sdr$pfs),
                  sd = as.numeric(sdr$`% of total area`))
  avrsdr<- data.frame(smp=avr$smp,
                      pfs=avr$pfs,
                      mean=avr$mean,
                      sd=sdr$sd,
                      `mean+sd`=avr$mean+sdr$sd,
                      `mean-sd`=avr$mean-sdr$sd,
                      check.names = F)
  
  postscript(file = file.path(analFolder,"SUB.POP_pfs",paste0(i,'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  par(xpd=T)
  sbc<-strsplit(i,'_')[[1]]
  plot(data.frame(jitter(as.integer(df_box$pfs)),
                  df_box$`% of total area`),
       bg='gray60',
       col='gray40',
       pch=21,
       ylab='% of parental area',
       xlab='pfs',
       
       bty='l',
       las=1,
       ylim=c(min(c(avrsdr$`mean-sd`,df_box$`% of total area`)),max(avrsdr$`mean+sd`,df_box$`% of total area`)))
  mtext(side = 3,adj=1,paste0('Cluster ',sbc[2],' within ',prettyLabels[sbc[1]]))
  for (ii in 1:nrow(avrsdr)){
    lines(matrix(c(avrsdr$pfs[ii]-0.5,avrsdr$pfs[ii]+0.5,avrsdr$mean[ii],avrsdr$mean[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(avrsdr$pfs[ii],avrsdr$pfs[ii],avrsdr$`mean+sd`[ii],avrsdr$`mean-sd`[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(avrsdr$pfs[ii]-0.2,avrsdr$pfs[ii]+0.2,avrsdr$`mean+sd`[ii],avrsdr$`mean+sd`[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
    lines(matrix(c(avrsdr$pfs[ii]-0.2,avrsdr$pfs[ii]+0.2,avrsdr$`mean-sd`[ii],avrsdr$`mean-sd`[ii]),
                 ncol = 2,
                 byrow = F),
          col='black')
  }
  
  xrng<-data.frame(pfs=par('usr')[c(1,2)])
  yrng<-predict(av_list_mp[[i]]$lm,xrng)
  
  lines(x=unlist(xrng),y=yrng,col='red',lty=3,lwd=3)
  legend('topleft',
         legend = c(paste0('Correlation = ',round(av_list_mp[[i]]$c.cor,3)),
                    paste0('Correlation p val = ',formatC(av_list_mp[[i]]$p.cor,format = 'e', digits = 2),' (',av_list_mp[[i]]$s.cor,')'),
                    paste0('Adjusted p val = ',
                           formatC(p.adjust(av_list_mp[[i]]$p.cor,n = length(unique(av_DF$class))),format = 'e', digits = 2),
                           ' (',gtools::stars.pval(round(p.adjust(av_list_mp[[i]]$p.cor,n = length(unique(av_DF$class))),3)),')')),
         bty = 'n',
         cex=0.8,
         inset=c(-0.05,-0.1))
  dev.off()
}


##### INTERACTIONS ######

superGroup_factor<-factor(sort(unique(expGate$primary_id[!expGate$dump])))
subgrGroup_factor<-factor(sort(unique(expGate$uMap_som[!expGate$dump])))
subgrGroupCollapse_factor<-factor(sort(apply(expand.grid(superGroup_factor,subgrGroup_factor),1,function(x)paste(x[1],x[2],sep='_'))))
intrctn<-pbapply::pblapply(setNames(studyTable$uid,studyTable$uid),function(u){
  
  subSet_data<-expData[expGate$uid==u & !expGate$dump,'geom']
  subSet_gate<-expGate[expGate$uid==u & !expGate$dump,]
  intrctn<-st_intersects(subSet_data,sparse = F)
  maskMatrix<-lower.tri(intrctn,diag=F)
  intrctn<-intrctn & maskMatrix
  intrctnComposition<-which(intrctn,arr.ind = T)
  tb<-table(unlist(intrctnComposition))
  tb.stats<-data.frame(uid=u,
                       sample=studyTable$sample[studyTable$uid==u],
                       median_interaction=median(tb))
  
  class_list<-lapply(1:nrow(intrctn),function(i){
    rowi<-which(intrctn[i,])
    coli<-which(intrctn[,i])
    listi<-c(rowi,coli)
    if (length(listi)!=0) {
      out<-list(o.super = subSet_gate$primary_id[i],
                t.super = subSet_gate$primary_id[listi],
                o.sub = factor(paste0(subSet_gate$primary_id[i],'_',subSet_gate$uMap_som[i]),levels=levels(subgrGroupCollapse_factor)),
                t.sub = factor(paste0(subSet_gate$primary_id[listi],'_',subSet_gate$uMap_som[listi]),levels=levels(subgrGroupCollapse_factor)))
    } else {
      out<-list(o.super = subSet_gate$primary_id[i],
                t.super = factor(0,levels=levels(superGroup_factor)),
                o.sub = factor(paste0(subSet_gate$primary_id[i],'_',subSet_gate$uMap_som[i]),levels=levels(subgrGroupCollapse_factor)),
                t.sub = factor(0,levels=levels(subgrGroupCollapse_factor)))
    }
    return(out)
  })
  
  listOut<-lapply(class_list,function(x) !is.null(x$o.super) | !is.null(x$o.sub) )
  class_list<-class_list[unlist(listOut)]
  
  extra.sub<-unique(unlist(lapply(class_list,'[[',"o.sub")))
  
  extra.sub<-levels(subgrGroupCollapse_factor)[!(levels(subgrGroupCollapse_factor) %in% extra.sub)]
  
  if (length(extra.sub)!=0){
    extra.sub<-lapply(extra.sub,function(o.sub){
      newLabels<-strsplit(o.sub,'_')[[1]][1]
      out<-list(o.super=factor(newLabels,levels=levels(superGroup_factor)),
                t.super=factor(0,levels=levels(superGroup_factor)),
                o.sub = factor(o.sub,levels = levels(subgrGroupCollapse_factor)),
                t.sub = factor(0,levels=levels(subgrGroupCollapse_factor)))
    })
    
    class_list<-append(class_list,extra.sub)
  }
  
  class_table<-lapply(class_list,function(x){
    out<-list(o.super = x$o.super,
              t.super = table(x$t.super),
              o.sub = x$o.sub,
              t.sub = table(x$t.sub))
  })
  
  # extra.super<-unique(unlist(lapply(class_table,'[[',"o.super")))
  # 
  # extra.super<-levels(superGroup_factor)[!(levels(superGroup_factor) %in% extra.super)]
  # 
  # if (length(extra.super)!=0){
  #   extra.super<-lapply(extra.super,function(o.super){
  #     newLabels<-paste(o.super,subgrGroup_factor,sep='_')
  #     out<-lapply(newLabels,function(o.sub){
  #       out<-list(o.super=factor(o.super,levels=levels(superGroup_factor)),
  #                 t.super=as.table(setNames(rep(0,length(levels(superGroup_factor))),superGroup_factor)),
  #                 o.sub = factor(o.sub,levels = levels(subgrGroupCollapse_factor)),
  #                 t.sub = as.table(setNames(rep(0,length(levels(subgrGroupCollapse_factor))),subgrGroupCollapse_factor)))
  #     })
  #   })
  #   extra.super<-unlist(extra.super,recursive = F)
  #   
  #   class_table<-append(class_table,extra.super)
  # }
  # 
  # extra.sub<-unique(unlist(lapply(class_table,'[[',"o.sub")))
  # 
  # extra.sub<-levels(subgrGroupCollapse_factor)[!(levels(subgrGroupCollapse_factor) %in% extra.sub)]
  # 
  # if (length(extra.sub)!=0){
  #   extra.sub<-lapply(extra.sub,function(o.sub){
  #     newLabels<-strsplit(o.sub,'_')[[1]][1]
  #     out<-list(o.super=factor(newLabels,levels=levels(superGroup_factor)),
  #               t.super=as.table(setNames(rep(0,length(levels(superGroup_factor))),superGroup_factor)),
  #               o.sub = factor(o.sub,levels = levels(subgrGroupCollapse_factor)),
  #               t.sub = as.table(setNames(rep(0,length(levels(subgrGroupCollapse_factor))),subgrGroupCollapse_factor)))
  #     })
  #   browser()
  #   # extra.sub<-unlist(extra.super,recursive = F)
  #   class_table<-append(class_table,extra.sub)
  # }
  # 
  
  superClass_abundance<-do.call(rbind.data.frame,lapply(class_table,function(x){
    out<-data.frame(uid=u,
                    sample=studyTable$sample[studyTable$uid==u],
                    origin=x$o.super)
    out<-cbind.data.frame(out,x$t.super)
    test<-names(out)
    names(out)[c(4,5)]<-c('target','N')
    return(out)
  }))
  
  # extraClass<-superGroup_factor[!(superGroup_factor %in% superClass_abundance$origin)]
  # 
  # if (length(extraClass)!=0){
  #   extraDF<-do.call(rbind.data.frame,lapply(extraClass,function(ec){
  #     out<-data.frame(uid=u,
  #                     sample=studyTable$sample[studyTable$uid==u],
  #                     origin=ec,
  #                     target=superGroup_factor,
  #                     N=0)
  #     return(out)
  #   }))
  #   
  #   superClass_abundance<-rbind.data.frame(superClass_abundance,extraDF)
  #   superClass_abundance<-aggregate(N~uid+sample+origin+target,superClass_abundance,mean)
  #   superClass_abundance<-superClass_abundance[order(superClass_abundance$origin),]
  # }
  
  subClass_abundance<-do.call(rbind.data.frame,lapply(class_table,function(x){
    out<-data.frame(uid=u,
                    sample=studyTable$sample[studyTable$uid==u],
                    origin=x$o.sub)
    out<-cbind.data.frame(out,x$t.sub)
    names(out)[c(4,5)]<-c('target','N')
    return(out)
  }))
  
  
  subClass_abundance<-aggregate(N~uid+sample+origin+target,subClass_abundance,mean)
  subClass_abundance<-subClass_abundance[order(superClass_abundance$origin),]
  
  class_entropy<-lapply(class_table,function(x){
    p.super<-x$t.super[x$t.super!=0]/sum(x$t.super[x$t.super!=0])
    l.super<-length(p.super)
    p.sub<-x$t.sub[x$t.sub!=0]/sum(x$t.sub[x$t.sub!=0])
    l.sub<-length(p.sub)
    
    out<-list(o.super = x$o.super,
              e.super = if(l.super==1) 0 else -sum(p.super*log(p.super,length(p.super))),
              o.sub = x$o.sub,
              e.sub = if(l.sub==1) 0 else -sum(p.sub*log(p.sub,length(p.sub))))
  })
  
  superClass_entropy<-do.call(rbind.data.frame,lapply(class_entropy,function(x){
    out<-data.frame(uid=u,
                    sample=studyTable$sample[studyTable$uid==u],
                    origin=x$o.super,
                    entropy=x$e.super)
    return(out)
  }))
  superClass_entropy<-aggregate(entropy~uid+sample+origin,superClass_entropy,mean)
  superClass_entropy<-superClass_entropy[order(superClass_entropy$origin),]
  
  subClass_entropy<-do.call(rbind.data.frame,lapply(class_entropy,function(x){
    out<-data.frame(uid=u,
                    sample=studyTable$sample[studyTable$uid==u],
                    origin=x$o.sub,
                    entropy=x$e.sub)
    return(out)
  }))
  subClass_entropy<-aggregate(entropy~uid+sample+origin,subClass_entropy,mean)
  subClass_entropy<-subClass_entropy[order(subClass_entropy$origin),]
  
  out<-list(median_interaction = tb.stats,
            abundance.super = superClass_abundance,
            entropy.super = superClass_entropy,
            abundance.sub = subClass_abundance,
            entropy.sub = subClass_entropy)
  return(out)
})

saveRDS(intrctn,file.path(analFolder,'bigTables','interactions.R'))



##### plot super.entropy #####
TEMP<-do.call(rbind.data.frame,lapply(intrctn,'[[','entropy.super'))
TEMP_correct<-TEMP
TEMP_presence<-table(expGate[!expGate$dump,c('primary_id','uid')])
for (uid in colnames(TEMP_presence)){
  for (prid in rownames(TEMP_presence)){
    if (TEMP_presence[prid,uid]==0) TEMP_correct$entropy[TEMP_correct$uid==uid & TEMP_correct$origin==prid]<-NA
  }
}
av_list_mp<-lapply(setNames(levels(TEMP_correct$origin),levels(TEMP_correct$origin)),function(x){
  av<-aov(entropy~sample,TEMP_correct[TEMP_correct$origin==x,])
  
  av_HSD<-TukeyHSD(av)
  
  out<-data.frame(origin = x,
                  sample1 = unlist(lapply(strsplit(rownames(av_HSD$sample),'-'),'[',1)),
                  sample2 = unlist(lapply(strsplit(rownames(av_HSD$sample),'-'),'[',2)),
                  p.adj= av_HSD$sample[,"p adj"],
                  s.adj= gtools::stars.pval(av_HSD$sample[,"p adj"]),
                  row.names = NULL)
  return(out)
})

av_list_mp<-do.call(rbind.data.frame,av_list_mp)

smpOrder<-unique(studyTable[,c('sample','bioGroup')])
smpOrder<-factor(smpOrder$sample[order(as.numeric(smpOrder$bioGroup),decreasing = F)],
                 levels=smpOrder$sample[order(as.numeric(smpOrder$bioGroup),decreasing = F)],
                 ordered = T)

dir.create(file.path(analFolder,"SUPER.ENTROPY"))
for (i in levels(TEMP_correct$origin)){
  subSet<-TEMP_correct$entropy[TEMP_correct$origin==i]
  names(subSet)<-TEMP_correct$sample[TEMP_correct$origin==i]
  subSet_untouched<-TEMP$entropy[TEMP$origin==i]
  names(subSet_untouched)<-TEMP$sample[TEMP$origin==i]
  smp<-factor(names(subSet),levels = levels(smpOrder),ordered = T)
  df_box<-data.frame(Sample=smp,
                     Entropy = subSet,
                     check.names = F)
  df_points<-data.frame(sample=smp,
                        Entropy = subSet_untouched,
                        check.names = F)
  df_symbols<-data.frame(samples=smp,
                         pch = unlist(lapply(subSet,function(x)if (is.na(x)) return(4) else return(21))),
                         check.names = F)
  postscript(file = file.path(analFolder,"SUPER.ENTROPY",paste0(i,'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  
  par(mfrow=c(2,1),mar=c(0,5,3,1),las=1,bty='l',xpd=T)
  
  boxplot(`Entropy`~Sample,
          data = df_box,
          # main=prettyLabels[i],
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n',
          ylim=c(-0.1,1))
  title(prettyLabels[i],line=2)
  points(data.frame(jitter(as.integer(df_points$sample)),df_points$Entropy),
         bg='gray60',
         col='gray40',
         pch=df_symbols$pch)
  boxplot(Entropy~Sample,
          data = df_box,
          main=i,
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n',
          add=T)
  
  legend('top',
         legend = c('valid interaction','not present'),
         pch = c(21,4),
         horiz = T,
         bty='n',
         pt.bg='gray60',
         col='gray40',
         inset=c(0,-0.1),
         cex = 0.8)
  
  ll<-par('usr')
  plot(NA,
       xlim=c(ll[1],ll[2]),
       ylim=c(ll[2],ll[1]),
       bty='n',
       xlab=NA,
       ylab=NA,
       xaxt='n',
       yaxt='n',
       xaxs='i')
  axis(2,
       at=as.integer(smpOrder),
       labels = paste0("(",smpOrder,") ",LETTERS[1:length(smpOrder)][as.integer(smpOrder)]),
       cex.axis=0.6,ylab=NA,
       col=NA,
       col.ticks = 'black')
  axis(3,
       at=as.integer(smpOrder),
       labels = LETTERS[1:length(smpOrder)],
       cex.axis=0.6,
       col=NA,
       col.ticks = 'black')
  subSet<-av_list_mp[av_list_mp$origin==i,]
  idx<-as.numeric(smpOrder)
  names(idx)<-levels(smpOrder)
  subSet[,"sample1"]<-idx[subSet[,"sample1"]]
  subSet[,"sample2"]<-idx[subSet[,"sample2"]]
  text(subSet[,"sample1"],subSet[,"sample2"],labels=subSet[,"s.adj"])
  for (iii in 1:nrow(subSet)){
    if (subSet$s.adj[iii]!=" " & subSet$s.adj[iii]!="."){
      ff<-subSet$sample1[iii]
      ss<-subSet$sample2[iii]
      
      jtr<-sample(seq(-0.3,0.3,0.01),1)
      sm<-c(ff+jtr,ll[3])
      pm<-c(ff+jtr,ss+jtr)
      em<-c(ll[1],ss+jtr)
      mm<-rbind(sm,pm,em)
      lines(mm,lty=3,col='gray60')
    }
  }
  dev.off()
}

#####plot sub.entropy ####

TEMP<-do.call(rbind.data.frame,lapply(intrctn,'[[','entropy.sub'))
TEMP_correct<-TEMP
TEMP_TEMP<-expGate[!expGate$dump,c('primary_id','uMap_som','uid')]
TEMP_TEMP<-data.frame(id=paste(TEMP_TEMP$primary_id,TEMP_TEMP$uMap_som,sep='_'),
                      uid=TEMP_TEMP$uid)
TEMP_presence<-table(TEMP_TEMP)
rm(TEMP_TEMP)
for (uid in colnames(TEMP_presence)){
  for (prid in rownames(TEMP_presence)){
    if (TEMP_presence[prid,uid]==0) TEMP_correct$entropy[TEMP_correct$uid==uid & TEMP_correct$origin==prid]<-NA
  }
}
av_list_mp<-lapply(setNames(unique(TEMP_correct$origin),unique(TEMP_correct$origin)),function(x){
  av<-aov(entropy~sample,TEMP_correct[TEMP_correct$origin==x,])
  
  av_HSD<-TukeyHSD(av)
  out<-data.frame(origin = x,
                  sample1 = unlist(lapply(strsplit(rownames(av_HSD$sample),'-'),'[',1)),
                  sample2 = unlist(lapply(strsplit(rownames(av_HSD$sample),'-'),'[',2)),
                  p.adj= av_HSD$sample[,"p adj"],
                  s.adj= gtools::stars.pval(av_HSD$sample[,"p adj"]),
                  row.names = NULL)
  return(out)
})

av_list_mp<-do.call(rbind.data.frame,av_list_mp)


smpOrder<-unique(percentage_super[,c('sample','pfs')])
smpOrder<-factor(smpOrder$sample[order(as.numeric(smpOrder$pfs),decreasing = F)],
                 levels=smpOrder$sample[order(as.numeric(smpOrder$pfs),decreasing = F)],
                 ordered = T)
dir.create(file.path(analFolder,"SUB.ENTROPY"))

for (i in unique(TEMP_correct$origin)){
  subSet<-TEMP_correct$entropy[TEMP_correct$origin==i]
  names(subSet)<-TEMP_correct$sample[TEMP_correct$origin==i]
  subSet_untouched<-TEMP$entropy[TEMP$origin==i]
  names(subSet_untouched)<-TEMP$sample[TEMP$origin==i]
  smp<-factor(names(subSet),levels = levels(smpOrder),ordered = T)
  df_box<-data.frame(Sample=smp,
                     Entropy = subSet,
                     check.names = F)
  df_points<-data.frame(sample=smp,
                        Entropy = subSet_untouched,
                        check.names = F)
  df_symbols<-data.frame(samples=smp,
                         pch = unlist(lapply(subSet,function(x)if (is.na(x)) return(4) else return(21))),
                         check.names = F)
  postscript(file = file.path(analFolder,"SUB.ENTROPY",paste0(i,'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  
  par(mfrow=c(2,1),mar=c(0,5,3,1),las=1,bty='l',xpd=T)
  
  boxplot(`Entropy`~Sample,
          data = df_box,
          # main=prettyLabels[i],
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n',
          ylim=c(-0.1,1))
  TEMP_title<-strsplit(as.character(i),'_')[[1]]
  title(paste0('Cluster ',TEMP_title[2],' within ',prettyLabels[TEMP_title[1]]),line=2)
  points(data.frame(jitter(as.integer(df_points$sample)),df_points$Entropy),
         bg='gray60',
         col='gray40',
         pch=df_symbols$pch)
  boxplot(Entropy~Sample,
          data = df_box,
          main=i,
          col=NA,
          pch=NA,
          bty='l',
          xlab=NA,
          xaxt='n',
          add=T)
  
  legend('top',
         legend = c('valid interaction','not present'),
         pch = c(21,4),
         horiz = T,
         bty='n',
         pt.bg='gray60',
         col='gray40',
         inset=c(0,-0.1),
         cex = 0.8)
  
  ll<-par('usr')
  plot(NA,
       xlim=c(ll[1],ll[2]),
       ylim=c(ll[2],ll[1]),
       bty='n',
       xlab=NA,
       ylab=NA,
       xaxt='n',
       yaxt='n',
       xaxs='i')
  axis(2,
       at=as.integer(smpOrder),
       labels = paste0("(",smpOrder,") ",LETTERS[1:length(smpOrder)][as.integer(smpOrder)]),
       cex.axis=0.6,ylab=NA,
       col=NA,
       col.ticks = 'black')
  axis(3,
       at=as.integer(smpOrder),
       labels = LETTERS[1:length(smpOrder)],
       cex.axis=0.6,
       col=NA,
       col.ticks = 'black')
  subSet<-av_list_mp[av_list_mp$origin==i,]
  idx<-as.numeric(smpOrder)
  names(idx)<-levels(smpOrder)
  subSet[,"sample1"]<-idx[subSet[,"sample1"]]
  subSet[,"sample2"]<-idx[subSet[,"sample2"]]
  text(subSet[,"sample1"],subSet[,"sample2"],labels=subSet[,"s.adj"])
  for (iii in 1:nrow(subSet)){
    if (subSet$s.adj[iii]!=" " & subSet$s.adj[iii]!="."){
      ff<-subSet$sample1[iii]
      ss<-subSet$sample2[iii]
      
      jtr<-sample(seq(-0.3,0.3,0.01),1)
      sm<-c(ff+jtr,ll[3])
      pm<-c(ff+jtr,ss+jtr)
      em<-c(ll[1],ss+jtr)
      mm<-rbind(sm,pm,em)
      lines(mm,lty=3,col='gray60')
    }
  }
  dev.off()
}

###### plot abundance super#######

TEMP<-do.call(rbind.data.frame,lapply(intrctn,'[[','abundance.super'))

TEMP.mean<-aggregate(N~sample+origin+target,TEMP,mean)
TEMP.sd<-aggregate(N~sample+origin+target,TEMP,sd)
smpOrder<-unique(studyTable[,c('sample','bioGroup')])
smpOrder<-factor(smpOrder$sample[order(as.numeric(smpOrder$bioGroup),decreasing = F)],
                 levels=smpOrder$sample[order(as.numeric(smpOrder$bioGroup),decreasing = F)],
                 ordered = T)

TEMP.big<-matrix(NA,
                 nrow=length(unique(TEMP$origin))*sqrt(length(levels(smpOrder))),
                 ncol=length(unique(TEMP$target))*sqrt(length(levels(smpOrder))),
                 dimnames = list(origin=rep(levels(superGroup_factor),3),
                                 target=rep(levels(superGroup_factor),3)))
ii=1
for(rr in 0:2){
  for(cc in 0:2){
    TEMP.sub<-TEMP.mean[TEMP.mean$sample==levels(smpOrder)[ii],]
    TEMP.sub<-lapply(unique(TEMP.sub$origin),function(oo){
      matrix(TEMP.sub[TEMP.sub$origin==oo,'N'],
             nrow=1,
             dimnames = list(origin=oo,target=TEMP.sub$target[TEMP.sub$origin==oo]))
    })
    TEMP.sub<-do.call(rbind,TEMP.sub)
    # TEMP.sub<-matrix(TEMP.sub,ncol=ncol(TEMP.sub))
    TEMP.big[(1:nrow(TEMP.sub))+nrow(TEMP.sub)*rr,
             (1:ncol(TEMP.sub))+ncol(TEMP.sub)*cc]<-TEMP.sub
    ii=ii+1
  }}

# colnames(TEMP.big)<-make.names(prettyLabels[colnames(TEMP.big)],unique = T)
# rownames(TEMP.big)<-make.names(prettyLabels[rownames(TEMP.big)],unique = T)

 rn<-rownames(TEMP.big)
 cn<-colnames(TEMP.big)
 rownames(TEMP.big)<-make.names(rn,unique = T,allow_ = T)
 colnames(TEMP.big)<-make.names(cn,unique = T,allow_ = T)

 grPalette<- setNames(RColorBrewer::brewer.pal(length(unique(prettyLabels)),name = 'Set2'),prettyLabels)
grAnnotationR<-data.frame(`Acceptor type`=factor(prettyLabels[rn],levels=prettyLabels[colnames(TEMP.sub)],ordered = T),
                         check.names = F,
                         row.names = rownames(TEMP.big),
                         check.rows = F)
grAnnotationC<-data.frame(`Donor type`=factor(prettyLabels[rn],levels=prettyLabels[colnames(TEMP.sub)],ordered = T),
                          check.names = F,
                          row.names = rownames(TEMP.big),
                          check.rows = F)
grAnnotation_color<-list(`Donor type`=grPalette[levels(grAnnotation$`Cell types`)],
                         `Acceptor type`=grPalette[levels(grAnnotation$`Cell types`)])

pheatmap::pheatmap(TEMP.big,
                   gaps_row = c(1:2)*nrow(TEMP.sub),
                   gaps_col = c(1:2)*ncol(TEMP.sub),
                   cellwidth = 15,
                   cellheight =15,
                   height = 6,
                   width = 7,
                   scale = 'none',
                   # annotation_row = grAnnotation,
                   filename = file.path(analFolder,'SUPER.ENTROPY','NOI.pdf'),
                   annotation_row = grAnnotationR,
                   annotation_col = grAnnotationC,
                   annotation_colors = grAnnotation_color,
                   show_rownames = F,
                   show_colnames = F,
                   cluster_rows = F,
                   cluster_cols = F)
# 
###### plot abundance sub#######

TEMP<-do.call(rbind.data.frame,lapply(intrctn,'[[','abundance.sub'))

TEMP.mean<-aggregate(N~sample+origin+target,TEMP,mean)
TEMP.sd<-aggregate(N~sample+origin+target,TEMP,sd)
smpOrder<-unique(studyTable[,c('sample','bioGroup')])
smpOrder<-factor(smpOrder$sample[order(as.numeric(smpOrder$bioGroup),decreasing = F)],
                 levels=smpOrder$sample[order(as.numeric(smpOrder$bioGroup),decreasing = F)],
                 ordered = T)

TEMP.big<-matrix(NA,
                 nrow=length(unique(TEMP.mean$origin))*sqrt(length(levels(smpOrder))),
                 ncol=length(unique(TEMP.mean$target))*sqrt(length(levels(smpOrder))),
                 dimnames = list(origin=rep(levels(subgrGroupCollapse_factor),3),
                                 target=rep(levels(subgrGroupCollapse_factor),3)))
ii=1
for(rr in 0:2){
  for(cc in 0:2){
    TEMP.sub<-TEMP.mean[TEMP.mean$sample==levels(smpOrder)[ii],]
    TEMP.sub<-lapply(unique(TEMP.sub$origin),function(oo){
      matrix(TEMP.sub[TEMP.sub$origin==oo,'N'],
             nrow=1,
             dimnames = list(origin=oo,target=TEMP.sub$target[TEMP.sub$origin==oo]))
    })
    TEMP.sub<-do.call(rbind,TEMP.sub)
    # TEMP.sub<-matrix(TEMP.sub,ncol=ncol(TEMP.sub))
    TEMP.big[(1:nrow(TEMP.sub))+nrow(TEMP.sub)*rr,
             (1:ncol(TEMP.sub))+ncol(TEMP.sub)*cc]<-TEMP.sub
    ii=ii+1
  }}

rn<-rownames(TEMP.big)
cn<-colnames(TEMP.big)
rownames(TEMP.big)<-make.names(rn,unique = T,allow_ = T)
colnames(TEMP.big)<-make.names(cn,unique = T,allow_ = T)

rn.super<-unlist(lapply(strsplit(rn,'_'),'[',1))
rn.sub<-unlist(lapply(strsplit(rn,'_'),'[',2))
grPalette.super<- setNames(RColorBrewer::brewer.pal(length(unique(prettyLabels)),name = 'Set2'),prettyLabels)
grPalette.sub<-setNames(RColorBrewer::brewer.pal(length(unique(expGate$uMap_som)),name = 'Set3'),levels(expGate$uMap_som))[1:8]
grAnnotationC<-data.frame(`Donor types`=factor(prettyLabels[rn.super],levels=prettyLabels[levels(expGate$primary_id)],ordered = T),
                         `Donor Clusters` = factor(rn.sub,levels = levels(expGate$uMap_som), ordered = T),
                         check.names = F,
                         row.names = rownames(TEMP.big),
                         check.rows = F)
grAnnotationR<-data.frame(`Acceptor types`=factor(prettyLabels[rn.super],levels=prettyLabels[levels(expGate$primary_id)],ordered = T),
                          `Acceptor Clusters` = factor(rn.sub,levels = levels(expGate$uMap_som), ordered = T),
                          check.names = F,
                          row.names = rownames(TEMP.big),
                          check.rows = F)
grAnnotation_color<-list(`Donor types`=grPalette.super[levels(grAnnotation$`Cell types`)],
                         `Donor Clusters` = grPalette.sub[levels(grAnnotation$Clusters)],
                         `Acceptor types`=grPalette.super[levels(grAnnotation$`Cell types`)],
                         `Acceptor Clusters` = grPalette.sub[levels(grAnnotation$Clusters)])


pheatmap::pheatmap(TEMP.big,
                   gaps_row = c(1:2)*nrow(TEMP.sub),
                   gaps_col = c(1:2)*ncol(TEMP.sub),
                   cellwidth = 8,
                   cellheight =8,
                   height = 21,
                   width = 22,
                   scale = 'none',
                   # annotation_row = grAnnotation,
                   filename = file.path(analFolder,'SUB.ENTROPY','NOI.pdf'),
                   annotation_row = grAnnotationR,
                   annotation_col = grAnnotationC,
                   annotation_colors = grAnnotation_color,
                   show_rownames = F,
                   show_colnames = F,
                   cluster_rows = F,
                   cluster_cols = F)

##### MAPS ######## 
dir.create(file.path(analFolder,'MAPS'))

for (u in unique(expGate$uid)){
  smp<-unique(expGate$sample[expGate$uid==u])
  uid<-u
  subSet_gate<-expGate[expGate$uid==u,]
  subSet_data<-expData$geom[expData$uid==u]
  postscript(file = file.path(analFolder,"MAPS",paste0(smp,'_',u,'.eps')),
             onefile = F,
             width = 20,
             height = 20,
             horizontal = F,
             paper = 'special',
             bg='white')
  plot.new()
  par(mfrow=c(3,3),mar=c(1,1,2,1))
  for (i in primary_id){
    lbl<-prettyLabels[i]
    ss<-subSet_data[!subSet_gate$dump & subSet_gate$primary_id==i]
    ss_color<-grPalette.sub[subSet_gate$uMap_som[!subSet_gate$dump & subSet_gate$primary_id==i]]
    bbcoord<-st_bbox(subSet_data)
    plot(NA,
         xlim=c(bbcoord[c(1,3)]),
         ylim=c(bbcoord[c(2,4)]),
         xaxt='n',
         yaxt='n',
         bty='n',
         xaxs='i',
         yaxs='i',asp=1)
    rect(bbcoord[1],bbcoord[2],bbcoord[3],bbcoord[4],col='black')
    plot(subSet_data,
         col=NA,
         border='gray50',
         lwd=0.1,
         key.pos = NULL,
         reset = F,
         add=T)
    title( main=prettyLabels[i],cex.main=3)
    plot(ss,col=ss_color,
         border='gray20',
         lwd=0.1,
         add=T,
         key.pos=NULL,
         reset = F)
  }
  ss<-subSet_data[subSet_gate$dump & subSet_gate$primary_id==i]
  plot(NA,
       xlim=c(bbcoord[c(1,3)]),
       ylim=c(bbcoord[c(2,4)]),
       xaxt='n',
       yaxt='n',
       bty='n',
       xaxs='i',
       yaxs='i',asp=1)
  rect(bbcoord[1],bbcoord[2],bbcoord[3],bbcoord[4],col='black')
  plot(subSet_data,
       col=NA,
       border='gray50',
       lwd=0.1,
       key.pos = NULL,
       reset = F,
       add=T)
  title( main='Rejected cells',cex.main=3)
  plot(ss,col='red',
       border='gray20',
       lwd=0.1,
       add=T,
       key.pos=NULL,
       reset = F)
  plot(NA,
       xlim=c(0,12),
       ylim=c(0,12),
       bty='n',
       ann = F,
       xaxt='n',
       yaxt='n')
  title(main=smp,cex.main=3)
  text(1,12,paste0('ROI: ',u),cex=2,adj=c(0,0))
  for(ii in 1:length(grPalette.sub)){
    rect(1,10-ii,2,10.9-ii,col = grPalette.sub[ii],border = 'gray20')
    text(3,10.5-ii,names(grPalette.sub)[ii],cex=2.5,adj=c(0,0.5))
  }
  rect(1,10-ii-1,2,10.9-ii-1,col = 'red',border = 'gray20')
  text(3,10.5-ii-1,'Rejected cells',cex=2.5,adj=c(0,0.5))
  dev.off()
}

# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# totComposition_super<-paste(expGate$primary_id[expGate$uid==u & !expGate$dump][intrctnComposition[,1]],
#                             expGate$primary_id[expGate$uid==u & !expGate$dump][intrctnComposition[,2]],
#                             sep='_x_')
# totComposition_super<-table(totComposition_super)
# 
# 
# totComposition_sub<-paste(paste(expGate$primary_id[expGate$uid==u & !expGate$dump][intrctnComposition[,1]],
#                                 expGate$uMap_som[expGate$uid==u & !expGate$dump][intrctnComposition[,1]],sep='_'),
#                           paste(expGate$primary_id[expGate$uid==u & !expGate$dump][intrctnComposition[,2]],
#                                 expGate$uMap_som[expGate$uid==u & !expGate$dump][intrctnComposition[,2]],sep='_'),sep='_x_')
# totComposition_sub<-table(totComposition_sub)
# 
# 
# norm_sub<-apply(subClusters_T_intrct,1,function(x){
#   if (x[1] %in% names(totComposition_sub)) {
#     actualInteraction<-as.numeric(totComposition_sub[x[1]])
#   } else actualInteraction<-0
#   out<-setNames(c(x[1],actualInteraction/as.numeric(x[2])),c('V','N'))
#   if (is.na(out[2])) out[2]<-0
#   return(out)
# })
# norm_sub<-as.data.frame(t(norm_sub))
# norm_sub$N<-as.numeric(norm_sub$N)
# 
# norm_super<-apply(superClusters_T_intrct,1,function(x){
#   if (x[1] %in% names(totComposition_super)) {
#     actualInteraction<-as.numeric(totComposition_super[x[1]])
#   } else actualInteraction<-0
#   out<-setNames(c(x[1],actualInteraction/as.numeric(x[2])),c('V','N'))
#   if (is.na(out[2])) out[2]<-0
#   return(out)
# })
# norm_super<-as.data.frame(t(norm_super))
# norm_super$N<-as.numeric(norm_super$N)
# 
# browser()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# out<-lapply(intrctn,function(x){
#   if (length(x)>1){
#     out<-data.frame(uid=u,
#                     sample=smp,
#                     id.origin = pIDs[x[1]],
#                     x.origin = subSet_ump[x[1],1],
#                     y.origin = subSet_ump[x[1],2],
#                     z.origin = primaryIDs[x[1]],
#                     id.target = pIDs[x[-1]],
#                     x.target = subSet_ump[x[-1],1],
#                     y.target = subSet_ump[x[-1],2],
#                     z.target = primaryIDs[x[-1]],
#                     row.names = NULL,
#                     stringsAsFactors = F)}
#   else {
#     out<-data.frame(uid=character(0),
#                     sample=character(0),
#                     id.origin = character(0),
#                     x.origin = numeric(0),
#                     y.origin = numeric(0),
#                     z.origin = character(0),
#                     id.target = character(0),
#                     x.target = numeric(0),
#                     y.target = numeric(0),
#                     z.target = character(0),
#                     row.names = NULL,
#                     stringsAsFactors = F)
#   }
#   return(out)
# })
# out<-do.call(rbind,out)
# })
# intrctn<-do.call(rbind,intrctn)
# 
# dir.create(file.path(analFolder,"interactions"))
# 
# combos<-expand.grid(primary_id,primary_id)
# combos_out<-apply(combos,1,function(x)x[1]==x[2])
# combos<-combos[!combos_out,]
# combosPalette<-rainbow(nrow(combos))
# for (smp in unique(expGate$sample)){
#   
#   # postscript(file = file.path(analFolder,"interactions",paste0(smp,'.eps')),
#   #            onefile = F,
#   #            width = 7.5,
#   #            height = 7.5,
#   #            horizontal = F,
#   #            paper = 'special',
#   #            bg='white')
#   pdf(file = file.path(analFolder,"interactions",paste0(smp,'.pdf')),
#       onefile = F,
#       width = 7.5,
#       height = 7.5,
#       paper = 'special',
#       bg='white')
#   
#   
#   plot(expGate[expGate$sample==smp & !expGate$dump,c('uMAP1','uMAP2')],
#        col='gray80',
#        cex=0.3,
#        main = paste0(smp,':',unique(studyTable$bioGroup[studyTable$sample==smp])))
#   
#   for (cc in 1:nrow(combos)){
#     subIntrctn<-intrctn[intrctn$sample==smp &
#                           intrctn$z.origin==combos[cc,1] &
#                           intrctn$z.target==combos[cc,2],]
#     radius<-nrow(subIntrctn)/
#       (prod(
#         length(
#           expGate$primary_id[expGate$sample==smp & 
#                                expGate$primary_id==combos[cc,1] & !expGate$dump]),
#         length(
#           expGate$primary_id[expGate$sample==smp & 
#                                expGate$primary_id==combos[cc,2] & !expGate$dump])))*3000
#     
#     ss<-sample(1:nrow(subIntrctn),nrow(subIntrctn)/100)
#     mdp_O<-c(mean(unlist(subIntrctn[ss,c(4)])),mean(unlist(subIntrctn[ss,c(5)])))
#     mdp<-c(mean(unlist(subIntrctn[ss,c(4,8)])),mean(unlist(subIntrctn[ss,c(5,9)])))
#     mdp_T<-c(mean(unlist(subIntrctn[ss,c(8)])),mean(unlist(subIntrctn[ss,c(9)])))
#     for (ii in ss){
#       points(subIntrctn[ii,4:5],cex=0.1,col=combosPalette[cc])
#       points(subIntrctn[ii,8:9],cex=0.1,col=combosPalette[cc])
#       ll<-matrix(c(mdp_O,
#                    mdp,
#                    mdp_T),ncol=2,byrow=T)
#       lines(ll,col=combosPalette[cc])
#     }
#     tcol<-col2rgb(combosPalette[cc],F)
#     tcol<-rgb(t(tcol)/255,alpha = 0.3)
#     plotrix::draw.circle(mdp[1],mdp[2],
#                          radius = radius,
#                          border=combosPalette[cc],
#                          col=tcol)
#   }
#   dev.off()
# }
# 
# 
# itr<-table(intrctn[,c('z.origin','z.target','sample')])
# 
# ttr<-table(expGate[!expGate$dump,c('primary_id','sample')])
# 
# itr_norm<-lapply(setNames(dimnames(ttr)[[2]],dimnames(ttr)[[2]]),function(x){
#   subItr<-itr[,,x]
#   subttr<-ttr[,x]
#   dama<-subItr
#   dama[,]<-1:length(dama)
#   out<-apply(dama,c(1,2),function(y){
#     xy<-which(dama==y,arr.ind = T)
#     out<-subItr[y]/prod(subttr[xy])*1000
#   })
#   out<-matrix(out,ncol=ncol(subItr),dimnames = dimnames(subItr))
# })
