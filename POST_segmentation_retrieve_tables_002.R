library(sf)
library(FlowSOM)
library(RUNIMCTEMP)
library(pbapply)
library(bezier)
##### declare vars ####

rootFolder<-"C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4"
analFolder<- "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP"

dirList<-list.dirs(rootFolder,recursive = F)
whichDir<-grepl('RUN_s',dirList,fixed = T)
dirList<-dirList[whichDir]
fileList<-paste0(dirList,"/analysis/Anal_001/archive/Expression.RDS")

trnsf_mod<- flowCore::biexponentialTransform(a=0.3,b=0.01,c=0.1,d=5,f=-0.3,w=10)

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

write_sf(expData,
         file.path(analFolder,'fullSet_biexp_trans.sqlite'))

rm(expTable)

gc()

##### retrieve study table #####

fileList<-paste0(dirList,"/archive/studyTable.RDS")

studyTable<-lapply(fileList,function(fl){
  out<-readRDS(fl)
  return(out)
})

studyTable<-do.call(rbind.data.frame,studyTable)

write.csv(studyTable,file.path(analFolder,'studyTable.csv'))

expData<- dplyr::full_join(expData,studyTable,by="uid")

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

saveRDS(som,
        file.path(analFolder,"som.rds"))

somClass<-GetClusters(som)

somClass<-formatC(somClass,width=3,flag='0')

somClass<-as.factor(somClass)

somClass<-cbind.data.frame(pIDs=expData$pIDs,somClass=somClass)

write.csv(somClass,
          file.path(analFolder,'som_class.csv'))

mixClass<-dplyr::inner_join(expData[,c('pIDs','area','primary_id'),drop=T],somClass,by='pIDs')

agg_mixClass<-aggregate(area~primary_id+somClass,mixClass,sum)

agg_mixClass$primary_id<-factor(agg_mixClass$primary_id)
agg_mixClass$somClass<-factor(agg_mixClass$somClass)

agg_mixClass<-lapply(setNames(levels(agg_mixClass$primary_id),levels(agg_mixClass$primary_id)),function(x){
  out<-lapply(setNames(levels(agg_mixClass$somClass),levels(agg_mixClass$somClass)),function(y){
    out<-agg_mixClass[agg_mixClass$primary_id==x & agg_mixClass$somClass==y,'area']
    if (length(out)==0) out<-0
    return(out)
  })
  out<-matrix(unlist(out),nrow = 1,dimnames = list(Pixel_level=x,Cell_level=levels(agg_mixClass$somClass)))
})
agg_mixClass<-do.call(rbind,agg_mixClass)

agg_mixClass_scale<-scale(agg_mixClass)
dst<-dist(t(agg_mixClass_scale),method = 'maximum')
hc<-hclust(dst,method = 'complete')
gr<-cutree(hc,k=7)
gr<-gr[hc$order]
grNames<-vector('character',length(gr))
names(grNames)<-names(gr)
grNames[gr==1]<-"Myeloid"
grNames[gr==2]<-"Epithelial"
grNames[gr==3]<-"CD8+"
grNames[gr==4]<-"CD4+"
grNames[gr==5]<-"Mesenchymal"
grNames[gr==6]<-"Unclassified"
grNames[gr==7]<-"CD20+"
grPalette<- setNames(RColorBrewer::brewer.pal(length(unique(grNames)),name = 'Set2'),unique(grNames))
grAnnotation<-data.frame(Composite=as.factor(grNames))
grAnnotation_color<-list(Composite=grPalette)

##### heat map som ####

pheatmap::pheatmap(t(agg_mixClass_scale),
                   cluster_rows = T,
                   cluster_cols = T,
                   cutree_rows = 7,
                   clustering_distance_rows = 'maximum',
                   clustering_method = 'complete',
                   annotation_row = grAnnotation,
                   annotation_colors = grAnnotation_color,
                   labels_col = c('CD20+','CD4+','CD8+','Epithelial','Myeloid','Mesenchymal','Unclassifed'),
                   filename = file.path(analFolder,'somMFIxforest.pdf'),
                   cellwidth = 13,
                   cellheight =13,
                   height = 22,
                   width = 6,
                   main = 'Coposite clusters',
                   scale = 'none')

somMFI<-dplyr::inner_join(expData[,c(1,2,4,9:39),drop=T],somClass,by='pIDs')

##### stripping bad cells #####

colnames(somMFI)[colnames(somMFI) %in% names(newCnames)]<-newCnames

grCode<-data.frame(somClass=names(gr),primary_id=gr)

# grNames[gr==1]<-"Myeloid"
# grNames[gr==2]<-"Epithelial"
# grNames[gr==3]<-"CD8+"
# grNames[gr==4]<-"CD4+"
# grNames[gr==5]<-"Mesenchymal"
# grNames[gr==6]<-"Unclassified"
# grNames[gr==7]<-"CD20+"

grCode$primary_id[grCode$primary_id==1]<-"MYLD"
grCode$primary_id[grCode$primary_id==2]<-"EPTL"
grCode$primary_id[grCode$primary_id==3]<-"CD8"
grCode$primary_id[grCode$primary_id==4]<-"CD4"
grCode$primary_id[grCode$primary_id==5]<-"SMA"
grCode$primary_id[grCode$primary_id==6]<-"UKNW"
grCode$primary_id[grCode$primary_id==7]<-"CD20"


##### calculate compound class  and uMAP#####

compoundClass<-apply(apply(grCode,1,function(x)somMFI$somClass==x[1] & somMFI$primary_id==x[2] & somMFI$primary_id!='MIX'),1,any)

somMFI_clean<-somMFI[compoundClass,]

ump<-umap::umap(expData[compoundClass,9:39,drop=T])

saveRDS(ump,file.path(analFolder,'ump.R'))

ump_total<- predict(ump,expData[,9:39,drop=T])

CC_table<-table(compoundClass)
CC_ratio<-CC_table[1]/CC_table[2]
ntot<-5000
CC_in<-sample(1:CC_table[2],ntot)
CC_out<-sample(1:CC_table[1],round(ntot*CC_ratio))
postscript(file = file.path(analFolder,'Umap_total.eps'),
           onefile = F,
           width = 7.5,
           height = 7.5,
           horizontal = F,
           paper = 'special',
           bg='white')

plot(ump_total,
     pch=16,
     cex=0.8,
     col='black',
     xlab='Umap1',ylab='Umap2')

legend('topleft',
       legend = c('Retained events','Rejected events'),
       pch=16,
       col = c('darkolivegreen4','firebrick'),
       cex=0.8,
       bty = 'n',
)
points(ump_total[compoundClass,][CC_in,],
       pch=16,
       cex=0.5,
       col='darkolivegreen4',
       xlab='Umap1',ylab='Umap2')

points(ump_total[!compoundClass,][CC_out,],
       pch = 16, 
       cex=0.5,
       col='firebrick')
dev.off()

postscript(file = file.path(analFolder,'Umap_primary_ID.eps'),
           onefile = F,
           width = 7.5,
           height = 7.5,
           horizontal = F,
           paper = 'special',
           bg='white')

grPalette_1<-grPalette[c(7,2,1,3,5,6,4)]

plot(ump_total[compoundClass,],
     pch=16,
     cex=0.5,
     col='black',
     xlab='Umap1',ylab='Umap2')

legend('topleft',
       legend = names(grPalette_1),
       pch=16,
       col = grPalette_1,
       cex=0.8,
       inset = c(0,0),
       bty = 'n')

for (i in 1:length(unique(expData$primary_id))){
  points(ump_total[compoundClass & expData$primary_id==unique(expData$primary_id)[i],],
         pch=16,
         cex=0.1,
         col=grPalette_1[i],
         xlab='Umap1',ylab='Umap2')
  
}
dev.off()


###### detailed heatmaps ######

colOfInterest<-list(CD4=c('CD45RA','CD45RO','CD27','CD127','CD25','FoxP3'),
                    CD8=c('CD45RA','CD45RO','CD27','CD127'),
                    MYLD=c('CD33','CD11b','CD14','CD15','CD68','CD74','CD16','CD38','Arginase_1','CD11c'),
                    CD20=c('CD20','CD11b','IgD','CD38'),
                    EPTL=c('Pan-keratin','CD68','eCadherin','Arginase_1'),
                    SMA=c(newCnames),
                    UKNW=c(newCnames))

mrkr<-unique(expData$primary_id)
colToUse<-c("area",newCnames)
somList<-vector('list',length(mrkr))
names(somList)<-mrkr
for (i in 1:length(mrkr)){
  
  subMatrix<-expMatrix[compoundClass & expData$primary_id==mrkr[i],newCnames]
  
  som<-ReadInput(input = subMatrix,
                 compensate = F,
                 transform = F,
                 scale = T,
                 silent = F)
  
  som<-BuildSOM(fsom = som,
                colsToUse = colOfInterest[[mrkr[i]]],
                xdim = 3,
                ydim = 3)
  
  saveRDS(som,
          file.path(analFolder,paste0("som_",mrkr[i],".rds")))
  
  somClass<-GetClusters(som)
  
  somClass<-formatC(somClass,width=3,flag='0')
  
  somClass<-as.factor(somClass)
  
  somClass<-cbind.data.frame(pIDs=expData$pIDs[compoundClass & expData$primary_id==mrkr[i]] ,somClass=somClass)
  
  somList[[mrkr[i]]]<-somClass
  
  write.csv(somClass,
            file.path(analFolder,paste0("somTable",mrkr[i],".csv")))
  
  somMFI<-dplyr::inner_join(expData[,c(1,2,4,9:39),drop=T],somClass,by='pIDs')
  
  colnames(somMFI)[colnames(somMFI) %in% names(newCnames)]<-newCnames
  
  colToUse<-c("area",newCnames)
  
  
  somMFI_median<-aggregate(somMFI[,colToUse],list(somMFI$somClass),median)
  TEMP<-somMFI_median[,1]
  somMFI_median<-as.matrix(somMFI_median[,-1])
  rownames(somMFI_median)<-TEMP
  
  colmin<-apply(expData[,c(9:39),drop=T],2,quantile,0.01)
  names(colmin)<-colToUse
  colmax<-apply(expData[,c(9:39),drop=T],2,quantile,0.99)
  names(colmax)<-colToUse
  somMFI_median<-do.call(cbind,lapply(setNames(colnames(somMFI_median),colnames(somMFI_median))
                                      ,function(x)(somMFI_median[,x]-colmin[x])/(colmax[x]-colmin[x])))
  
  pheatmap::pheatmap(t(somMFI_median),
                     cluster_rows = F,
                     cluster_cols = T,
                     filename = file.path(analFolder,paste0("som_",mrkr[i],".pdf")),
                     cellwidth = 13,
                     cellheight =13,
                     height = 8,
                     width = 5,
                     scale = 'none',
                     main = mrkr[i])
  
}

saveRDS(somList,file.path(analFolder,'Sub_cluster_Som_list.R'))

limitList = list(CD20=c(-10,-8,-6,-4.5),
                 CD4 = c(-12,-7,-5,0),
                 CD8 = c(-16,-12.5,-5,0),
                 EPTL = c(-6,9,2.5,16),
                 MYLD = c(-7.5,0,-5,3),
                 SMA = c(0,5,-13,-8),
                 UKNW = c(0,10,-7,2.5))

for (i in 1:length(mrkr)){
  
  somClass<-somList[[i]]
  uPalette<- setNames(c(RColorBrewer::brewer.pal(length(levels(somClass$somClass)),name = 'Set2'),'dodgerblue3'),levels(somClass$somClass))
  
  postscript(file = file.path(analFolder,paste0('Umap_',mrkr[i],'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  
  plot(ump_total,
       pch=16,
       cex=0.6,
       col='black',
       xlim=c(limitList[[i]][1:2]),
       ylim=c(limitList[[i]][3:4]),
       xlab='Umap1',ylab='Umap2')
  
  legend('topleft',
         legend = names(uPalette),
         pch=16,
         col = uPalette,
         cex=0.8,
         inset = c(0,0),
         bty = 'n')
  
  for (ii in unique(somClass$somClass)){
    
    pIDs<-somClass$pIDs[somClass$somClass==ii]
    focusPoint<-expData$pIDs %in% pIDs
    points(ump_total[focusPoint,],
           pch=16,
           cex=0.5,
           col=uPalette[ii],
           xlab='Umap1',ylab='Umap2')
  }
  dev.off()
}

# subSetList<-list(
#   CD4 = list(`Treg` = c("003","006"),
#              `Treg/CD25low` = c("009"),
#              `EM` = c("007","008"),
#              `EM-CD45RA+` = c("001","005"),
#              `Naive` = c("002"),
#              `EFF` = c("004")),
#   CD8 = list(`CM` = c("009"),
#              `EM` = c("008"),
#              `Naive` = c("007"),
#              `EFF` = c("001","002","003","004","005","006")),
#   CD20 = list(`CD38+` = c("003"),
#               `IgD+` = c("002"),
#               `CD38-CD24-` = c("001","004","005","006","007","008","009")),
#   EPTL = list(`PANK+` = c("007"),
#               `PANK-` = c("001","002","003","004","005","006","008","009")),
#   MYLD = list(`Granulo` = c("004","007"),
#               `Macro` = c("001"),
#               `Mono` = c("002","005"),
#               `DC` = c("006"),
#               `MDSC` = c("003","009"),
#               `Mono-ARG+` = c("008")),
#   SMA = list(`DNA+` = c("006","009"),
#              `DNA-` = c("001","002","003","004","005","007","008")),
#   UKNW = list(`CD27+CD45RA+`=c("008","009"),
#               `CD27+CD45RA+`=c("001","002","003","004","005","006","007"))
# )

##### SUPER POP percentages ####
totArea<-aggregate(area~uid,expData,sum)

colToKeep<-c('uid','sample','bioGroup','pIDs','area')
baseFrame<-data.frame(somClass=sort(unique(somClass$somClass)),area=0)

aggArea<-lapply(setNames(unique(expData$uid),unique(expData$uids)),function(u){
  expSubSet<-expData[expData$uid==u,colToKeep,drop=T]
  smp<-unique(expSubSet$sample)
  pfs<-unique(expSubSet$bioGroup)
  totArea<-sum(expSubSet$area)
  out<-lapply(setNames(names(somList),names(somList)),function(ss){
    expSSubSet<-dplyr::inner_join(expSubSet,somList[[ss]],'pIDs')
    superArea<-sum(expSSubSet$area)
    out<-data.frame(uid=u,
                    sample = smp,
                    pfs = pfs,
                    super.cluster = ss,
                    sub.cluster = baseFrame$somClass,
                    sample.area = totArea,
                    super.area = superArea,
                    sub.area = baseFrame$area)
    if (nrow(expSSubSet)==0){
      return(out)
    } else {
      somArea<-aggregate(area~somClass,expSSubSet,sum)
      out$sub.area[out$sub.cluster %in% somArea$somClass]<-somArea$area
      return(out)}
  })
  out<-do.call(rbind.data.frame,out)
  out<-cbind.data.frame(out[,1:6],retained.area=sum(out$sub.area),out[,7:8])
  return(out)
})

aggArea<-do.call(rbind.data.frame,aggArea)

percentage_super<-aggArea[,-c(5,9)]
percentage_super<-unique(percentage_super)

av_DF<-data.frame(pfs=percentage_super$pfs,
                  smp=percentage_super$sample,
                  class=percentage_super$super.cluster,
                  area=percentage_super$super.area/percentage_super$sample.area*100)

av<-aov(area~smp:class,av_DF)

av_HSD<-TukeyHSD(av)

rnms<-rownames(av_HSD$`smp:class`)
av_first_smp<-unlist(lapply(rnms,function(x){
  rgx<-regexpr("^([0-9]+)\\:[A-Z0-9\\_]+.",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_first_som<-unlist(lapply(rnms,function(x){
  rgx<-regexpr("^[0-9]+\\:([A-Z0-9\\_]+).",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_second_smp<-unlist(lapply(rnms,function(x){
  rgx<-regexpr(".+\\-([0-9]+):[A-Z0-9\\_]+$",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_second_som<-unlist(lapply(rnms,function(x){
  rgx<-regexpr(".+\\-[0-9]+:([A-Z0-9\\_]+)$",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_results<-data.frame(av_first_smp,
                       av_first_som,
                       av_second_smp,
                       av_second_som,
                       p.adj=av_HSD$`smp:class`[,"p adj"],
                       s.adj= gtools::stars.pval(av_HSD$`smp:class`[,"p adj"]))

av_results<-av_results[apply(av_results,1,function(x) x['av_first_som']==x['av_second_som']),]

smpOrder<-unique(percentage_super[,c('sample','pfs')])
smpOrder<-factor(smpOrder$sample[order(smpOrder$pfs,decreasing = F)],
                 levels=smpOrder$sample[order(smpOrder$pfs,decreasing = F)],
                 ordered = T)

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
          main=i,
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
  subSet<-av_results[av_results$av_first_som==i,]
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

##### SUB POP percentages ####


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


av<-aov(area~smp:class,av_DF)

av_HSD<-TukeyHSD(av)

rnms<-rownames(av_HSD$`smp:class`)
av_first_smp<-unlist(lapply(rnms,function(x){
  rgx<-regexpr("^([0-9]+)\\:[A-Z0-9\\_]+.",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_first_som<-unlist(lapply(rnms,function(x){
  rgx<-regexpr("^[0-9]+\\:([A-Z0-9\\_]+).",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_second_smp<-unlist(lapply(rnms,function(x){
  rgx<-regexpr(".+\\-([0-9]+):[A-Z0-9\\_]+$",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_second_som<-unlist(lapply(rnms,function(x){
  rgx<-regexpr(".+\\-[0-9]+:([A-Z0-9\\_]+)$",x,perl = T)
  substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
}))

av_results<-data.frame(av_first_smp,
                       av_first_som,
                       av_second_smp,
                       av_second_som,
                       p.adj=av_HSD$`smp:class`[,"p adj"],
                       s.adj= gtools::stars.pval(av_HSD$`smp:class`[,"p adj"]))

av_results<-av_results[apply(av_results,1,function(x) x['av_first_som']==x['av_second_som']),]

smpOrder<-unique(percentage_super[,c('sample','pfs')])
smpOrder<-factor(smpOrder$sample[order(smpOrder$pfs,decreasing = F)],
                 levels=smpOrder$sample[order(smpOrder$pfs,decreasing = F)],
                 ordered = T)

for (i in unique(av_DF$class)){
  subSet<-av_DF$area[av_DF$class==i]
  names(subSet)<-av_DF$smp[av_DF$class==i]
  smp<-factor(names(subSet),levels = levels(smpOrder),ordered = T)
  df_box<-data.frame(Sample=smp,
                     `% of total area` = subSet,
                     check.names = F)
  postscript(file = file.path(analFolder,"SUB.POP",paste0(i,'.eps')),
             onefile = F,
             width = 7.5,
             height = 7.5,
             horizontal = F,
             paper = 'special',
             bg='white')
  
  par(mfrow=c(2,1),mar=c(0,5,2,1),las=1,bty='l')
  
  boxplot(`% of total area`~Sample,
          data = df_box,
          main=i,
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
  subSet<-av_results[av_results$av_first_som==i,]
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

##### INTERACTIONS ######
intrctn<-lapply(setNames(studyTable$uid,studyTable$uid),function(u){
  
  subSet_data<-expData[expData$uid==u & compoundClass,'geom']
  primaryIDs<-expData[expData$uid==u & compoundClass,'primary_id',drop=T]
  pIDs<-expData[expData$uid==u & compoundClass,'pIDs',drop=T]
  smp<-unique(studyTable$sample[studyTable$uid==u])
  subSet_ump<-ump_total[expData$uid==u & compoundClass,]
  intrctn<-st_intersects(subSet_data)
  
  out<-lapply(intrctn,function(x){
    if (length(x)>1){
      out<-data.frame(uid=u,
                      sample=smp,
                      id.origin = pIDs[x[1]],
                      x.origin = subSet_ump[x[1],1],
                      y.origin = subSet_ump[x[1],2],
                      z.origin = primaryIDs[x[1]],
                      id.target = pIDs[x[-1]],
                      x.target = subSet_ump[x[-1],1],
                      y.target = subSet_ump[x[-1],2],
                      z.target = primaryIDs[x[-1]],
                      row.names = NULL,
                      stringsAsFactors = F)}
    else {
      out<-data.frame(uid=character(0),
                      sample=character(0),
                      id.origin = character(0),
                      x.origin = numeric(0),
                      y.origin = numeric(0),
                      z.origin = character(0),
                      id.target = character(0),
                      x.target = numeric(0),
                      y.target = numeric(0),
                      z.target = character(0),
                      row.names = NULL,
                      stringsAsFactors = F)
    }
    return(out)
  })
  out<-do.call(rbind,out)
})
intrctn<-do.call(rbind,intrctn)

plot(ump$layout,col='gray80',cex=0.3)
points(ump$layout,col='white',cex=0.1)
subIntrctn<-intrctn[intrctn$sample=='40040006' &
                      intrctn$z.origin=='EPTL' &
                      intrctn$z.target=='SMA',]
ss<-sample(1:nrow(subIntrctn),nrow(subIntrctn)/100)
mdp<-c(mean(unlist(subIntrctn[ss,c(4,8)])),mean(unlist(subIntrctn[ss,c(5,9)])))
t<-seq(0,1,length=100)
for (ii in ss){
  ll<-matrix(c(subIntrctn[ii,4:5],
               mdp,
               subIntrctn[ii,8:9]),ncol=2,byrow=T)
  if (any(diff(unlist(ll[,1]))<0)){
    spl<-spline(ll[,2:1])
    spl<-spl[2:1]
    names(spl)<-c('x','y')
  } else {
    spl<-spline(ll)
  }
  lines(spl)
 }
# av_long<- matrix(unlist(av_DF$area),ncol=9*7,byrow = T)
# colnames(av_long)<-unique(av_DF$class)
# rownames(av_long)<-unique(expData$uid)
# 
# av_long[is.na(av_long)]<-0
# 
# colnames(av_long)[19:27]
# pc<- prcomp(av_long)
# tPalette<-setNames(c(RColorBrewer::brewer.pal(length(unique(studyTable$sample)),name = 'Set2'),'dodgerblue3'),unique(studyTable$sample))
# tbg<-colorRampPalette(c('blue','white','red'))(255)
# tbgv<-as.numeric(studyTable$bioGroup)
# tbgv<-round((tbgv-min(tbgv))/(max(tbgv)-min(tbgv))*254)+1
# plot(pc$x[,c(1,3)],
#      col='black',
#      bg='black',
#      pch=21,
#      cex=2.2)
# points(pc$x[,c(1,3)],
#      col=tPalette[studyTable$sample],
#      bg=tPalette[studyTable$sample],
#      pch=21,
#      cex=2)
# points(pc$x[,c(1,3)],
#      col=tbg[tbgv],
#      bg=tbg[tbgv],
#      pch=21,
#      cex=0.8)
# 
# 
# library(factoextra)
#   
# fviz_pca_ind(pc,
#              label='none',
#              col.ind = as.numeric(studyTable$bioGroup))
# fviz_pca_var(pc)
# av_DF<-data.frame(pfs=aggArea$pfs,
#                  smp=aggArea$sample,
#                  class=paste0(aggArea$super.cluster,'_',aggArea$sub.cluster),
#                  area=aggArea$sub.area/aggArea$sample.area*100)
# 
# av<-aov(area~smp:class,av_DF)
# 
# av_HSD<-TukeyHSD(av)
# 
# rnms<-rownames(av_HSD$`smp:class`)
# av_first_smp<-unlist(lapply(rnms,function(x){
#   rgx<-regexpr("^([0-9]+)\\:[A-Z0-9\\_]+.",x,perl = T)
#   substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
# }))
# 
# av_first_som<-unlist(lapply(rnms,function(x){
#   rgx<-regexpr("^[0-9]+\\:([A-Z0-9\\_]+).",x,perl = T)
#   substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
# }))
# 
# av_second_smp<-unlist(lapply(rnms,function(x){
#   rgx<-regexpr(".+\\-([0-9]+):[A-Z0-9\\_]+$",x,perl = T)
#   substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
# }))
# 
# av_second_som<-unlist(lapply(rnms,function(x){
#   rgx<-regexpr(".+\\-[0-9]+:([A-Z0-9\\_]+)$",x,perl = T)
#   substr(x,attr(rgx,'capture.start'),attr(rgx,'capture.start')+attr(rgx,'capture.length')-1)
# }))
# 
# av_results<-data.frame(av_first_smp,
#                        av_first_som,
#                        av_second_smp,
#                        av_second_som,
#                        p.adj=av_HSD$`smp:class`[,"p adj"],
#                        s.adj= gtools::stars.pval(av_HSD$`smp:class`[,"p adj"]))
# 
# TEMP<-av_results[apply(av_results,1,function(x) x['av_first_som']==x['av_second_som']),]
# TEMP[order(TEMP$p.adj,decreasing = F),]
# 
# 
# percentages_super<-lapply(setNames(names(subSetList),names(subSetList)),function(i){
#   
#   subSet<-expData[expData$pIDs %in% somList[[i]]$pIDs,c('uid','pIDs','area'),drop=T]
#   subSet<-dplyr::left_join(subSet,somList[[i]],'pIDs')
#   partialArea<-aggregate(area~uid,subSet,sum)
#   out<-dplyr::left_join(partialArea,totArea,by='uid',suffix=c(".super",".total"))
#   newCluster<-vector('character',nrow(subSet))
#   for (ii in 1:length(subSetList[[i]])){
#     newCluster[subSet$somClass %in% subSetList[[i]][[ii]]]<-names(subSetList[[i]])[[ii]]
#   }
#   browser()
#   subSet<-cbind.data.frame(subSet,as.factor(newCluster))
#   
#   
#   PartialArea<-aggregate(area~newCluster+uid,subSet,sum)
#   PartialArea<- tidyr::pivot_wider(data = PartialArea,
#                                    values_from = 'area',
#                                    names_from = 'newCluster')
#   out<-dplyr::left_join(out,PartialArea,by='uid')
#   out<-cbind.data.frame(Super_class=i,out)
#   out[is.na(out)]<-0
#   return(out)
# })
# 
# # 
# # 
# # plot(ump$layout,cex=0.1)
# # for(i)
# #   points(ump$layout[somMFI$somClass %in% c('052','062'),],col='red',cex=0.1)
# # points(ump$layout[expData$primary_id[compoundClass]=='SMA',],col='cyan',cex=0.1)
# # 
# # 
# # 
# # somMFI_median<-aggregate(somMFI[compoundClass,c("area",newCnames)],list(somMFI[compoundClass,'somClass']),median)
# # somMFI_median<-as.matrix(somMFI_median[,-1])
# # rownames(somMFI_median)<-levels(somClass$somClass)
# # somMFI_median<-apply(somMFI_median,2,function(x)(x-min(x))/(max(x)-min(x)))
# # pheatmap::pheatmap(as.matrix(somMFI_median[hc$order,]),
# #                    cluster_rows = F,
# #                    cluster_cols = T,
# #                    
# #                    filename = file.path(analFolder,'somMFI.tiff'),
# #                    cellwidth = 13,
# #                    cellheight =13,
# #                    height = 20,
# #                    width = 10,
# #                    annotation_row = data.frame(pixel_level=grNames),
# #                    scale = 'none')
# # ######
# # ######
# # ######
# # ######not stripping bad cells ######
# # ######
# # ######
# # somMFI_median<-aggregate(somMFI[,c("area",newCnames)],list(somMFI[,'somClass']),median)
# # somMFI_median<-as.matrix(somMFI_median[,-1])
# # rownames(somMFI_median)<-levels(somClass$somClass)
# # somMFI_median<-apply(somMFI_median,2,function(x)(x-min(x))/(max(x)-min(x)))
# # pheatmap::pheatmap(as.matrix(somMFI_median[hc$order,]),
# #                    cluster_rows = F,
# #                    cluster_cols = T,
# #                    
# #                    filename = file.path(analFolder,'somMFI.tiff'),
# #                    cellwidth = 13,
# #                    cellheight =13,
# #                    height = 20,
# #                    width = 10,
# #                    annotation_row = data.frame(pixel_level=grNames),
# #                    scale = 'none')
# # 
# # # totalVariance<-apply(somMFI[,c("area",newCnames)],2,sd)
# # # partialVariance<-aggregate(somMFI[,c("area",newCnames)],list(somMFI$primary_id),sd)
# # # partialVariance<-lapply(setNames(1:nrow(partialVariance),partialVariance$Group.1),function(i)(totalVariance-partialVariance[i,-1])/totalVariance*100)
# # # partialVariance<-do.call(rbind,partialVariance)
# # # tvar<-lapply(setNames(unique(somMFI$primary_id),unique(somMFI$primary_id)),function(x){
# # #   out<-lapply(setNames(c("area",newCnames),c("area",newCnames)),function(y){
# # #     out<-lapply(setNames(unique(somMFI$primary_id),unique(somMFI$primary_id)),function(z){
# # #       out<-var.test(somMFI[somMFI$primary_id==x,y],
# # #                     somMFI[somMFI$primary_id==z,y],ratio=1,'g')
# # #       out<-out$p.value
# # #       return(out)
# # #     })
# # #     out<-do.call(c,out)
# # #   })
# # #   out<-do.call(rbind,out)
# # # })
# # # tvar<-lapply(tvar,function(x)apply(x,1,function(y)sum(y<0.05)/ncol(x)*100))
# # # colOfInterest<-lapply(tvar,function(x)names(x)[x>=60])
# # 
# # colOfInterest<-list(CD4=c('CD45RA','CD45RO','CD27','CD127','CD25','FoxP3'),
# #                     CD8=c('CD45RA','CD45RO','CD27','CD127'),
# #                     MYLD=c('CD33','CD11b','CD14','CD15','CD68','CD74','CD16','CD38','Arginase_1','CD11c'),
# #                     CD20=c('CD20','CD11b','IgD','CD38'),
# #                     EPTL=c('Pan-keratin'),
# #                     SMA=c(newCnames),
# #                     UKNW=c(newCnames))
# # 
# # for (i in names(colOfInterest)){
# #   colToUse<-c("area",newCnames)
# #   colmin<-apply(somMFI_clean[,colToUse],2,min)
# #   colmax<-apply(somMFI_clean[,colToUse],2,max)
# #   somSubset<-grCode$somClass[grCode$primary_id==i]
# #   somMFI_median<-aggregate(somMFI_clean[somMFI_clean$somClass %in% somSubset,c("area",newCnames)],
# #                            list(somMFI_clean[somMFI_clean$somClass %in% somSubset,'somClass']),median)
# #   TEMP<-somMFI_median[,1]
# #   somMFI_median<-as.matrix(somMFI_median[,-1])
# #   rownames(somMFI_median)<-TEMP
# #   somMFI_median<-do.call(cbind,lapply(setNames(colToUse,colToUse),function(x)(somMFI_median[,x]-colmin[x])/(colmax[x]-colmin[x])))
# #   
# #   mtx<-somMFI_median[,colOfInterest[[i]]]
# #   dst<-dist(mtx)
# #   hct<-hclust(dst)
# #   
# #   pheatmap::pheatmap(as.matrix(somMFI_median[hct$order,c("area",newCnames)]),
# #                      cluster_rows = F,
# #                      cluster_cols = F,
# #                      filename = file.path(analFolder,paste0('somMFI_',i,'.tiff')),
# #                      cellwidth = 13,
# #                      cellheight =13,
# #                      height = if(nrow(somMFI_median)<10) nrow(somMFI_median)/2 else nrow(somMFI_median)/3,
# #                      width = 15,
# #                      main = i,
# #                      scale = 'none')
# # }
# # 
# # 
# # ump<-umap::umap(expData[compoundClass,9:39,drop=T])
# # 
# # saveRDS(ump,file.path(analFolder,'ump.R'))
# # 
# # plot(ump$layout,cex=0.1)
# # points(ump$layout[somMFI$somClass %in% c('052','062'),],col='red',cex=0.1)
# # points(ump$layout[expData$primary_id[compoundClass]=='SMA',],col='cyan',cex=0.1)
