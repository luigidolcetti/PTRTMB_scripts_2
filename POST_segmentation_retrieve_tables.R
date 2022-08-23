library(sf)
library(FlowSOM)
library(RUNIMCTEMP)

rootFolder<-"C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4"
analFolder<- "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP"

dirList<-list.dirs(rootFolder,recursive = F)
whichDir<-grepl('RUN_s',dirList,fixed = T)
dirList<-dirList[whichDir]
fileList<-paste0(dirList,"/analysis/Anal_001/archive/Expression.RDS")

trnsf_mod<- flowCore::biexponentialTransform(a=0.3,b=0.01,c=0.1,d=5,f=-0.3,w=10)

expTable<-lapply(fileList[[1]],function(fl){
  
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
dst<-dist(t(agg_mixClass_scale))
hc<-hclust(dst)
gr<-cutree(hc,k=7)
gr<-gr[hc$order]
grNames<-vector('character',length(gr))
names(grNames)<-names(gr)
grNames[gr==6]<-"Epithelial"
grNames[gr==2]<-"aSMA"
grNames[gr==1]<-"Unclassified"
grNames[gr==7]<-"CD8"
grNames[gr==5]<-"CD20"
grNames[gr==4]<-"Myeloid"
grNames[gr==3]<-"CD4"


pheatmap::pheatmap(t(agg_mixClass_scale[,hc$order]),
                   cluster_rows = F,
                   cluster_cols = T,
                   annotation_row = data.frame(pixel_level=grNames),
                   scale = 'none')


somMFI<-dplyr::inner_join(expData[,c(1,2,4,9:39),drop=T],somClass,by='pIDs')


##### stripping bad cells #####
colnames(somMFI)[colnames(somMFI) %in% names(newCnames)]<-newCnames

grCode<-data.frame(somClass=names(gr),primary_id=gr)
grCode$primary_id[grCode$primary_id==6]<-"EPTL"
grCode$primary_id[grCode$primary_id==2]<-"SMA"
grCode$primary_id[grCode$primary_id==1]<-"UKNW"
grCode$primary_id[grCode$primary_id==7]<-"CD8"
grCode$primary_id[grCode$primary_id==5]<-"CD20"
grCode$primary_id[grCode$primary_id==4]<-"MYLD"
grCode$primary_id[grCode$primary_id==3]<-"CD4"

compoundClass<-apply(apply(grCode,1,function(x)somMFI$somClass==x[1] & somMFI$primary_id==x[2]),1,any)

somMFI_clean<-somMFI[compoundClass,]

somMFI_median<-aggregate(somMFI[compoundClass,c("area",newCnames)],list(somMFI[compoundClass,'somClass']),median)
somMFI_median<-as.matrix(somMFI_median[,-1])
rownames(somMFI_median)<-levels(somClass$somClass)
somMFI_median<-apply(somMFI_median,2,function(x)(x-min(x))/(max(x)-min(x)))
pheatmap::pheatmap(as.matrix(somMFI_median[hc$order,]),
                   cluster_rows = F,
                   cluster_cols = T,
                   
                   filename = file.path(analFolder,'somMFI.tiff'),
                   cellwidth = 13,
                   cellheight =13,
                   height = 20,
                   width = 10,
                   annotation_row = data.frame(pixel_level=grNames),
                   scale = 'none')
######
######
######
######not stripping bad cells ######
######
######
somMFI_median<-aggregate(somMFI[,c("area",newCnames)],list(somMFI[,'somClass']),median)
somMFI_median<-as.matrix(somMFI_median[,-1])
rownames(somMFI_median)<-levels(somClass$somClass)
somMFI_median<-apply(somMFI_median,2,function(x)(x-min(x))/(max(x)-min(x)))
pheatmap::pheatmap(as.matrix(somMFI_median[hc$order,]),
                   cluster_rows = F,
                   cluster_cols = T,
                   
                   filename = file.path(analFolder,'somMFI.tiff'),
                   cellwidth = 13,
                   cellheight =13,
                   height = 20,
                   width = 10,
                   annotation_row = data.frame(pixel_level=grNames),
                   scale = 'none')

# totalVariance<-apply(somMFI[,c("area",newCnames)],2,sd)
# partialVariance<-aggregate(somMFI[,c("area",newCnames)],list(somMFI$primary_id),sd)
# partialVariance<-lapply(setNames(1:nrow(partialVariance),partialVariance$Group.1),function(i)(totalVariance-partialVariance[i,-1])/totalVariance*100)
# partialVariance<-do.call(rbind,partialVariance)
# tvar<-lapply(setNames(unique(somMFI$primary_id),unique(somMFI$primary_id)),function(x){
#   out<-lapply(setNames(c("area",newCnames),c("area",newCnames)),function(y){
#     out<-lapply(setNames(unique(somMFI$primary_id),unique(somMFI$primary_id)),function(z){
#       out<-var.test(somMFI[somMFI$primary_id==x,y],
#                     somMFI[somMFI$primary_id==z,y],ratio=1,'g')
#       out<-out$p.value
#       return(out)
#     })
#     out<-do.call(c,out)
#   })
#   out<-do.call(rbind,out)
# })
# tvar<-lapply(tvar,function(x)apply(x,1,function(y)sum(y<0.05)/ncol(x)*100))
# colOfInterest<-lapply(tvar,function(x)names(x)[x>=60])

colOfInterest<-list(CD4=c('CD45RA','CD45RO','CD27','CD127','CD25','FoxP3'),
                    CD8=c('CD45RA','CD45RO','CD27','CD127'),
                    MYLD=c('CD33','CD11b','CD14','CD15','CD68','CD74','CD16','CD38','Arginase_1','CD11c'),
                    CD20=c('CD20','CD11b','IgD','CD38'),
                    EPTL=c('Pan-keratin'),
                    SMA=c(newCnames),
                    UKNW=c(newCnames))

for (i in names(colOfInterest)){
  colToUse<-c("area",newCnames)
  colmin<-apply(somMFI_clean[,colToUse],2,min)
  colmax<-apply(somMFI_clean[,colToUse],2,max)
  somSubset<-grCode$somClass[grCode$primary_id==i]
  somMFI_median<-aggregate(somMFI_clean[somMFI_clean$somClass %in% somSubset,c("area",newCnames)],
                           list(somMFI_clean[somMFI_clean$somClass %in% somSubset,'somClass']),median)
  TEMP<-somMFI_median[,1]
  somMFI_median<-as.matrix(somMFI_median[,-1])
  rownames(somMFI_median)<-TEMP
  somMFI_median<-do.call(cbind,lapply(setNames(colToUse,colToUse),function(x)(somMFI_median[,x]-colmin[x])/(colmax[x]-colmin[x])))
  
  mtx<-somMFI_median[,colOfInterest[[i]]]
  dst<-dist(mtx)
  hct<-hclust(dst)
  
  pheatmap::pheatmap(as.matrix(somMFI_median[hct$order,c("area",newCnames)]),
                     cluster_rows = F,
                     cluster_cols = F,
                     filename = file.path(analFolder,paste0('somMFI_',i,'.tiff')),
                     cellwidth = 13,
                     cellheight =13,
                     height = if(nrow(somMFI_median)<10) nrow(somMFI_median)/2 else nrow(somMFI_median)/3,
                     width = 15,
                     main = i,
                     scale = 'none')
}


ump<-umap::umap(expData[compoundClass,9:39,drop=T])

saveRDS(ump,file.path(analFolder,'ump.R'))

plot(ump$layout,cex=0.1)
points(ump$layout[somMFI$somClass %in% c('052','062'),],col='red',cex=0.1)
points(ump$layout[expData$primary_id[compoundClass]=='SMA',],col='cyan',cex=0.1)
