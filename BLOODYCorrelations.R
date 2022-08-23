dir.create("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/BloodyCorrelation")

bigTable <- read.csv("D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/bigTables/areaTable.csv")
suPperc <- bigTable[,c(1,2,3,4,5)]
suPperc <- suPperc [!duplicated(suPperc[,c(1,2,4)]),]
suPperc$super.area <- suPperc$super.area/suPperc$sample.area*100
suPperc <- aggregate(super.area~sample+super.cluster,suPperc,mean)
suPperc <- as.data.frame(tidyr::pivot_wider(data = suPperc,
                              names_from = 'super.cluster',
                              values_from = 'super.area'))
suPperc <- suPperc[order(suPperc$sample),]
write.csv(suPperc,
          "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/BloodyCorrelation/percentage_SuperPop.csv",
          row.names = FALSE)


suBpercXparental <- bigTable[,c(1:7)]
suBpercXparental$sub.cluster <- paste(suBpercXparental$super.cluster,formatC(suBpercXparental$sub.cluster,width = 3,flag='0'),sep = '_')
suBpercXparental$sub.area <- suBpercXparental$sub.area / suBpercXparental$super.area*100
suBpercXparental <- aggregate(sub.area~sample+sub.cluster,suBpercXparental,mean)
suBpercXparental <- as.data.frame(tidyr::pivot_wider(data = suBpercXparental,
                                            names_from = 'sub.cluster',
                                            values_from = 'sub.area'))
suBpercXparental[is.na(suBpercXparental)]<-0
suBpercXparental <- suBpercXparental[order(suBpercXparental$sample),]
write.csv(suBpercXparental,
          "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/BloodyCorrelation/percentage_SubPopXParental.csv",
          row.names = FALSE)

suBpercXtotal <- bigTable[,c(1:7)]
suBpercXtotal$sub.cluster <- paste(suBpercXtotal$super.cluster,formatC(suBpercXtotal$sub.cluster,width = 3,flag='0'),sep = '_')
suBpercXtotal$sub.area <- suBpercXtotal$sub.area / suBpercXtotal$sample.area*100
suBpercXtotal <- aggregate(sub.area~sample+sub.cluster,suBpercXtotal,mean)
suBpercXtotal <- as.data.frame(tidyr::pivot_wider(data = suBpercXtotal,
                                                     names_from = 'sub.cluster',
                                                     values_from = 'sub.area'))
suBpercXtotal[is.na(suBpercXtotal)]<-0
suBpercXtotal <- suBpercXtotal[order(suBpercXtotal$sample),]
write.csv(suBpercXtotal,
    "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/BloodyCorrelation/percentage_SubPopXTotal.csv",
    row.names = FALSE)


bloodyTable<-bigTable <- read.csv("D:/Luigi/ALL_POSSIBLE_DAICHI_4/IMCSamples_data_PFS_BOR.csv")

bloodyTable <- bloodyTable[bloodyTable$Patient.No.from.exosome.file %in% suBpercXparental$sample,]
bloodyTable <- bloodyTable[order(bloodyTable$Patient.No.from.exosome.file),]

bloodyCols <- c(15:40,47:72)
supCols <- 2:length(suPperc)

ctc <- lapply(bloodyCols,function(i){
  out <- lapply(supCols, function(ii){
    out <-cor.test(bloodyTable[,i],suPperc[,ii])
    out <- data.frame(blood = colnames(bloodyTable)[i],
                      tissue = colnames(suPperc)[ii],
                      correlation = unname(out$estimate),
                      pVal = out$p.value,
                      pVal.s ='',
                      adj.p = 0,
                      adj.s = '')
  })
  out<-do.call(rbind.data.frame,out)
})
ctc <- do.call (rbind.data.frame, ctc)
ctc$pVal.s <- gtools::stars.pval(ctc$pVal)
ctc$adj.p <- p.adjust(ctc$pVal)
ctc$adj.s <- gtools::stars.pval(ctc$adj.p)
ctc <- ctc[order(ctc$pVal,decreasing = F),]

write.csv(ctc,
          "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/BloodyCorrelation/superPop.csv",
          row.names = FALSE)


bloodyCols <- c(15:40,47:72)
supCols <- 2:length(suBpercXparental)

ctc <- lapply(bloodyCols,function(i){
  out <- lapply(supCols, function(ii){
    out <-cor.test(bloodyTable[,i],suBpercXparental[,ii])
    out <- data.frame(blood = colnames(bloodyTable)[i],
                      tissue = colnames(suBpercXparental)[ii],
                      correlation = unname(out$estimate),
                      pVal = out$p.value,
                      pVal.s ='',
                      adj.p = 0,
                      adj.s = '')
  })
  out<-do.call(rbind.data.frame,out)
})
ctc <- do.call (rbind.data.frame, ctc)
ctc$pVal.s <- gtools::stars.pval(ctc$pVal)
ctc$adj.p <- p.adjust(ctc$pVal)
ctc$adj.s <- gtools::stars.pval(ctc$adj.p)
ctc <- ctc[order(ctc$pVal,decreasing = F),]

write.csv(ctc,
          "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/BloodyCorrelation/subPopXparental.csv",
          row.names = FALSE)

bloodyCols <- c(15:40,47:72)
  supCols <- 2:length(suBpercXtotal)

ctc <- lapply(bloodyCols,function(i){
  out <- lapply(supCols, function(ii){
    out <-cor.test(bloodyTable[,i],suBpercXtotal[,ii])
    out <- data.frame(blood = colnames(bloodyTable)[i],
                      tissue = colnames(suBpercXtotal)[ii],
                      correlation = unname(out$estimate),
                      pVal = out$p.value,
                      pVal.s ='',
                      adj.p = 0,
                      adj.s = '')
  })
  out<-do.call(rbind.data.frame,out)
})
ctc <- do.call (rbind.data.frame, ctc)
ctc$pVal.s <- gtools::stars.pval(ctc$pVal)
ctc$adj.p <- p.adjust(ctc$pVal)
ctc$adj.s <- gtools::stars.pval(ctc$adj.p)
ctc <- ctc[order(ctc$pVal,decreasing = F),]

write.csv(ctc,
          "D:/Luigi/ALL_POSSIBLE_DAICHI_4/POST_SEGMENTATION_REVAMP_003-RERUN3/BloodyCorrelation/subPopXtotal.csv",
          row.names = FALSE)
