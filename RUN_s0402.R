startTime<-Sys.time()
library(RUNIMCTEMP)

out<-try({
  myStudy<-initStudy(fn_studyName = 'RUN_s0402',
                     fn_rootFolder = "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4",
                     fn_rawDataFolder = "D:/Luigi/DAICHI_FULL_SET",
                     fn_whichFiles = c(75:92),
                     fn_whichColumns = 'named')
  
  newAnalysis(myStudy,'Anal_001')
  
  myStudy$currentAnalysis$classifier<-retrieve(
    "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4/RUN_model_creation_reduced/analysis/Anal_001/archive/Classifier.RDS"
  )
  
  myStudy$currentAnalysis$classificationDirectives<-retrieve(
    "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4/RUN_model_creation_reduced/analysis/Anal_001/archive/ClassificationDirectives.RDS"
  )
  
  myStudy$currentAnalysis$filters<-retrieve(
    "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4/RUN_model_creation_reduced/analysis/Anal_001/archive/filters.RDS"
  )
  
  myStudy$currentAnalysis$extractionDirectives<-retrieve(
    "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4/RUN_model_creation_reduced/analysis/Anal_001/archive/ExtractionDirectives.RDS"
  )
  
  myStudy$currentAnalysis$trainingFeatures<-retrieve(
    "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4/RUN_model_creation_reduced/analysis/Anal_001/archive/TrainingFeatures.RDS"
  )
  
  polyList<-list.files(
    "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4/backup_polygons/polygons",
    full.names = T
  )
  
  file.copy(from = polyList,
            to = file.path(myStudy$currentAnalysis$folder,'training/polygons/'))
  
  deployFilters(x = myStudy,cl = 7)
  
  nsmpls<-length(st_samples(myStudy))
  
  if (nsmpls>6) nsmpls<-6
  
  classify(myStudy,method = 'randomForest',cl = nsmpls)
  
  localCorrection(myStudy,
                  cl = nsmpls,
                  matrixExtent = 3)
  
  classify(myStudy, cl = nsmpls,method = 'randomOnions')
  
  addSegmentationDirectives(myStudy,method = 'pandaMap')
  
  myStudy$currentAnalysis$segmentationDirectives
  
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$lowerQuantile = 0
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$upperQuantile = 1
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$lowerAreaLimit = 5
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$numberOfCores = 12
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$ClampDetectionDirection = '4'
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$overlapExtent = 50
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$movingWindowDimension =180:220
  myStudy$currentAnalysis$segmentationDirectives@methodParameters$nOfCutBrakes = 12
  
  avalLayer<-layerNames(myStudy)
  
  segment(myStudy,
          labelLayer = avalLayer$class[c(11:17)])
  
  
  
  topLayers<-lapply(st_uids(myStudy),function(uids){
    out<-RUNIMCTEMP:::.topLayer(x = myStudy$currentAnalysis$exprs@exprs$pandaMap,
                                uid = uids)})
  
  topLayers<-do.call(dplyr::bind_rows,topLayers)
  
  ept_add_secondary(myStudy,topLayers,'toplayers')
  
  ept_extractPixelStatistics(x = myStudy,
                             geometry = 'toplayers',
                             name = 'toplayers_mean',
                             raster = 'raw',
                             fun = 'mean')
  
  ept_extractPixelStatistics(x = myStudy,
                             geometry = 'toplayers',
                             name = 'toplayers_median',
                             raster = 'raw',
                             fun = 'median')
  
  clnmNms<-colnames(myStudy$currentAnalysis$exprs@exprs$toplayers_mean)[2:31]
  
  ept_apply(x = myStudy,
            name = 'toplayers_mean',
            newName = 'toplayers_mean_transformed',
            columns = clnmNms,
            fun = scales::modulus_trans(0)$transform)
  
  ept_apply(x = myStudy,
            name = 'toplayers_median',
            newName = 'toplayers_median_transformed',
            columns = clnmNms,
            fun = scales::modulus_trans(0)$transform)
  
  sample_encoded<-setNames(c('40040006','70060002','50020007','40020008','70010006','40020013','70020001','60050001','50020005'),
                           c('s0394','s0395','s0396','s0397','s0398','s0399','s0400','s0401','s0402'))
  
  fps<-setNames(c('3.1','5.2','8.2','11.0','6.5','11.2','1.6','2.8','1.4'),
                c('s0394','s0395','s0396','s0397','s0398','s0399','s0400','s0401','s0402'))
  
  IMC_text<-studyTable(myStudy)[,'IMC_text_file',drop=T]
  
  IMC_code<-substr(IMC_text,1,5)
  
  IMC_ROI<-regexpr('._(ROI_[0-9]{3}).',IMC_text,perl = T)
  IMC_ROI<-substr(IMC_text,attr(IMC_ROI,'capture.start'),attr(IMC_ROI,'capture.start')+attr(IMC_ROI,'capture.length')-1)
  
  st_samples(myStudy)<-sample_encoded[IMC_code]
  st_rois(myStudy)<-IMC_ROI
  st_bioGroups(myStudy)<-fps[IMC_code]
  
  archive(myStudy)
  
  rootFolder<-"C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4"
  targetFolder<-file.path(rootFolder,'Check_img')
  dir.create(targetFolder,showWarnings = F)
  
  uids<-st_uids(myStudy)
    
  stTbl<-studyTable(myStudy)
  
    for (ii in uids){
      
      rst<-myStudy$currentAnalysis$classification[[ii]]$label_clean
      ext<-raster::extent(rst)
      
      tiff(filename = file.path(targetFolder,
                                paste0(stTbl[ii,'IMC_text_file'],'.tiff')),
           width = (ext[2]-ext[1])*2,
           height = (ext[4]-ext[3])*2,
           units = 'px',
           compression = 'lzw',
           antialias = NULL)
      raster::plot(rst,xaxt='n',yaxt='n',legend=F,col=c('white',RColorBrewer::brewer.pal(8,'Set3')))
      TEMP_rst<-myStudy$currentAnalysis$exprs@exprs$toplayers$geom[myStudy$currentAnalysis$exprs@exprs$toplayers$uid==ii]
      plot(TEMP_rst,lwd=0.1,add=T)
      dev.off()
    
  }
  
  targetFolder<-file.path(rootFolder,'Check_img_nuc')
  dir.create(targetFolder,showWarnings = F)
  
  for ( ii in uids){
    
      nuc<-RUNIMCTEMP::quantNorm(myStudy$raster[[ii]]$x191ir.dna1.ir191di.,0.95)
      col<-RUNIMCTEMP::quantNorm(myStudy$raster[[ii]]$x158gd.ecadherin.gd158di.,0.95)
      vim<-RUNIMCTEMP::quantNorm(myStudy$raster[[ii]]$x143nd.vimentin.nd143di.,0.95)
      newrst<-raster::stack(list(nuc,col,vim))
      ext<-raster::extent(newrst)
      
      tiff(filename = file.path(targetFolder,
                                paste0(stTbl[ii,'IMC_text_file'],'.tiff')),
           width = (ext[2]-ext[1])*2,
           height = (ext[4]-ext[3])*2,
           units = 'px',
           compression = 'lzw',
           antialias = NULL)
      raster::plotRGB(newrst,1,3,2,scale=1)
      TEMP_rst<-myStudy$currentAnalysis$exprs@exprs$toplayers$geom[myStudy$currentAnalysis$exprs@exprs$toplayers$uid==ii]
      plot(TEMP_rst,lwd=0.1,add=T,col=NA,border='cyan')
      dev.off()
  }
  
  
  fileToUnlink<-list.files(
    file.path(myStudy$currentAnalysis$folder,
              'Temp'),
    full.names = T
  )
  
  unlink(fileToUnlink)
  
})

stopTime<-Sys.time()

print(stopTime-startTime)
bot <- telegram::TGBot$new(token = telegram::bot_token('R_lgd_bot'))
bot$set_default_chat_id(telegram::user_id('me'))

if (inherits(out,'try-error')){
  bot$sendMessage(paste(myStudy$name,' start: ',startTime,'end: ',stopTime,'...Failed'))
} else {
  bot$sendMessage(paste(myStudy$name,' start: ',startTime,'end: ',stopTime,'...Done'))
}