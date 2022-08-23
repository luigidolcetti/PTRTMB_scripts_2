rootFolder<-"C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4"
analFolder<-"POST_segmentation_tables"

dir.create(file.path(rootFolder,analFolder))

dirList<-list.dirs(rootFolder,recursive = F)
whichDir<-grepl('RUN_s',dirList,fixed = T)
dirList<-dirList[whichDir]
fileList<-paste0(dirList,"/archive/studyTable.RDS")

studyTable<-lapply(fileList,function(fl){
  out<-readRDS(fl)
  return(out)
})

studyTable<-do.call(rbind.data.frame,studyTable)

write.csv(studyTable,file.path(rootFolder,analFolder,'studyTable.csv'))
