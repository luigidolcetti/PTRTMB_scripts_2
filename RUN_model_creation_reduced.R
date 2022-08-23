startTime<-Sys.time()
library(RUNIMCTEMP)

out<-try({
TEMP<-list.files("D:/Luigi/DAICHI_FULL_SET")
TEMPNUMBERS<-c(6,23,27,30,33,34,45,47,48,71,74,86)
# TEMP[TEMPNUMBERS]

myStudy<-initStudy(fn_studyName = 'RUN_model_creation_reduced',
                   fn_rootFolder = "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_4",
                   fn_rawDataFolder = "D:/Luigi/DAICHI_FULL_SET",
                   fn_whichFiles = TEMPNUMBERS,
                   fn_whichColumns = 'named')

newAnalysis(myStudy,'Anal_001')

archive(myStudy)

# runIMC(myStudy)

addFilter(x = myStudy,
          filter = 'vanvliet',
          parameters = list(list(sigma=1,order=0,axis='x'),
                            list(sigma=3,order=0,axis='y')),
          channels = ch_Rnames(myStudy),
          append = F)
addFilter(x = myStudy,
          filter = 'vanvliet',
          parameters = list(list(sigma=1,order=1,axis='y'),
                            list(sigma=3,order=1,axis='x')),
          channels = ch_Rnames(myStudy),
          append = T)
addFilter(x = myStudy,
          filter = 'vanvliet',
          parameters = list(list(sigma=1,order=2,axis='x'),
                            list(sigma=3,order=2,axis='y')),
          channels = ch_Rnames(myStudy),
          append = T)

addFilter(x = myStudy,
          filter = 'deriche',
          parameters = list(list(sigma=1,order=0,axis='y'),
                            list(sigma=3,order=0,axis='x')),
          channels = ch_Rnames(myStudy),
          append = T)
addFilter(x = myStudy,
          filter = 'deriche',
          parameters = list(list(sigma=1,order=1,axis='x'),
                            list(sigma=3,order=1,axis='y')),
          channels = ch_Rnames(myStudy),
          append = T)
addFilter(x = myStudy,
          filter = 'deriche',
          parameters = list(list(sigma=1,order=2,axis='y'),
                            list(sigma=3,order=2,axis='x')),
          channels = ch_Rnames(myStudy),
          append = T)

addFilter(myStudy,
          'blur_anisotropic',
          list(list(amplitude=3),list(amplitude=1)),
          channels = ch_Rnames(myStudy),append=T)

deployFilters(x = myStudy,cl = 7)

addExtractionDirectives(myStudy,c(1,0),'core',append=F)

polyList<-list.files(
  "C:/Users/luigi/Documents/R_projects/ALL_POSSIBLE_DAICHI_3/backup_polygons/polygons",
  full.names = T
)

file.copy(from = polyList,
          to = file.path(myStudy$currentAnalysis$folder,'training/polygons/'))

extractTrainingFeatures(myStudy)

availableTF<-tf_featureList(myStudy)
availableLabels<-tf_labelList(myStudy)
availabeleNumbers<-table(myStudy$currentAnalysis$trainingFeatures@value$label)

eventNumbers<-rep(3000,length(availableLabels))
names(eventNumbers)<-availableLabels

addClassificationDirectives(x = myStudy,
                            method = 'randomForest',
                            methodParameters = list(responseVariable = 'label',
                                                    predictiveFeatures = availableTF,
                                                    PvalueTreshold = 0.25,
                                                    eventNumbers = eventNumbers,
                                                    ntree=1000,
                                                    mtry=25,
                                                    importance=F,
                                                    nodesize=1,
                                                    do.trace=500))

makeClassificationModel(myStudy,
                        method = 'randomForest')

eventNumbers<-c(3200,8900,7400,21000,17000,7100,4900)
names(eventNumbers)<-tf_parLabelList(myStudy)[2:8]

addClassificationDirectives(x = myStudy,
                            method = 'randomOnions',
                            methodParameters = list(responseVariable = 'SLI',
                                                    classificationLyr = 'label_clean',
                                                    predictiveFeatures = availableTF,
                                                    prefix = 'sko_',
                                                    eventNumbers=eventNumbers,
                                                    labels = tf_parLabelList(myStudy)[c(2:8)],
                                                    importance = F,
                                                    mtry = 150,
                                                    ntree = 1000,
                                                    do.trace=500))

makeClassificationModel(x=myStudy,method = 'randomOnions')

archive(myStudy)

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