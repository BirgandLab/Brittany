#library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package
library(xts)
library(hydroTSM)
#library(signal)
library(tseries)

set.seed(1) #should have some repeatability, if we need


runBSModel<-function(chem,numComp,dates,calibData,modelData,type,flow,subsetRatio=0,bsIterations=10,keepSpecData,fitEval,fitFile,fitFileOut){
    #subset data to the windows specified above
    lotOfData<-subsetAll(calibData,modelData,flow,chem,cbind(startDates,stopDates),keepSpecData)
    #run model
    indices<-bsIndMat(lotOfData$calibData$realTime,bsIterations)
        #rows have realTime data for each "simulation"  
    for (r in 1:bsIterations){
          #reassemble calib data
            c<-as.matrix(indices[r,])
            colnames(c)<-"realTime"
          #merge index, with calibDataset, keeping index
              #make a data frame of the calibDataData
              cD<-cbind(lotOfData$calibData$realTime,lotOfData$calibData$ChemData,lotOfData$calibData$fingerPrints)
              colnames(cD)[1:2]<-c("realTime",chem)
              
            cc<-merge(c,cD,by.x="realTime",by.y="realTime",all.x=TRUE,sort=FALSE)
          #restructure it to resemble the data list I have been using
            bsRun<-list(realTime=cc[,1],ChemData=cc[,2],fingerPrints=cc[,3:NCOL(cc)])
          
          #then feed it into the model function
            
            
            putout<-modelExecution(numComp,subsetRatio,bsRun,lotOfData$specData,fitEval,fitFile,fitFileOut)
            
            if(r==1){
                      modelQuality<-as.data.frame(matrix(nrow=bsIterations,ncol=dim(putout$modelQuality)[2]+5))
                      modelQuality[1,]<-c(putout$modelQuality,chem,numComp,r,subsetRatio,type)
                      colnames(modelQuality)<-c(colnames(putout$modelQuality),"chem","numComp","iter","subset","type")
                      
                      predictions<-as.data.frame(matrix(nrow=dim(putout$PredictedConcentrations)[1],ncol=bsIterations+1))
                      predictions[,1:2]<-putout$PredictedConcentrations
                      
                      Observed<-as.data.frame(lotOfData$calibData$realTime)#all possible observations
                      colnames(Observed)<-"realTime"
                      
                      Observed<-merge(Observed,putout$ObservedandPredicted,by.x="realTime",by.y="realTime",all.x=TRUE)
                      
                    }else{
                      modelQuality[r,]<-c(putout$modelQuality,chem,numComp,r,subsetRatio,type)
                      predictions[,r+1]<-putout$PredictedConcentrations
                      Observed<-merge(Observed,putout$ObservedandPredicted,by.x="realTime",by.y="realTime",all.x=TRUE)
                    }
  }
  return(list(modelQuality=modelQuality,predictions=predictions,calibPoints=Observed))  
}



#set environment
#*******************************************************************************************************************/
fitPath<-"C:/Users/FBlab/Desktop/FITEVAL2_win/FITEVAL2_win/"    #location on downstairs box
OutputPath<-"C:/Users/FBlab/Desktop/work_here/Data/Brittany/Output/"   #downstairs box
flowPath<-"C:/Users/FBlab/Desktop/work_here/Data/Brittany/Flow/"
FingerPrintPath<-"C:/Users/FBlab/Desktop/work_here/Data/Brittany/FingerPrints/" #downstairs box
CalibrationFingerPrintsPath<-"C:/Users/FBlab/Desktop/work_here/Data/Brittany/CalibrationFingerPrints/"


fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
filePaths=cbind(fitEval,fitFile,fitFileOut)

ModelFilename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")
excludeRows<-c(29,153,155,471,482,486,514,525,530,548,615,616,617,652,658,660,668,671,673,735,790,817,831,845,902,911)
if(!exists("britFuncLoaded",1)){
  source("C:/Users/FBlab/Desktop/Brittany/BrittanyFunctions.R", echo=FALSE)
  source("C:/Users/FBlab/Desktop/Brittany/bootstrapIndexMatrix.R", echo=FALSE)
}
if(!exists("flow",1)){
  source("C:/Users/FBlab/Desktop/Brittany/loadData.R", echo=FALSE)
}

#for year 2010-2011
startDates<-c(as.POSIXct(paste("2010-11-18 00:00:00",sep=""),tz="UTC"),
              as.POSIXct(paste("2011-04-25 00:00:00",sep=""),tz="UTC")
)
stopDates<-c(as.POSIXct(paste("2011-03-01 00:00:00",sep=""),tz="UTC"),
             as.POSIXct(paste("2011-05-20 00:00:00",sep=""),tz="UTC")
)

dates<-c(min(startDates),max(stopDates))

#bsIndMat(time,reps)
NO3<-runBSModel(chem="NNO3",
               numComp=5,
               dates=dates,
               calibData=prunedTC1D,
               modelData=Derivative,
               type="turbComp",
               flow=Flow,
               subset=1,
               bsIterations=15,
               keepSpecData=1,
               fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)



ak<-runBSModel(chem="DIC",
            numComp=5,
            dates=dates,
            calibData=prunedTC1D,
            modelData=Derivative,
            type="turbComp",
            flow=Flow,
            subset=1,
            bsIterations=15,
            keepSpecData=1,
            fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)



colurs<-rainbow(10)

b<-ak
#b<-a
#bootstrap a model to look for error in predictions based on calibration with existing full dataset
rnge<-c(min(b$predictions),max(b$predictions[,2])*1.2)
plot(b$predictions[,1],b$predictions[,2],col=colurs[1],ylim=rnge,type="l")
for(i in 2:15){
  points(b$predictions[,i],col=colurs[i],type="l")
  
}