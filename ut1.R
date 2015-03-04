  #  ModelFilename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")

    #load the calibration fingerprints with their associated lab data
   # originalmyData<-loadDataFile(CalibrationFingerPrintsPath,"OriginalBrittanyOutliersRemoved.csv")   

    #returns a list with objects  $realTime $ChemData and $fingerPrints 


#startDates<-c(as.POSIXct(paste("2011-10-18 00:00:00",sep=""),tz="UTC"),
#              as.POSIXct(paste("2012-3-2 00:00:00",sep=""),tz="UTC")
#              )
#stopDates<-c(as.POSIXct(paste("2012-2-8 00:00:00",sep=""),tz="UTC"),
#             as.POSIXct(paste("2012-06-01 00:00:00",sep=""),tz="UTC")
#             )



startDates<-c(as.POSIXct(paste("2012-10-18 00:00:00",sep=""),tz="UTC"),
             as.POSIXct(paste("2013-3-2 00:00:00",sep=""),tz="UTC")
             )
stopDates<-c(as.POSIXct(paste("2013-3-2 00:00:00",sep=""),tz="UTC"),
            as.POSIXct(paste("2013-07-01 00:00:00",sep=""),tz="UTC")
            )


startDate<-startDates[y]
stopDate<-stopDates[y]

flow<-subset(Flow,"flow",startDates,stopDates)
chemical<-1
numComp<-5
k<-1
      fileType<-Types[k]   #specify file TYPE 1-REGFP 2-1STDER 3-TURBCOMP 4-1STDERTURBCOMP 
      subsetRate=0
  
      calibration<-subsetSpecData(fileType,"calibration",startDates,stopDates,chemN[chemical])
      specDataToModel<-subsetSpecData(fileType,"fingerPrints",startDates,stopDates)
      
      modelOutput<-modelExecution(chemical,numComp,subsetRate,calibration,specDataToModel,fitEval,fitFile,fitFileOut)
        
       
          plot(as.POSIXct(flow$realTime,origin="1970-01-01 00:00:00"),
               flow$flow,col="black",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
               xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1])),
          )
            
          par(new=TRUE)
  
          plot(specDataToModel$realTime,rowSums(specDataToModel$fingerPrints),type="l",col="magenta")
          
          par(new=TRUE)
          plotMax<-max(modelOutput$ObservedandPredicted[,1])*1.5
          plotMin<-min((min(modelOutput$ObservedandPredicted[,1])-min(modelOutput$ObservedandPredicted[,1])*0.1),0)
          plot(as.POSIXct(modelOutput$PredictedConcentrations[,1],origin="1970-01-01 00:00:00"),
               modelOutput$PredictedConcentrations[,2],type="l",col=rainbow(10)[i],
               main=chemicals[chemical],
               ylim=c(plotMin,plotMax),
               xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE)),
               xlab="date",
               ylab=paste(chemicals[chemical],"concentration",sep=" ")
          )
          points(modelOutput$ObservedandPredicted[,3],
                 modelOutput$ObservedandPredicted[,1],
                 col="green",cex=0.8,pch=18)
          
          points(modelOutput$ObservedandPredicted[,3],
                 modelOutput$ObservedandPredicted[,2],
                 col="blue",cex=.8,pch=20)
          abline(h=0)
          
          mtext(paste(fileType,numComp,startDate,sep="   "))
                
