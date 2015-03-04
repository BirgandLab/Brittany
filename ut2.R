prunedO<-loadDataFile(CalibrationFingerPrintsPath,"OriginalBrittanyP.csv")   

calibration2<-subsetSpecData("prunedO","calibration",startDates,stopDates,chemN[chemical])


calibration<-subsetSpecData("original","calibration",startDates,stopDates,chemN[chemical])


startDates<-c(as.POSIXct(paste("2011-10-18 00:00:00",sep=""),tz="UTC"),
              as.POSIXct(paste("2012-3-2 00:00:00",sep=""),tz="UTC")
)
stopDates<-c(as.POSIXct(paste("2012-2-8 00:00:00",sep=""),tz="UTC"),
             as.POSIXct(paste("2012-06-01 00:00:00",sep=""),tz="UTC")
)

startDate<-min(startDates)
stopDate<-max(stopDates)

for(chemical in 1:9){
  for(numComp in 3:9){
    #numComp<-4    #specify number of components to use in the model
    #for(k in 1:length(Types)){
     # fileType<-Types[1]   #specify file TYPE 1-REGFP 2-1STDER 3-TURBCOMP 4-1STDERTURBCOMP 
      subsetRate=0
      #   #mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)
      iD[counter,1]<-paste((as.POSIXlt(startDate))$year+1900,"-",(as.POSIXlt(stopDate))$year+1900,sep="")
      iD[counter,2]<-fileType
      iD[counter,3]<-numComp
      iD[counter,4]<-subsetRate
      iD[counter,30]<-Chem[chemN[chemical]]
      #change here to send vector of startDates and stopDate
      
     
      specDataToModel<-subsetSpecData(fileType,"fingerPrints",startDates,stopDates)
      
      #calibration<-subsetSpecData("fileType","calibration",startDates,stopDates,chemN[chemical])
            #do the subsetHere  
                      logicalIndexofInclusion<-matrix(nrow=length(modifiedmyData$realTime),ncol=length(startDates))
                      
                      for(i in 1:length(startDates)){
                        logicalIndexofInclusion[,i]<-(modifiedmyData$realTime>startDates[i]&modifiedmyData$realTime<stopDates[i])
                      }
                      
                      keep<-as.logical(rowSums(logicalIndexofInclusion,na.rm=TRUE))
                      
                      
                        realTime<-modifiedmyData$realTime[keep]
                        fingerprint<-modifiedmyData$fingerPrint[keep,]
                        ChemData<-modifiedmyData$ChemData[keep,]
                        
                        fp<-cbind(realTime,fingerprint,as.matrix(ChemData[,chemN[chemical]]))
                        goodData<-fp[complete.cases(fp[,2:dim(fp)[2]]),] #removes all the rows for which there is a NA--keeping time in there
                        realTime<-goodData[,1]                      #pull components back out
                        ChemData<-goodData[,dim(goodData)[2]]
                        fingerprints<-goodData[,-1] 
                        fingerprints<-fingerprints[,-dim(fingerprints)[2]]
                        
                        
                        
                        calibMod<-(list(realTime=realTime,fingerPrints=fingerprints,ChemData=ChemData))
                                          
          
      
      
      for(i in 1:1){
        #returns a list
        #   output$modelQuality a named vector of summary characteristics for the model and its verification, 
        #     including NSE, r2 and fitEval for the model and, if applicable validation datasets
        #   output$PredictedConcentrations
       # modelOutput<-modelExecution(chemical,numComp,subsetRate,calibration,specDataToModel,fitEval,fitFile,fitFileOut)
        modelOutput<-modelExecution(chemical,numComp,subsetRate,calibMod,specDataToModel,fitEval,fitFile,fitFileOut)
        #          totmin<-max((modelOutput$PredictedConcentration[,1])-min(modelOutput$PredictedConcentration[,1]))/60
        #         PredictedAnnualConcentrationTS<-as.data.frame(approx(modelOutput$PredictedConcentration[,1],
        #                              modelOutput$PredictedConcentration[,2],n=totmin+1)) #create a time series of [ ] at 1 minute intervals for the whole dataset
        #         colnames(PredictedAnnualConcentrationTS)<-c("realTime","concentration")                    
        # PredictedAnnualConcentrationTS$y[(PredictedAnnualConcentratinTS$y<0)]<-0
        #         PredictedEventConcentrationTS<-(PredictedAnnualConcentrationTS$realTime>as.numeric(eventStart)&
        #                                           PredictedAnnualConcentrationTS$realTIme<as.numeric(eventStop))
        
        iD[counter,5:30]<-modelOutput$modelQuality
        
        #if(i==1){
        png(file = paste(OutputPath,
                         Chem[chemN[chemical]],
                         "_",
                         fileType,
                         "_pruned_",
                         numComp,
                         "_",
                         (as.POSIXlt(startDate))$year+1900,
                         ".png",sep=""),
            width=1000,height=500
            
        )
        
        if (sum(is.na(iD[counter,5:30]))!=26) {
          plot(as.POSIXct(flow$realTime,origin="1970-01-01 00:00:00"),
               flow$flow,col="black",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
               xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1])),
          )
          
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
        }
        
        dev.off()
        #  }
        #                                       else{
        #                                                   points(as.POSIXct(PredictedAnnualConcentrationTS$realTime,origin="1970-01-01 00:00:00"),
        #                                                           PredictedAnnualConcentrationTS$concentration,type="l",col=rainbow(10)[i]
        #                                                           )
        #                                             }
        
        
        counter<-counter+1
      }#end #repititions
    }#end for each fileTYpe
  }#end for #componenets
}#end for each chemical