#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#
densityDependentSubset<-function(ChemConc,realTime,fingerPrint,subsetRatio,Replace){
              #ChemConc a vector of chemical concentrations for each fingerprint
              #fingerPrint uv/Vis values
              #realTime the time at which samples were taken
              #subsetRatio is the proportion of the total samples to use
              #replace is a logical value (TRUE or FALSE) specifying whether when subsampling is done with replacement or not
              
            #FP<-cbind(myData$fingerPrints,myData$realTime,myData$chemConc)
              
              FP<-cbind(fingerPrint,realTime,ChemConc)
              #SORTED IN ASCENDING ORDER BY CHemConc
              FP<-FP[order(FP[,dim(FP)[2]]),] 
              ChemConc<-FP[,dim(FP)[2]]
               
              Histogram<-hist(ChemConc,breaks="FD",plot=FALSE)
              #histSturges<-hist(ChemConc,breaks="Sturges")
              #histScott<-hist(ChemConc,breaks="Scott",plot=FALSE)
              #histFD<-hist(ChemConc,breaks="FD",plot=FALSE)
              
              sampleRatio<-subsetRatio
              count<-0
              
              #for each of the bins in the histogram
              for (q in seq(1,(length(Histogram$breaks)-1),1)){  
                    #find the fingerprints that fit the historgram bin
                      candidates<-FP[(ChemConc>=Histogram$breaks[q])&(ChemConc<=Histogram$breaks[q+1]),]  
                    #calculate the density dependent number of samples from this bin  
                      sampleNum<-sum(ceiling(Histogram$counts[q]*sampleRatio))   
                    #pick the samples from the candidates
                      selected<-candidates[sample(1:dim(candidates)[1],sampleNum,replace=Replace),]   
                    #if it is the first time through, make a container variable to accumulate the calibration data
                      if(count==0){
                            calibration<-selected
                          }
                    #if it isn't the first time through, just add the data to the end
                      if(count>0){
                            calibration<-rbind(calibration,selected)
                          }
                    #augment counter for some reason
                      count<-count+1
                }#end for q in sequence
              
              
               # return (calibration$RT calibration$ChemConc calibration$FP)
              calibrationRealTime<-calibration[,(dim(calibration)[2]-1)]  #realtime is in the second to last column
              calibrationChemConc<-calibration[,dim(calibration)[2]]      #chemConc is in the last column
              calibration<-calibration[,-dim(calibration)[2]]             #remove the last column
              calibration<-calibration[,-dim(calibration)[2]]           #remove the last column again
              
              return(list(realTime=calibrationRealTime,ChemConc=calibrationChemConc,fingerPrint=calibration))
              
              
            } #end density dependent subset




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

loadDataFile<-function(filepath,filename){
              data<-read.table(file=paste(filepath,filename,sep=""),sep=",",header=TRUE,skip=0)
              
              #parse the date
              Date<-substr(data$DateScan, 1, 10)        #the date
              T<-substr(data$DateScan, 12, 19)       #the time
              T[T==""]="00:00:00"                    #ah, but where it should be midnight, often it was " "
              D<-paste(Date,T,sep=" ")               #put back in the fixed values
              D<-strptime(D, '%d/%m/%Y %H:%M:%S',tz="UTC")    #convert it to a _flawless_ time opject
              data$DateScan<-D                       #store it back in the matrix
              
              realTime<-(data$DateScan)                                       #pull out time
              ChemData<-data[,(dim(data)[2]-16):(dim(data)[2])]            #pull out the analyte data
              fingerPrints<-data[,-(dim(data)[2]-20):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
              fingerPrints<-fingerPrints[,-1:-2]                      #remove status and datetime
              
              return(list(realTime=realTime, ChemData=ChemData, fingerPrints=fingerPrints))
            }





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#a function that you pass 2 matrices 1 with scan fingerprints. the other with chem conc data
# the function will remove datalines with NANs and incomplete records
#the function will run a PLSR model on the data
#so we need to specify number of points to use
#type of subsampling  uniform, density dependent, windowed or not
#number of points to use to build the model
#it will return model parameters, a matrix of measured and predicted values, and a vector of gooness of fit parameters calculated from the matlab routine
numberOfComponentsToUse<-function(fingerPrint,ChemConc,ylabel,subTitle){
              #fingerPrint<-allFingerPrints
            
              fingerPrint<-cbind(fingerPrint,as.matrix(ChemConc))                          #bind the two matrices together to determine complete cases
              fingerPrint<-fingerPrint[complete.cases(fingerPrint[,2:dim(fingerPrint)[2]]),] # removes all the rows for which there is a NA--keeping time in there
              
              ChemConc<-as.matrix(fingerPrint[,dim(fingerPrint)[2]])
              if (length(ChemConc)<40) {
                    maxNCompsToTest<-length(ChemConc)*0.8;
                    }else{maxNCompsToTest<-30}
              fingerPrint<-fingerPrint[,-dim(fingerPrint)[2]]
              
              Fat<-RMSEP(plsr(ChemConc~data.matrix(fingerPrint),ncomp=maxNCompsToTest,validation="CV"))
              c<-as.matrix(Fat$val[2,1,])
              ex<-1:(length(c))
              nComps<-min(which(c==min(c)))
              plot(ex,c,xlab=c("model components"),ylab=ylabel,main=(""),type="b")
              mtext(subTitle)
  
  return(nComps=nComps)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
PLSRFitAndTest<-function(fingerPrint,ChemConc,realTime,numParameters,fitEval,fitFile,fitFileOut,subsample) {
  
  #assemble so complete cases keeps time data organized
  fingerPrint<-cbind(realTime,fingerPrint,as.matrix(ChemConc))                          #bind the two matrices together to determine complete cases
  fingerPrint<-fingerPrint[complete.cases(fingerPrint[,2:dim(fingerPrint)[2]]),] # removes all the rows for which there is a NA--keeping time in there
 

  
  #pull apart again
  ChemConc<-as.matrix(fingerPrint[,dim(fingerPrint)[2]])                       #strip off the chemConc again
  
    if(subsample==0){
      Stats<-matrix(nrow=1,ncol=5)
      colnames(Stats)<-c("n","r2","rmse","nrmse","slope")
    }
  
    if(subsample > 0){#enter density dependent subsampling with Subsample Value guiding the process
          #How about using native HIST function with Sturgess, scott and DF methods for breaking.
          #Then, if possible take 1 point from each bin
          #histKyle<-hist(ChemConc,breaks=steps, plot=FALSE)
          #histSturges<-hist(ChemConc,breaks="Sturges")
          #histScott<-hist(ChemConc,breaks="Scott",plot=FALSE)
          Stats<-matrix(nrow=1,ncol=15)
          colnames(Stats)<-c("n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope")
          Histogram<-hist(ChemConc,breaks="FD",plot=FALSE)
            
          sampleRatio<-subsample
          count<-0
         # fingerPrint1<-fingerPrint[order(fingerPrint[,dim(fingerPrint)[2]]),] #SORTED IN ASCENDING ORDER
          for (q in seq(1,(length(Histogram$breaks)-1),1)){
            #for (q in seq((length(Histogram$breaks)-1)*0.25,(length(Histogram$breaks)-1),1)){
            inTheBin<-fingerPrint[(ChemConc>=Histogram$breaks[q])&(ChemConc<=Histogram$breaks[q+1]),]
            sampleNum<-sum(ceiling(Histogram$counts[q]*sampleRatio))
            someNumber<-inTheBin[sample(1:dim(inTheBin)[1],sampleNum,replace=TRUE),]
            if(count==0){
              calibration<-someNumber
            }
            if(count>0){
              calibration<-rbind(calibration,someNumber)
            }
            count<-count+1
          }
          remove(someNumber)
      
  
            calibrationChemConc<-(calibration[,dim(calibration)[2]])
            calibrationRealTime<-as.matrix(calibration[,1]) 
            calibrationFingerPrint<-calibration[,-1]
            calibrationFingerPrint<-calibrationFingerPrint[,-dim(calibrationFingerPrint)[2]]
            
            RealTime<-as.matrix(fingerPrint[,1])                        #complete cases realtime
            fingerPrint<-fingerPrint[,-dim(fingerPrint)[2]]  #complete cases fingerprints (remove chemConc)
            fingerPrint<-fingerPrint[,-1]                    #complete cases fingerprints (remove RealTime)
            
            #RUN PLSR MODEL FOR calibration SUBSAMPLE
            calibrationFit<-plsr(calibrationChemConc~data.matrix(calibrationFingerPrint),ncomp=numParameters,validation="CV")
          
            #predict concentrations for the calibration (uninteresting)
            calibPredict1<-predict(calibrationFit,data.matrix(calibrationFingerPrint),ncomp=numParameters,type=c("response"))
                fitQualityModel<-OB(calibrationChemConc,calibPredict1,fitEval,fitFile,fitFileOut)
                op1<-cbind(calibrationChemConc,calibPredict1,calibrationRealTime)
                fitStats<-summary(lm(calibrationChemConc~calibPredict1))
                confInt<-predict(lm(calibrationChemConc~calibPredict1),interval='prediction',level=0.95)
                Stats[1]<-fitStats$df[2]
                Stats[2]<-fitStats$r.squared
                Stats[3]<-rmse(obs=as.matrix(calibrationChemConc),sim=as.matrix(calibPredict1))
                Stats[4]<-nrmse(obs=as.matrix(calibrationChemConc),sim=as.matrix(calibPredict1))
                Stats[5]<-fitStats$coefficients[2,1] #slope
            #predict concentrations for all available lab data
            calibPredict2<-predict(calibrationFit,data.matrix(fingerPrint),ncomp=numParameters,type=c("response"))
                fitQualityModelFull<-OB(ChemConc,calibPredict2,fitEval,fitFile,fitFileOut)
                op2<-cbind(ChemConc,calibPredict2,RealTime)
            
            fitStats<-summary(lm(ChemConc~calibPredict2))
            confInt<-predict(lm(ChemConc~calibPredict2),interval='prediction',level=0.95)
            Stats[6]<-fitStats$df[2]
            Stats[7]<-fitStats$r.squared
            Stats[8]<-rmse(obs=as.matrix(ChemConc),sim=as.matrix(calibPredict2))
            Stats[9]<-nrmse(obs=as.matrix(ChemConc),sim=as.matrix(calibPredict2))
            Stats[10]<-fitStats$coefficients[2,1] #slope
     }#end if subsample>0
  RealTime<-as.matrix(fingerPrint[,1])                        #complete cases realtime
  fingerPrint<-fingerPrint[,-dim(fingerPrint)[2]]  #complete cases fingerprints (remove chemConc)
  fingerPrint<-fingerPrint[,-1]                    #complete cases fingerprints (remove RealTime)
  
  
#run PLSR Model for all available data
 Fit<-plsr(ChemConc~data.matrix(fingerPrint),ncomp=numParameters,validation="CV")  #PLSR model to predict chemConc with cross validation
  #make predictions 
 Predict<-predict(Fit,data.matrix(fingerPrint),ncomp=numParameters,type=c("response"))
  #assess quality    
 fitQualityFull<-OB(ChemConc,Predict,fitEval,fitFile,fitFileOut)
      op3<-cbind(ChemConc,Predict,RealTime)
         if(subsample>0){
                       fitStats<-summary(lm(ChemConc~Predict))
                       confInt<-predict(lm(ChemConc~Predict),interval='prediction',level=0.95)
                       Stats[11]<-fitStats$df[2]
                       Stats[12]<-fitStats$r.squared
                       Stats[13]<-rmse(obs=ChemConc,sim=as.matrix(Predict))
                       Stats[14]<-nrmse(obs=ChemConc,sim=as.matrix(Predict))
                       Stats[15]<-fitStats$coefficients[2,1] #slope
                       }
         if(subsample==0){
                       fitStats<-summary(lm(ChemConc~Predict))
                       confInt<-predict(lm(ChemConc~Predict),interval='prediction',level=0.95)
                       Stats[1]<-fitStats$df[2]
                       Stats[2]<-fitStats$r.squared
                       Stats[3]<-rmse(obs=ChemConc,sim=as.matrix(Predict))
                       Stats[4]<-nrmse(obs=ChemConc,sim=as.matrix(Predict))
                       Stats[5]<-fitStats$coefficients[2,1] #slope
                     }
    
  
#if we did not subsample, return the full model result
if (subsample==0){
  return (list(Fit=Fit,
               fitQuality=fitQualityFull,
               ObservedAndPredicted=op3,
               Stats=Stats
               )
          )
   }
  if (subsample>0){
    return (list(Fit=Fit,
                 calibrationFit=calibrationFit,
                 fitQuality=c(fitQualityModel,fitQualityModelFull,fitQualityFull),
                 OaP1=op1,
                 OaP2=op2,
                 ObservedAndPredicted=op3,
                 Stats=Stats
                 )
            )
    
  }
 
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

OB<-function(observed,predicted,fitEval,fitFile,fitFileOut){
                ObsAndPred<-cbind(observed,predicted)
                write.table(ObsAndPred,file=fitFile, append = FALSE,row.names=FALSE,col.names=FALSE)
                system(paste(fitEval,fitFile,sep=" "),wait=TRUE,show.output.on.console=FALSE)
                
                a<-readLines(fitFileOut)
                #unlink(fitFileOut)
                
                #in line a[7] there are data for very good fits, after tha : and before the %
                fitQuality<-c(0,0,0,0)
                names(fitQuality)<-c("veryGoood","good","acceptable","bad")
                
                for (i in 7:10){
                  b<-toString(a[i])
                  b<-unlist(strsplit(b,":"))[2]
                  fitQuality[i-6]<-as.numeric(unlist(strsplit(b,"%"))[1])
                  #print(fitQuality[i-6])
                }
                
return(fitQuality)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
PLSRFitAndTestNoNSE<-function(fingerPrint,ChemConc,realTime,numParameters,fitEval,fitFile,fitFileOut,subsample) {
  
  #assemble so complete cases keeps time data organized
  fingerPrint<-cbind(realTime,fingerPrint,as.matrix(ChemConc))                          #bind the two matrices together to determine complete cases
  fingerPrint<-fingerPrint[complete.cases(fingerPrint[,2:dim(fingerPrint)[2]]),]        #removes all the rows for which there is a NA--keeping time in there
  
  
  
  #pull apart again
  ChemConc<-as.matrix(fingerPrint[,dim(fingerPrint)[2]])                       #strip off the chemConc again
  
  if(subsample==0){
    Stats<-matrix(nrow=1,ncol=5)
    colnames(Stats)<-c("n","r2","rmse","nrmse","slope")
  }
  
  if(subsample > 0){#enter density dependent subsampling with Subsample Value guiding the process
    #How about using native HIST function with Sturgess, scott and DF methods for breaking.
    #Then, if possible take 1 point from each bin
    #histKyle<-hist(ChemConc,breaks=steps, plot=FALSE)
    #histSturges<-hist(ChemConc,breaks="Sturges")
    #histScott<-hist(ChemConc,breaks="Scott",plot=FALSE)
    Stats<-matrix(nrow=1,ncol=15)
    colnames(Stats)<-c("n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope")
    Histogram<-hist(ChemConc,breaks="FD",plot=FALSE)
    
    sampleRatio<-subsample
    count<-0
    # fingerPrint1<-fingerPrint[order(fingerPrint[,dim(fingerPrint)[2]]),] #SORTED IN ASCENDING ORDER
    for (q in seq(1,(length(Histogram$breaks)-1),1)){
      #for (q in seq((length(Histogram$breaks)-1)*0.25,(length(Histogram$breaks)-1),1)){
      inTheBin<-fingerPrint[(ChemConc>=Histogram$breaks[q])&(ChemConc<=Histogram$breaks[q+1]),]
      sampleNum<-sum(ceiling(Histogram$counts[q]*sampleRatio))
      someNumber<-inTheBin[sample(1:dim(inTheBin)[1],sampleNum,replace=TRUE),]
      if(count==0){
        calibration<-someNumber
      }
      if(count>0){
        calibration<-rbind(calibration,someNumber)
      }
      count<-count+1
    }
    remove(someNumber)
    
    
    calibrationChemConc<-(calibration[,dim(calibration)[2]])
    calibrationRealTime<-as.matrix(calibration[,1]) 
    calibrationFingerPrint<-calibration[,-1]
    calibrationFingerPrint<-calibrationFingerPrint[,-dim(calibrationFingerPrint)[2]]
    
    RealTime<-as.matrix(fingerPrint[,1])                        #complete cases realtime
    fingerPrint<-fingerPrint[,-dim(fingerPrint)[2]]  #complete cases fingerprints (remove chemConc)
    fingerPrint<-fingerPrint[,-1]                    #complete cases fingerprints (remove RealTime)
    
    #RUN PLSR MODEL FOR calibration SUBSAMPLE
    calibrationFit<-plsr(calibrationChemConc~data.matrix(calibrationFingerPrint),
                         ncomp=numParameters,validation="CV")
    
    #predict concentrations for the calibration (uninteresting)
    calibPredict1<-predict(calibrationFit,data.matrix(calibrationFingerPrint),
                           ncomp=numParameters,type=c("response"))
        fitQualityModel<-OB(calibrationChemConc,calibPredict1,fitEval,fitFile,fitFileOut)
        op1<-cbind(calibrationChemConc,calibPredict1,calibrationRealTime)
        fitStats<-summary(lm(calibrationChemConc~calibPredict1))
        confInt<-predict(lm(calibrationChemConc~calibPredict1),interval='prediction',level=0.95)
    Stats[1]<-fitStats$df[2]
    Stats[2]<-fitStats$r.squared
    Stats[3]<-rmse(obs=as.matrix(calibrationChemConc),sim=as.matrix(calibPredict1))
    Stats[4]<-nrmse(obs=as.matrix(calibrationChemConc),sim=as.matrix(calibPredict1))
    Stats[5]<-fitStats$coefficients[2,1] #slope
    
    #predict concentrations for all available lab data
    calibPredict2<-predict(calibrationFit,data.matrix(fingerPrint),ncomp=numParameters,type=c("response"))
        fitQualityModelFull<-OB(ChemConc,calibPredict2,fitEval,fitFile,fitFileOut)
        op2<-cbind(ChemConc,calibPredict2,RealTime)
    
        fitStats<-summary(lm(ChemConc~calibPredict2))
        confInt<-predict(lm(ChemConc~calibPredict2),interval='prediction',level=0.95)
    Stats[6]<-fitStats$df[2]
    Stats[7]<-fitStats$r.squared
    Stats[8]<-rmse(obs=as.matrix(ChemConc),sim=as.matrix(calibPredict2))
    Stats[9]<-nrmse(obs=as.matrix(ChemConc),sim=as.matrix(calibPredict2))
    Stats[10]<-fitStats$coefficients[2,1] #slope
  }#end if subsample>0
  RealTime<-as.matrix(fingerPrint[,1])                        #complete cases realtime
  fingerPrint<-fingerPrint[,-dim(fingerPrint)[2]]  #complete cases fingerprints (remove chemConc)
  fingerPrint<-fingerPrint[,-1]                    #complete cases fingerprints (remove RealTime)
  
  
  #run PLSR Model for all available data
  Fit<-plsr(ChemConc~data.matrix(fingerPrint),ncomp=numParameters,validation="CV")  #PLSR model to predict chemConc with cross validation
  Predict<-predict(Fit,data.matrix(fingerPrint),ncomp=numParameters,type=c("response"))
  #fitQualityFull<-OB(ChemConc,Predict,fitEval,fitFile,fitFileOut)
  op3<-cbind(ChemConc,Predict,RealTime)
  
  if(subsample>0){
    fitStats<-summary(lm(ChemConc~Predict))
    confInt<-predict(lm(ChemConc~Predict),interval='prediction',level=0.95)
    Stats[11]<-fitStats$df[2]
    Stats[12]<-fitStats$r.squared
    Stats[13]<-rmse(obs=ChemConc,sim=as.matrix(Predict))
    Stats[14]<-nrmse(obs=ChemConc,sim=as.matrix(Predict))
    Stats[15]<-fitStats$coefficients[2,1] #slope
  }
  if(subsample==0){
    fitStats<-summary(lm(ChemConc~Predict))
    confInt<-predict(lm(ChemConc~Predict),interval='prediction',level=0.95)
    Stats[1]<-fitStats$df[2]
    Stats[2]<-fitStats$r.squared
    Stats[3]<-rmse(obs=ChemConc,sim=as.matrix(Predict))
    Stats[4]<-nrmse(obs=ChemConc,sim=as.matrix(Predict))
    Stats[5]<-fitStats$coefficients[2,1] #slope
  }
  
  if (subsample==0){
    return (list(Fit=Fit,ObservedAndPredicted=op3,Stats=Stats))
  }
  if (subsample>0){
    return (list(Fit=Fit,calibrationFit=calibrationFit,OaP1=op1,OaP2=op2,ObservedAndPredicted=op3,Stats=Stats))
    
  }
  
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Load fingerprint data and return them in a list%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
loadFingerPrints<-function(FingerPrintPath,filename,type){
    if(type=="original") {
          Row<-1
          data0<-read.table(file=paste(FingerPrintPath,filename[Row,1],sep=""),sep="\t",header=TRUE,skip=1,dec=".")
          data1<-read.table(file=paste(FingerPrintPath,filename[Row,2],sep=""),sep="\t",header=TRUE,skip=1,dec=",")
          data2<-read.table(file=paste(FingerPrintPath,filename[Row,3],sep=""),sep="\t",header=TRUE,skip=1,dec=",")
          data3<-read.table(file=paste(FingerPrintPath,filename[Row,4],sep=""),sep="\t",header=TRUE,skip=1,dec=",")
          }
    if(type=="1stDerivative")Row<-2
    if(type=="TurbidityCompensated")Row<-3  
    if(type=="1stDerivativeTurbidityCompensated")Row<-4

    if(Row>1){
              data0<-read.table(file=paste(FingerPrintPath,filename[Row,1],sep=""),sep="\t",header=TRUE,skip=1,dec=".")
              data1<-read.table(file=paste(FingerPrintPath,filename[Row,2],sep=""),sep="\t",header=TRUE,skip=1,dec=".")
              data2<-read.table(file=paste(FingerPrintPath,filename[Row,3],sep=""),sep="\t",header=TRUE,skip=1,dec=".")
              data3<-read.table(file=paste(FingerPrintPath,filename[Row,4],sep=""),sep="\t",header=TRUE,skip=1,dec=".")
              }
  #put them all together
  data<-rbind(data0,data1,data2,data3)
  
  #remove the individual files
  remove(data0,data1,data2,data3)
  
  #get rid of the files with bad probe status
  data<-data[data[,2]=="Ok",]
  
  #convert time to a dateTime Object
  realTime<-strptime(data[,1], "%Y.%m.%d  %H:%M:%S",tz="UTC")                     
  
  #remove  columns of NaNs in 742.5-750 nM bins
  fingerPrints<-data[,-(dim(data)[2]-3):-(dim(data)[2])]
  
  #remove first couple columns
  fingerPrints<-fingerPrints[,-1:-2]
  
  
  
  realTime<-realTime[complete.cases(fingerPrints[,2:dim(fingerPrints)[2]])]  
  fingerPrints<-fingerPrints[complete.cases(fingerPrints[,2:dim(fingerPrints)[2]]),]
  
  
  #align data to closest minute
  realTime<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)
  
  #left with fingerPrints and realTime
  #FingerPrints<-fingerPrints
  #RealTime<-realTime
  
  #get rid of the old container
  #remove(data)
  return(list(fingerPrints=fingerPrints,realTime=realTime))
}



