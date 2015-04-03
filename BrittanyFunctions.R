#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#
britFuncLoaded<-1

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

loadDataFile<-function(filepath,filename,excluderows=0,timezone="UTC"){
              data<-read.table(file=paste(filepath,filename,sep=""),sep=",",header=TRUE,skip=0)
              
              realTime<-strptime(data[,1],'%d/%m/%Y %H:%M:%S',tz=timezone)

              data$DateScan<-realTime                       #store it back in the matrix
            #exclude the rows specified in the environment setup--all observations from each of those rows is removed
            #it might be smarter to do it by realtime
              if(excluderows>0){
                keepThese<-!realTime%in%realTime[excluderows]
                }
              else{
                keepThese=seq(1:dim(data)[1])
              }
               
              realTime<-(data$DateScan[keepThese])                                       #pull out time
              ChemData<-data[keepThese,(dim(data)[2]-16):(dim(data)[2])]            #pull out the analyte data
              fingerPrints<-data[keepThese,-(dim(data)[2]-20):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
              fingerPrints<-fingerPrints[,-1:-2]                      #remove status and datetime
               realTime<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)
            
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
                fitQualityModelFull<-OB(ChemConcFull,calibPredict2,fitEval,fitFile,fitFileOut)
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
                system(paste(fitEval,fitFile,"jpg",sep=" "),wait=TRUE,show.output.on.console=FALSE)
                
                a<-readLines(fitFileOut)
                #unlink(fitFileOut)
                
                #in line a[7] there are data for very good fits, after tha : and before the %
                fitQuality<-c(0,0,0,0)
                names(fitQuality)<-c("veryGood","good","acceptable","bad")
                
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
    
  #pull apart again, and keep a copy for validating a subset model
  ChemConc<-as.matrix(fingerPrint[,dim(fingerPrint)[2]])                       #strip off the chemConc again
 # ChemConcFull<-ChemConc
 # RealTimeFull<-as.matrix(fingerPrint[,1])
 # fingerPrintFull<-fingerPrint[,-dim(fingerPrint)[2]]             #complete cases fingerprints (remove chemConc)
 # fingerPrintFull<-FingerPrintFull[,-1]
  
  if(subsample==0){  #if we are not going to take a subset of the data, follow this protocol
    Stats<-matrix(nrow=1,ncol=9)
    colnames(Stats)<-c("n","r2","rmse","nrmse","slope","vg","g","a","b")
    RealTime<-as.matrix(fingerPrint[,1])                        #complete cases realtime
    fingerPrint<-fingerPrint[,-dim(fingerPrint)[2]]             #complete cases fingerprints (remove chemConc)
    fingerPrint<-fingerPrint[,-1]                               #complete cases fingerprints (remove RealTime)
    
  }
  
  if(subsample > 0){#enter density dependent subsampling with Subsample Value guiding the process
                  Stats<-matrix(nrow=1,ncol=15)
                  colnames(Stats)<-c("n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope")
                  Histogram<-hist(ChemConc,breaks="FD",plot=FALSE)
                  sampleRatio<-(subsample/length(ChemConc))
                  #sampleRatio<-subsample
                  count<-0
                  # fingerPrint1<-fingerPrint[order(fingerPrint[,dim(fingerPrint)[2]]),] #SORTED IN ASCENDING ORDER
                  for (q in seq(1,(length(Histogram$breaks)-1),1)){
                    #for each histogram bin, make a list of fingerPrints
                    inTheBin<-fingerPrint[(ChemConc>=Histogram$breaks[q])&(ChemConc<=Histogram$breaks[q+1]),]
                    #print(inTheBin)
                    #estimate the number of samples we need from that bin (rounding up here)
                    sampleNum<-sum(ceiling(Histogram$counts[q]*sampleRatio))
                    #print(sampleNum)
                    #collect that number of samples from the bin
                    someNumber<-inTheBin[sample(1:dim(inTheBin)[1],sampleNum,replace=FALSE),]
                   # print(someNumber)
                    if(count==0){
                      calibration<-someNumber
                    }
                    if(count>0){
                      calibration<-rbind(calibration,someNumber)
                    }
                    count<-count+1
                  #  readline()
                  }
                  remove(someNumber)
          #restructure the calibration dataset to use in PLSR     
             #pull of the the chem concentration
                  calibrationChemConc<-(calibration[,dim(calibration)[2]])
             #and the realTime    
                  calibrationRealTime<-as.matrix(calibration[,1]) 
             #and drop their columns      
                  calibrationFingerPrint<-calibration[,-1]
                  calibrationFingerPrint<-calibrationFingerPrint[,-dim(calibrationFingerPrint)[2]]
                  
                  
          #restructure the full dataset to be used in PLSR and validation stages
              #pull of the the realTime and ChemConce data
                  ChemConc<-as.matrix(fingerPrint[,dim(fingerPrint)[2]])
                  RealTime<-as.matrix(fingerPrint[,1])  
              #and drop their columns     
                  fingerPrint<-fingerPrint[,-dim(fingerPrint)[2]]  
                  fingerPrint<-fingerPrint[,-1]                    
                  
#RUN PLSR MODEL FOR calibration SUBSAMPLE
                  calibrationFit<-plsr(calibrationChemConc~data.matrix(calibrationFingerPrint),
                                       ncomp=numParameters,validation="CV")
          
          #predict concentrations for the calibration (uninteresting)
                  calibPredict1<-predict(calibrationFit,data.matrix(calibrationFingerPrint),
                                         ncomp=numParameters,type=c("response"))
          #calculate the fit quality calling on FitEval
              #    fitQualityModel<-OB(calibrationChemConc,calibPredict1,fitEval,fitFile,fitFileOut)
                        
          #combine the observed and predicted values
                      op1<-cbind(calibrationChemConc,calibPredict1,calibrationRealTime)
                      colnames(op1)<-c("observed","predicted","realTime")
          #calculate some other fit statistics 
                      fitStats<-summary(lm(calibrationChemConc~calibPredict1))
                      confInt<-predict(lm(calibrationChemConc~calibPredict1),interval='prediction',level=0.95)
          #prepare output       
                  Stats[1]<-fitStats$df[2]
                  Stats[2]<-fitStats$r.squared
                  Stats[3]<-rmse(obs=as.matrix(calibrationChemConc),sim=as.matrix(calibPredict1))
                  Stats[4]<-nrmse(obs=as.matrix(calibrationChemConc),sim=as.matrix(calibPredict1))
                  Stats[5]<-fitStats$coefficients[2,1] #slope
                  
    #predict concentrations using the model PLSR Fit with the full dataset of fingerprints
                        #this is what I call validating the model
                  calibPredict2<-predict(calibrationFit,data.matrix(fingerPrint),ncomp=numParameters,type=c("response"))

          #calculate the fit quality calling on FitEval
             #     fitQualityModelFull<-OB(ChemConc,calibPredict2,fitEval,fitFile,fitFileOut)
          #combine the observed and predicted values
                      op2<-cbind(ChemConc,calibPredict2,RealTime)
          #find the values that were used in the ccalibration model by date (column 3)
          #this lets us keep track later on of which subset of values were used in making the model
                      indx<-op2[,3]%in%op1[,3]
          #add the index variable to the observed and predicted dataset
                      op2<-cbind(op2,indx)
                      colnames(op2)<-c("observed","predicted","realTime","calib")
          #calculate some summary statistics    using all the available ChemConc and all predicted values  
                      fitStats<-summary(lm(ChemConc~calibPredict2))
                      confInt<-predict(lm(ChemConc~calibPredict2),interval='prediction',level=0.95)
          #prepare output
                  Stats[6]<-fitStats$df[2]
                  Stats[7]<-fitStats$r.squared
                  Stats[8]<-rmse(obs=as.matrix(ChemConc),sim=as.matrix(calibPredict2))
                  Stats[9]<-nrmse(obs=as.matrix(ChemConc),sim=as.matrix(calibPredict2))
                  Stats[10]<-fitStats$coefficients[2,1] #slope
  }#end if subsample>0
  
  #run PLSR Model for all available data
  Fit<-plsr(ChemConc~data.matrix(fingerPrint),ncomp=numParameters,validation="CV")  #PLSR model to predict chemConc with cross validation
  #predicte chem concentrations
  Predict<-predict(Fit,data.matrix(fingerPrint),ncomp=numParameters,type=c("response"))
  #fitQualityFull<-OB(ChemConc,Predict,fitEval,fitFile,fitFileOut)
  #collect data for output
  op3<-cbind(ChemConc,Predict,RealTime)

  colnames(op3)<-c("observed","predicted","realTime")
 
  
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
    return (list(FitAll=Fit,Fit=calibrationFit,OaP1=op1,OaP2=op3,ObservedAndPredicted=op2,Stats=Stats))
    
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




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

selectSpecData<-function(fileType,dataType){
      if(fileType=="prunedO"){
        data<-original
        myData<-prunedO
        
      }
      if(fileType=="original"){
        data<-original
        myData<-originalmyData
      }
      
      if(fileType=="1stDer"){
        data<-Derivative
        myData<-DermyData
      }

      if(fileType=="pruned1D"){
        data<-Derivative
        myData<-pruned1D
      }
      
      
      if(fileType=="turbComp"){
        data<-TurbidityCompensated
        myData<-TCmyData
      }
      if(fileType=="prunedTC"){
        data<-TurbidityCompensated
        myData<-prunedTC
      }
      
      if(fileType=="1stDerTurbComp"){
        data<-FstDerivativeTurbidityCompensated
        myData<-TC1DmyData
      }
      if(fileType=="prunedTC1D"){
        data<-FstDerivativeTurbidityCompensated
        myData<-prunedTC1D
      }
      
     
  if(dataType=="fingerPrints"){
    output<-data #subset(data,dataType,startDates,stopDates)
  }
  
  if(dataType=="calibration"){
    output<-myData #subset(myData,dataType,startDates,stopDates,chem)
  }
return(output)  
}

# 
# subset<-function(data,type,startDates,stopDates,chem=NULL){
#   
#   logicalIndexofInclusion<-matrix(nrow=length(data$realTime),ncol=length(startDates))
#   
#   for(i in 1:length(startDates)){
#     logicalIndexofInclusion[,i]<-(data$realTime>startDates[i]&data$realTime<stopDates[i])
#   }
#   
#   keep<-as.logical(rowSums(logicalIndexofInclusion,na.rm=TRUE))
#   
#   
#   if (type=="calibration"){
#     realTime<-data$realTime[keep]
#     fingerprint<-data$fingerPrint[keep,]
#     ChemData<-data$ChemData[keep,]
#     
#     fp<-cbind(realTime,fingerprint,as.matrix(ChemData[,chem]))
#     goodData<-fp[complete.cases(fp[,2:dim(fp)[2]]),] #removes all the rows for which there is a NA--keeping time in there
#     realTime<-goodData[,1]                      #pull components back out
#     ChemData<-goodData[,dim(goodData)[2]]
#     fingerprints<-goodData[,-1] 
#     fingerprints<-fingerprints[,-dim(fingerprints)[2]]
#     
#     
#     
#     return(list(realTime=realTime,fingerPrints=fingerprints,ChemData=ChemData))
#     #returnStatement
#   }
#   
#   if (type=="fingerPrints"){
#     realTime<-data$realTime[keep]
#     fingerprint<-data$fingerPrint[keep,]
#     return(list(realTime=realTime,fingerPrints=fingerprint))
#     
#   }
#   
#   if(type=="flow"){
#     realTime<-data$realTime[keep]
#     flow<-data$flow[keep]
#     return(list(realTime=realTime,flow=flow))
#   }
#   
# }



modelExecution<-function(numComp,subsampleRate,calibData,dataToModel,fitEval,fitFile,fitFileOut){
  modelQuality<-matrix(ncol=26)
  colnames(modelQuality)<-c("M_n","M_NSE",
                            "M_VG","M_G","M_A","M_B",
                            "V_VG","V_G","V_A","V_B",
                            "C_n","C_r2","C_RMSE","C_NRMSE","C_slope",
                            "V_n","V_r2","V_RMSE","V_NRMSE","V_slope",
                            "F_n","F_r2","F_RMSE","F_NRMSE","F_slope",
                            "nsubset")
  
  #for any model--first subset the data, and pass subsetted data to the execution function
  if(length(calibData$ChemData)>15){
    ModelConcentration<-PLSRFitAndTestNoNSE(calibData$fingerPrints,
                                            calibData$ChemData,
                                            calibData$realTime,
                                            numComp,
                                            fitEval,fitFile,fitFileOut,subsampleRate)                  
    #USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
    PredictedConcentration<-predict(ModelConcentration$Fit,
                                    as.matrix(dataToModel$fingerPrints),
                                    ncomp=numComp,
                                    type=c("response"))  
    PredictedConcentration<-cbind(as.numeric(dataToModel$realTime),PredictedConcentration)
    
    #calculate stats for goodness of fit etc
 #$$$$   
    
      modelQuality[1]<-ModelConcentration$Stats[1]
    if(subsampleRate>0){
      modelQuality[2]<-ModelConcentration$Stats[3]
      modelQuality[3:6]<-OB(ModelConcentration$OaP1[,1],ModelConcentration$OaP1[,2],fitEval,fitFile,fitFileOut)
      modelQuality[7:10]<-OB(ModelConcentration$ObservedAndPredicted[,1],ModelConcentration$ObservedAndPredicted[,2],fitEval,fitFile,fitFileOut)
    }
    if(subsampleRate==0){
      modelQuality[2]<-ModelConcentration$Stats[3]
      modelQuality[3:6]<-OB(ModelConcentration$ObservedAndPredicted[,1],ModelConcentration$ObservedAndPredicted[,2],fitEval,fitFile,fitFileOut)
    }
    modelQuality[11:25]<-ModelConcentration$Stats
    modelQuality[26]<-subsampleRate
    #return summary statistics, and chemographic prediction vector (in case someone wants to keep it)
    #return(PredictedConcentration)
    #this can be repeated 1x or
    return(list(modelQuality=modelQuality,PredictedConcentrations=PredictedConcentration,ObservedandPredicted=ModelConcentration$ObservedAndPredicted))
    
  }
  else return(list(modelQuality=modelQuality,PredictedConcentrations=c(NaN,NaN),ObservedandPredicted=c(NaN,NaN,NaN)))
  
}


#foobar<-modelExecution(chemical,numComp,subsetRate,calibration,specDataToModel,fitEval,fitFile,fitFileOut)





#in this one chem is a character string that has to match the name in the data column

subsetAll<-function(calibData=0,modelData=0,flow=0,chem,dateWindows,keepSpecData=0,rangeOrwindows="windows"){
  #for each type of data (calibData, spec data to model, and flow data), check to see if it is present,
  # if it is, subset it to the date windows specified
  #bad data rows have already been excluded
  calib<-0
  model<-0
  floutput<-0
  
  #if rangeOrWindows is defined as "range" only the outer limits will be used for calibrationData
  #obviously, this should force the program to "keepSpecData=1" as well.
  
  #if keepSpecData is set to 1, spec data will only be hacked down to the min and max of the date range
  #windows of data omitted from the calib data will be kept in the spec data
  
  if(rangeOrwindows=="range"){
    keepSpecData=1;
  }
  
  
  #if there are calibdata coming in, subset them
  if(length(as.data.frame(calibData))>1){  
                        if(rangeOrwindows=="windows"){
                          keep<-(calibData$realTime>min(dateWindows)&calibData$realTime<max(dateWindows))
                    
                        }
                  else{
                        logicalIndexofInclusion<-matrix(nrow=length(calibData$realTime),ncol=dim(dateWindows)[1])
                        for(i in 1:dim(dateWindows)[1]){
                          logicalIndexofInclusion[,i]<-(calibData$realTime>dateWindows[i,1]&calibData$realTime<dateWindows[i,2])
                        }
                        
                        keep<-as.logical(rowSums(logicalIndexofInclusion,na.rm=TRUE))
                        }
                  
                      calibData<-cbind(calibData$realTime[keep],calibData$fingerPrint[keep,],subset(calibData$ChemData[keep,],select=chem))
                      colnames(calibData)[1]<-"realTime"    
                      
                      calibData<-calibData[complete.cases(calibData[,2:dim(calibData)[2]]),] #removes all the rows for which there is a NA--keeping time in there
                      crealTime<-calibData[,1]                      #pull components back out
                      cChemData<-calibData[,dim(calibData)[2]]
                      cfingerprints<-calibData[,-1] 
                      cfingerprints<-cfingerprints[,-dim(cfingerprints)[2]]
                      crealTime<-align.time(crealTime,60)
                      calibData$realTime<-crealTime
                      calib<-list(realTime=crealTime,fingerPrints=cfingerprints,ChemData=cChemData)
                }
  
    if(length(as.data.frame(modelData))>1){
                  if(keepSpecData==0){
                        logicalIndexofInclusion<-matrix(nrow=length(modelData$realTime),ncol=dim(dateWindows)[1])
                            
                          for(i in 1:dim(dateWindows)[1]){
                              logicalIndexofInclusion[,i]<-(modelData$realTime>dateWindows[i,1]&modelData$realTime<dateWindows[i,2])
                            }
                        
                        keep<-as.logical(rowSums(logicalIndexofInclusion,na.rm=TRUE))
              
                        }
                
                  else{  
                        keep<-(modelData$realTime>min(dateWindows) & modelData$realTime<max(dateWindows))
                        }
                
                        modelData<-cbind(modelData$realTime[keep],modelData$fingerPrint[keep,])
                        colnames(modelData)[1]<-"realTime"    
                        
                        modelData<-modelData[complete.cases(modelData[,2:dim(modelData)[2]]),] #removes all the rows for which there is a NA--keeping time in there
                        mrealTime<-modelData[,1]                      #pull components back out
                        mfingerprints<-modelData[,-1] 
                        mrealTime<-align.time(mrealTime,60)
                        model<-list(realTime=mrealTime,fingerPrints=mfingerprints)  
                }
          
  #check to see if there are flow data. if so, subset them
  if(length(as.data.frame(flow))>1){  
                  keep<-(flow$realTime>dateWindows[1,1] & flow$realTime<dateWindows[dim(dateWindows)[1],2])
                  flowTime<-flow$realTime[keep]
                  flowFlow<-flow$flow[keep]
                  floutput<-list(realTime=flowTime,flow=flowFlow)
                }
  return(list(calibData=calib,specData=model,flow=floutput))
  


}


runModel<-function(chem,numComp,dates,calibData,modelData,type,flow,subsetRatio=0,iterations=10,keepSpecData=1,fitEval,fitFile,fitFileOut){
  #subset data to the windows specified above
  lotOfData<-subsetAll(calibData=calibData,modelData=modelData,flow=flow,chem=chem,dateWindows=dates,keepSpecData=keepSpecData)
  #run model
  for (r in 1:iterations){
    
    putout<-modelExecution(numComp=numComp,subsample=100,
                           calibData=lotOfData$calibData,
                           dataToModel=lotOfData$specData,
                           fitEval,fitFile,fitFileOut)
    

    
    
    if(r==1){
      modelQuality<-as.data.frame(matrix(nrow=iterations,ncol=dim(putout$modelQuality)[2]+5))
      modelQuality[1,]<-c(putout$modelQuality,chem,numComp,r,subsetRatio,type)
      colnames(modelQuality)<-c(colnames(putout$modelQuality),"chem","numComp","iter","subset","type")
      
      predictions<-as.data.frame(matrix(nrow=dim(putout$PredictedConcentrations)[1],ncol=iterations+1))
      predictions[,1:2]<-putout$PredictedConcentrations
      
      Observed<-as.data.frame(lotOfData$calibData$realTime)
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





