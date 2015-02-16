#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
if(fileType==1){
  fingerPrints<-original$fingerPrint
  realTime<-original$realTime
  myData<-originalmyData
}

if(fileType==2){
  fingerPrints<-Derivative$fingerPrint
  realTime<-Derivative$realTime
  myData<-DermyData
}

if(fileType==3){
  fingerPrints<-TurbidityCompensated$fingerPrint
  realTime<-TurbidityCompensated$realTime
  myData<-TCmyData
}

if(fileType==4){
  fingerPrints<-FstDerivativeTurbidityCompensated$fingerPrint
  realTime<-FstDerivativeTurbidityCompensated$realTime
  myData<-TC1DmyData
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ANNUAL MODEL  SUBSET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#SUBSETDATA
          startYear<-c(2010)
          stopYear<-c(2011)
       
        #IDENTIFTY TIME WINDOW FOR THE "ANNUAL MODEL"
            startDate<-as.POSIXct(paste(startYear[1],"-10-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
        #startDate<-as.POSIXct(paste(startYear[1],"-11-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
            stopDate<-as.POSIXct(paste(stopYear[1],"-7-1 00:00:00",sep=""),tz="UTC")    #june 15 2011

        #select the time range for the event scale load estimation
            eventStart<-as.POSIXct("2010-12-01 07:13:38",tz="UTC")
            eventStop<-as.POSIXct("2010-12-15 04:33:31",tz="UTC")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
        #SUBSET THE FINGERPRINT DATA THAT WE ARE HOPING TO MODEL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
        #rough first pass cut
            subsetRealTime<-(realTime[realTime>startDate&realTime<stopDate])
            subsetFingerPrints<-fingerPrints[(realTime>startDate&realTime<stopDate) ,]
        
        #identify the window of GOOOD data (excuding the fouling)
            goodStart<-as.POSIXct(paste("2010-11-18 00:00:00",sep=""),tz="UTC")
            goodStop<-as.POSIXct(paste("2011-05-1 00:00:00",sep=""),tz="UTC")
            
            good2Start<-as.POSIXct(paste("2011-05-1 00:00:00",sep=""),tz="UTC")
            good2Stop<-as.POSIXct(paste("2011-06-15 00:00:00",sep=""),tz="UTC")

            #use the realtime to get a logical vector of data that fit into the window of interest
            bar<-(realTime>goodStart&realTime<goodStop | 
                    realTime>good2Start&realTime<good2Stop )

          #here the 2 indicates that it is the second round of time-subsetting
            subsetRealTime2<-realTime[bar]
            subsetFingerPrints2<-fingerPrints[bar,]

  #but it isn't exactly within the range we specified, so, lets get the correct range here
            startDate<-min(useUs[!is.na(useUs)])
            stopDate<-max(useUs[!is.na(useUs)])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
          #subset the CALIBRATION data by the time window
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
          calibrationRealTime<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
          calibrationFingerPrints<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
          calibrationAnalytes<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])
          
          #use the realtime to get a logical vector of data that fit into the window of interest
          bar<-(myData$realTime>goodStart&myData$realTime<goodStop | 
                  myData$realTime>good2Start&myData$realTime<good2Stop )
          
          calibrationRealTime2<-myData$realTime[bar]
          calibrationFingerPrints2<-myData$fingerPrints[bar,]
          calibrationAnalytes2<-myData$ChemData[bar,]



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#subset flow for annual model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
        AnnualFlowTime<-flow$realTime[flow$realTime>startDate & flow$realTime<stopDate]
        AnnualFlowFlow<-flow$flow[flow$realTime>startDate & flow$realTime<stopDate]
        
        begin<-startDate
        end<-stopDate

        totmin<-as.numeric(end-begin)*60*24

        #calculate the flow values for all 1 minute intervals within the year of interest
        AnnualFlowApprox<-approx(AnnualFlowTime,AnnualFlowFlow,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset

        #subsetFlow for Load calculation
        eventFlowTime<-AnnualFlowApprox$x[AnnualFlowApprox$x>eventStart & AnnualFlowApprox$x<eventStop]
        eventFlowFlow<-AnnualFlowApprox$y[AnnualFlowApprox$x>eventStart & AnnualFlowApprox$x<eventStop]


  #REMOVE ALL NAN POINTS from 2 different datasets
  #first a group that are not filtered for quality
        fp<-cbind(calibrationRealTime,calibrationFingerPrints,as.matrix(calibrationAnalytes[,chemN[chemical]]))
        foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] #removes all the rows for which there is a NA--keeping time in there
        calibrationRealTime<-foo[,1]                      #pull components back out
        calibrationAnalytes<-foo[,dim(foo)[2]]
        calibrationFingerPrints<-foo[,-1] 
        calibrationFingerPrints<-calibrationFingerPrints[,-dim(calibrationFingerPrints)[2]]

    #second a group that have been subsetted more agressively
        fp<-cbind(calibrationRealTime2,calibrationFingerPrints2,as.matrix(calibrationAnalytes2[,chemN[chemical]]))
        foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] #removes all the rows for which there is a NA--keeping time in there
        calibrationRealTime2<-foo[,1]                      #pull components back out
        calibrationAnalytes2<-foo[,dim(foo)[2]]
        calibrationFingerPrints2<-foo[,-1] 
        calibrationFingerPrints2<-calibrationFingerPrints2[,-dim(calibrationFingerPrints2)[2]]
        

