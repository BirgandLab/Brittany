
if(fool==1){
  fingerPrints<-FingerPrints
  realTime<-RealTime<-realTime
  myData<-OrigmyData
}

if(fool==2){
  fingerPrints<-DerFingerPrints
  realTime<-DerRealTime<-realTime
  myData<-DermyData
}

if(fool==3){
  fingerPrints<-TCFingerPrints
  realTime<-TCRealTime<-realTime
  myData<-TCmyData
}

if(fool==4){
  fingerPrints<-TC1DFingerPrints
  realTime<-TC1DRealTime
  myData<-TC1DmyData
}



#SUBSETDATA
startYear<-c(2010)
stopYear<-c(2011)
#IDENTIFTY TIME WINDOW FOR THE "ANNUAL MODEL"
startDate<-as.POSIXct(paste(startYear[1],"-10-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
#startDate<-as.POSIXct(paste(startYear[1],"-11-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
stopDate<-as.POSIXct(paste(stopYear[1],"-7-1 00:00:00",sep=""),tz="UTC")    #june 15 2011


#SUBSET THE FINGERPRINT DATA THAT WE ARE HOPING TO MODEL
useUs<-(realTime[realTime>startDate&realTime<stopDate])
useUsfp<-fingerPrints[(realTime>startDate&realTime<stopDate) ,]



#identify the window of GOOOD data (excuding the fouling)
goodStart<-as.POSIXct(paste("2010-11-18 00:00:00",sep=""),tz="UTC")
goodStop<-as.POSIXct(paste("2011-05-1 00:00:00",sep=""),tz="UTC")

good2Start<-as.POSIXct(paste("2011-05-1 00:00:00",sep=""),tz="UTC")
good2Stop<-as.POSIXct(paste("2011-06-15 00:00:00",sep=""),tz="UTC")

bar<-(realTime>goodStart&realTime<goodStop | 
        realTime>good2Start&realTime<good2Stop )
useUs2<-realTime[bar]
useUsfp2<-fingerPrints[bar,]


#but it isn't exactly within the range we specified, so, lets get the correct range here
startDate<-min(useUs[!is.na(useUs)])
stopDate<-max(useUs[!is.na(useUs)])



#select the time range for the event scale load estimation
eventStart<-as.POSIXct("2010-12-01 07:13:38",tz="UTC")
eventStop<-as.POSIXct("2010-12-15 04:33:31",tz="UTC")

#subset the CALIBRATION data by the time window
useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])

bar<-(myData$realTime>goodStart&myData$realTime<goodStop | 
        myData$realTime>good2Start&myData$realTime<good2Stop )

useUsRealTIme2<-myData$realTime[bar]
useUsFP2<-myData$fingerPrints[bar,]
useUsChem2<-myData$ChemData[bar,]

#SUBSET THE FLOW DATA FOR ANNUAL MODEL
uuFlowTime<-flow$realTime[flow$realTime>startDate & flow$realTime<stopDate]
uuFlowFlow<-flow$flow[flow$realTime>startDate & flow$realTime<stopDate]
begin<-startDate
end<-stopDate
totmin<-as.numeric(end-begin)*60*24
#calculate the flow values for all 1 minute intervals within the year of interest
flowApprox<-approx(uuFlowTime,uuFlowFlow,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset

#subsetFlow for Load calculation
eventFlowTime<-flowApprox$x[flowApprox$x>eventStart & flowApprox$x<eventStop]
eventFlowFlow<-flowApprox$y[flowApprox$x>eventStart & flowApprox$x<eventStop]







#REMOVE ALL NAN POINTS
fp<-cbind(useUsRealTIme,useUsFP,as.matrix(useUsChem[,chemN[chemical]]))
foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] #removes all the rows for which there is a NA--keeping time in there
useUsRealTIme<-foo[,1]                      #pull components back out
useUsChem<-foo[,dim(foo)[2]]
useUsFP<-foo[,-1] 
useUsFP<-useUsFP[,-dim(useUsFP)[2]]

fp<-cbind(useUsRealTIme2,useUsFP2,as.matrix(useUsChem2[,chemN[chemical]]))
foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] 
useUsRealTIme2<-foo[,1]    
useUsChem2<-foo[,dim(foo)[2]]
useUsFP2<-foo[,-1]
useUsFP2<-useUsFP2[,-dim(useUsFP2)[2]]

#figure out how many components to use in plsr model
#numComp<-numberOfComponentsToUse(useUsFP,useUsChem$MES)
#numComp<-17 #MAYBE JUST FORCE IT TO USE A REASONABLE NUMBER
#numComp<-NUMCOMP[chemical]#specify the number in a vector (above)
