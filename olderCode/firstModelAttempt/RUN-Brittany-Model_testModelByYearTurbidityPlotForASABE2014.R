library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package
library(xts)
source("C:/Users/FBlab/Documents/GitHub/Brittany/BrittanyFunctions.R", echo=TRUE)

#******Specify file paths and names
#FPpath<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\FichiersFP\\" #Specify folder where data is located
#path<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\"        #upstairs box
#fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/"  #fiteval_out.txt location on upstairs box

fitPath<-"C:/Users/FBlab/Desktop/FITEVAL2_win/FITEVAL2_win/"    #location on downstairs box
path<-"C:/Users/FBlab/Desktop/work_here/Data/"   #downstairs box
FPpath<-"C:/Users/FBlab/Desktop/work_here/Data/FichiersFP/" #downstairs box
flow<-read.table(file="C:/Users/FBlab/Desktop/work_here/Data/Brittany/flow/flowKervidyNaizin_2010-2012.txt",skip=1)
Date<-flow$V1        #the date
T<-flow$V2     #the time
D<-paste(Date,T,sep=" ")               #put back in the fixed values
D<-strptime(D, '%d/%m/%Y %H:%M:%S',tz="UTC")    #convert it to a _flawless_ time opject
flow$realTime<-D    
flow$flow<-flow$V3
flow<-flow[,-1:-3]
flow<-flow[complete.cases(flow$flow),]
flow<-flow[complete.cases(flow$realTime),]# removes all the rows for which there is a NA--keeping time in there
remove(D,T,Date)


fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
#filename<-c("S2010-2011TurbComp.fp" ,"S2011-2012TurbComp.fp","S2012-2013TurbComp.fp","S2013-2014TurbComp.fp")
filename<-c("S2010-20111stDer.fp", "S2011-20121stDer.fp","S2012-20131stDer.fp","S2013-20141stDer.fp")
#LOAD DATA FROM BRITTANY FOR PROJECTION
      data0<-read.table(file=paste(FPpath,filename[1],sep=""),sep="\t",header=TRUE,skip=1)
      data1<-read.table(file=paste(FPpath,filename[2],sep=""),sep="\t",header=TRUE,skip=1)
      data<-rbind(data0,data1)
      remove(data0,data1)

#parse the date
      realTime<-strptime(data$Date.Time, "%Y.%m.%d  %H:%M:%S",tz="UTC")                       #store it back in the matrix
      fingerPrints<-data[,-(dim(data)[2]-3):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
      fingerPrints<-fingerPrints[,-1:-2]  
      remove(data)

#align data to closest minute
foo<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)

#meh<-approx(flow[71462:71562,1],flow[71462:71562,2],method="linear",n=200)

#filename<-c("S2010-20111stDerTurbComp.fp" ,"S2011-20121stDerTurbComp.fp","S2012-20131stDerTurbComp.fp","S2013-20141stDerTurbComp.fp")
#LOAD DATA FOR MODEL AND MAKE PLSR MODEL BASED ON WHAT WE KNOW ABOUT WHATEVER
              ModelFilename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")
              Chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
                      "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity");
            Turbidity<-17
            startYear<-c(2010)
            stopYear<-c(2011)
#IDENTIFTY TIME WINDOW
          startDate<-as.POSIXct(paste(startYear[1],"-08-01 00:00:00",sep=""),tz="UTC")
          stopDate<-as.POSIXct(paste(stopYear[1],"-10-01 00:00:00",sep=""),tz="UTC")
#LOAD DATA
            myData<-loadDataFile(path,ModelFilename[2]) 

#SUBSET THE FINGERPRINT DATA THAT WE ARE HOPING TO MODEL
            useUs<-(realTime[realTime>startDate&realTime<stopDate])
            useUsfp<-fingerPrints[(realTime>startDate&realTime<stopDate) ,]
#but it isn't exactly within the range we specified, so, lets get the correct range here
      startDate<-min(useUs)
      stopDate<-max(useUs)
                                      #subset the CALIBRATION data by the time window
                                      #useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
                                      #useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])                                      #useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])
                                      
                                          


j<-1 #the number of times to run through each model style
n=64 #NUMBER OF CALIBRATION POINTS
numComp<-11

#
for (awesome in 1:5){

RNSE<-matrix(nrow=j,ncol=4)
for(i in 1:j){
#OPTION 1 FOR SUBSETTING/CALIBRATION
              #REFRESH DATA FOR SUBSETTING
                        useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
                        useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
                        useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])
                        #CALCULATE NUMBER OF COMPONENTS TO USE
                        #numComp<-numberOfComponentsToUse(useUsFP,useUsChem$Turb)
                       # numComp<-17 #MAYBE JUST FORCE IT TO USE A REASONABLE NUMBER
              #REMOVE ALL NAN POINTS
                        fp<-cbind(useUsRealTIme,useUsFP,as.matrix(useUsChem[,17]))
                        foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] # removes all the rows for which there is a NA--keeping time in there
                        useUsRealTIme<-foo[,1]
                        useUsChem<-foo[,dim(foo)[2]]
                        useUsFP<-foo[,-1]
                        useUsFP<-useUsFP[,-dim(useUsFP)[2]]
             
              #SUBSET CALIBRATION DATA RANDOMLY
                          calibrationPoints<-sample(1:length(useUsChem),n,replace=TRUE)
                          #numComp<-11#numberOfComponentsToUse(useUsFP[calibrationPoints,],useUsChem[calibrationPoints])
                          hist(useUsChem[calibrationPoints])
                          #GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
                          TurbModelRandom<-PLSRFitAndTest(useUsFP[calibrationPoints,],useUsChem[calibrationPoints],useUsRealTIme[calibrationPoints],numComp,fitEval,fitFile,fitFileOut,0)
                          #USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
                          PredictRandom<-predict(TurbModelRandom$Fit,data.matrix(useUsfp),ncomp=numComp,type=c("response"))
              
              RNSE[i,1:4]<-TurbModelRandom$fitQuality
}
#OPTION 2 FOR CALIBRATION/SUBSETTING
LPNSE<-matrix(nrow=j,ncol=4)
for(i in 1:j){
            #REFRESH DATA FOR SUBSETTING
                          useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
                          useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
                          useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])
            #REMOVE ALL NAN POINTS
                          fp<-cbind(useUsRealTIme,useUsFP,as.matrix(useUsChem[,17]))
                          foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] # removes all the rows for which there is a NA--keeping time in there
                          useUsRealTIme<-foo[,1]
                          useUsChem<-foo[,dim(foo)[2]]
                          useUsFP<-foo[,-1]
                          useUsFP<-useUsFP[,-dim(useUsFP)[2]]

              #IDENTIFY CUTOFF
                          cutOff<-40
                       #   n=64
              #KEEP THOSE GREATER THAN value
                          useUsPoints<-which(useUsChem<cutOff)
                          useUsRealTIme<-useUsRealTIme[useUsPoints]
                          useUsChem<-useUsChem[useUsPoints]
                          useUsFP<-useUsFP[useUsPoints,]
              #GRAB A RANDOM SUBSET OF N POINTS
                          calibrationPoints<-sample(1:length(useUsChem),n,replace=TRUE)
                         # numComp<-11#numberOfComponentsToUse(useUsFP[calibrationPoints,],useUsChem[calibrationPoints])
                          hist(useUsChem[calibrationPoints])
              #GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
                          TurbModelLowerPortion<-PLSRFitAndTest(useUsFP[calibrationPoints,],useUsChem[calibrationPoints],useUsRealTIme[calibrationPoints],numComp,fitEval,fitFile,fitFileOut,0)
              #USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
                          PredictLowerPortion<-predict(TurbModelLowerPortion$Fit,data.matrix(useUsfp),ncomp=numComp,type=c("response"))
            LPNSE[i,1:4]<-TurbModelLowerPortion$fitQuality
}
#KEEP THOSE GREATER THAN 50
UPNSE<-matrix(nrow=j,ncol=4)
for(i in 1:j){
              #REFRESH DATA FOR SUBSETTING
                          useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
                          useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
                          useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])
              #REMOVE ALL NAN POINTS
                          fp<-cbind(useUsRealTIme,useUsFP,as.matrix(useUsChem[,17]))
                          foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] # removes all the rows for which there is a NA--keeping time in there
                          useUsRealTIme<-foo[,1]
                          useUsChem<-foo[,dim(foo)[2]]
                          useUsFP<-foo[,-1]
                          useUsFP<-useUsFP[,-dim(useUsFP)[2]]
                          
              #IDENTIFY CUTOFF
                          cutOff<-40
                          #n=64
              #KEEP THOSE GREATER THAN value
                          useUsPoints<-which(useUsChem>=cutOff)
                          useUsRealTIme<-useUsRealTIme[useUsPoints]
                          useUsChem<-useUsChem[useUsPoints]
                          useUsFP<-useUsFP[useUsPoints,]
              #GRAB A RANDOM SUBSET OF N POINTS
                          calibrationPoints<-sample(1:length(useUsChem),n,replace=TRUE)
                          #numComp<-11#numberOfComponentsToUse(useUsFP[calibrationPoints,],useUsChem[calibrationPoints])
                          hist(useUsChem[calibrationPoints])
              #GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
                          TurbModelUpperPortion<-PLSRFitAndTest(useUsFP[calibrationPoints,],useUsChem[calibrationPoints],useUsRealTIme[calibrationPoints],numComp,fitEval,fitFile,fitFileOut,0)
              #USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
                          PredictUpperPortion<-predict(TurbModelUpperPortion$Fit,data.matrix(useUsfp),ncomp=numComp,type=c("response"))
              UPNSE[i,1:4]<-TurbModelUpperPortion$fitQuality
}
#Try the Density Dependent Sample Routine
DDNSE<-matrix(nrow=j,ncol=4)
for(i in 1:j){
              #REFRESH DATA FOR SUBSETTING
                          useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
                          useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
                          useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])
              #REMOVE ALL NAN POINTS
                          fp<-cbind(useUsRealTIme,useUsFP,as.matrix(useUsChem[,17]))
                          foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] # removes all the rows for which there is a NA--keeping time in there
                          useUsRealTIme<-foo[,1]
                          useUsChem<-foo[,dim(foo)[2]]
                          useUsFP<-foo[,-1]
                          useUsFP<-useUsFP[,-dim(useUsFP)[2]]
              
              #IDENTIFY CUTOFF
                         # n<-64
              #KEEP THOSE GREATER THAN value
              subsetRatio<-round(n/length(useUsChem),2)
              subset<-densityDependentSubset(useUsChem,useUsRealTIme,useUsFP,subsetRatio,TRUE)
              useUsRealTIme<-subset$realTime
              useUsChem<-subset$Chem
              useUsFP<-subset$fingerPrint
              #GRAB A RANDOM SUBSET OF N POINTS
              calibrationPoints<-sample(1:length(useUsChem),n,replace=TRUE)
             # numComp<-11#numberOfComponentsToUse(useUsFP[calibrationPoints,],useUsChem[calibrationPoints])
              hist(useUsChem[calibrationPoints],breaks="FD")
              #GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
              TurbModelDD<-PLSRFitAndTest(useUsFP[calibrationPoints,],useUsChem[calibrationPoints],useUsRealTIme[calibrationPoints],numComp,fitEval,fitFile,fitFileOut,0)
              #USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
              PredictDD<-predict(TurbModelDD$Fit,data.matrix(useUsfp),ncomp=numComp,type=c("response"))
              DDNSE[i,1:4]<-TurbModelDD$fitQuality
}



#SUBSET THE FLOW DATA
uuFlowTime<-flow$realTime[flow$realTime>startDate & flow$realTime<stopDate]
uuFlowFlow<-flow$flow[flow$realTime>startDate & flow$realTime<stopDate]
    begin<-startDate
    end<-stopDate
totmin<-as.numeric(end-begin)*60*24

flowApprox<-approx(uuFlowTime,uuFlowFlow,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset

AP_R<-approx(useUs,PredictRandom,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset
AP_R$y[(AP_R$y<0)]<-0

AP_UP<-approx(useUs,PredictUpperPortion,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset
AP_UP$y[(AP_UP$y<0)]<-0

AP_LP<-approx(useUs,PredictLowerPortion,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset
AP_LP$y[(AP_LP$y<0)]<-0

AP_DD<-approx(useUs,PredictDD,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset
AP_DD$y[(AP_DD$y<0)]<-0


#load=flow*turbidity
date<-as.POSIXct(flowApprox$x,origin="1970-1-1 0:0:0",tz="UTC")
load_R<-AP_R$y*flowApprox$y
load_UP<-AP_UP$y*flowApprox$y
load_LP<-AP_LP$y*flowApprox$y
load_DD<-AP_DD$y*flowApprox$y

plot(date,cumsum(load_UP),col="black",type="l")
points(date,cumsum(load_R),col="red",type="l")
points(date,cumsum(load_LP),col="orange",type="l")
points(date,cumsum(load_DD),col="green",type="l")


}