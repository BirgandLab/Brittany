library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package
library(xts)
library(hydroTSM)


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
#filename<-c("S2010-2011.fp" ,"S2011-2012.fp","S2012-2013.fp","S2013-2014.fp")
#filename<-c("S2010-2011TurbComp.fp" ,"S2011-2012TurbComp.fp","S2012-2013TurbComp.fp","S2013-2014TurbComp.fp")
#filename<-c("S2010-20111stDer.fp", "S2011-20121stDer.fp","S2012-20131stDer.fp","S2013-20141stDer.fp")
filename<-c("2010-20111stDerTurbComp.fp", "2011-20121stDerTurbComp.fp","2012-20131stDerTurbComp.fp","2013-20141stDerTurbComp.fp")
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
realTime<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)

#filename<-c("S2010-20111stDerTurbComp.fp" ,"S2011-20121stDerTurbComp.fp","S2012-20131stDerTurbComp.fp","S2013-20141stDerTurbComp.fp")
#LOAD DATA FOR MODEL AND MAKE PLSR MODEL BASED ON WHAT WE KNOW ABOUT WHATEVER
              ModelFilename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")
              Chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
                      "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity");
            #Turbidity<-17
            startYear<-c(2010)
            stopYear<-c(2011)
#IDENTIFTY TIME WINDOW
          startDate<-as.POSIXct(paste(startYear[1],"-11-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
          stopDate<-as.POSIXct(paste(stopYear[1],"-6-15 00:00:00",sep=""),tz="UTC")    #june 15 2011
#select the time range for the event scale load estimation
            eventStart<-as.POSIXct("2010-12-01 07:13:38",tz="UTC")
            eventStop<-as.POSIXct("2010-12-15 04:33:31",tz="UTC")


#SUBSET THE FINGERPRINT DATA THAT WE ARE HOPING TO MODEL
            useUs<-(realTime[realTime>startDate&realTime<stopDate])
            useUsfp<-fingerPrints[(realTime>startDate&realTime<stopDate) ,]
#but it isn't exactly within the range we specified, so, lets get the correct range here
      startDate<-min(useUs)
      stopDate<-max(useUs)
                                      #subset the CALIBRATION data by the time window
                                      #useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
                                      #useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])                                      #useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])

#SUBSET THE FLOW DATA
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



#parameters to use when iterativelly loading and testing data
#FOR MES

chemicals<-c("DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL","NO2","NTOT","NTOTFILT")#,"SILICA")
mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2)#,2)#number corresponding to element of ModelFilename to use
chemN<-c(8,7,12,5,10,11,6,17,1,2,14,15)#,16)#number in chemistry vector
NUMCOMP<-c(9,5,5,5,10,10,9,9,9,9,2,8)#,8)#component number to use
#changed ptot and ppo43 to 10 components
calibPoints<-c(13,17,24,36,40,46,52,60,73,91,104)#,122,183,245,365) #NUMBER OF CALIBRATION POINTS
#changed min ncomp to 13 to accomodate 10 point component models
#12 is 1/30 days
#17 is 1/21 days 1/3 weeks
#24 is 1/15 days 1/2 weeks
#36 is 1/1o days
#40 is 1/9 days
#46 is 1/8 days
#52 is 1/7 days
#60 is 1/6 days
#73 is 1/5 days
#91 is 1/4 days
#104 is 2/week
#122 is 1/3 days
#183 is 1/2 days
#245 is 1/1.5 days
#365 is 1/1 day




j<-100 #the number of times to run through each model style
RNSE<-matrix(nrow=j,ncol=length(calibPoints))
ENSE<-matrix(nrow=j,ncol=length(calibPoints))
load_R<-matrix(nrow=j,ncol=length(calibPoints))
load_E<-matrix(nrow=j,ncol=length(calibPoints))

for (chemical in 5:6){#1:length(chemicals)){
  myData<-loadDataFile(path,ModelFilename[mfN[chemical]])       
       
        #OPTION 1 FOR SUBSETTING/CALIBRATION
        #REFRESH DATA FOR SUBSETTING
        useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
        useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
        useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])
        #CALCULATE NUMBER OF COMPONENTS TO USE
        #numComp<-numberOfComponentsToUse(useUsFP,useUsChem$MES)
        # numComp<-17 #MAYBE JUST FORCE IT TO USE A REASONABLE NUMBER
        #REMOVE ALL NAN POINTS
        #change chem here?
        fp<-cbind(useUsRealTIme,useUsFP,as.matrix(useUsChem[,chemN[chemical]]))
        
        numComp<-NUMCOMP[chemical]      
  
        foo<-fp[complete.cases(fp[,2:dim(fp)[2]]),] # removes all the rows for which there is a NA--keeping time in there
        useUsRealTIme<-foo[,1]
        useUsChem<-foo[,dim(foo)[2]]
        useUsFP<-foo[,-1]
        useUsFP<-useUsFP[,-dim(useUsFP)[2]]
       
for (awesome in 1:length(calibPoints)){
          n<-calibPoints[awesome]
            for(i in 1:j){
               
                   #SUBSET CALIBRATION DATA RANDOMLY
                                      if(n<=length(useUsChem)){
                                      calibrationPoints<-sample(1:length(useUsChem),n,replace=FALSE)
                                      }
                                      if(n>length(useUsChem)){
                                        calibrationPoints<-sample(1:length(useUsChem),n,replace=TRUE)
                                      }
                                      #numComp<-11#numberOfComponentsToUse(useUsFP[calibrationPoints,],useUsChem[calibrationPoints])
                                      #hist(useUsChem[calibrationPoints])
                    #GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
                            ModelRandom<-PLSRFitAndTestNoNSE(useUsFP[calibrationPoints,],useUsChem[calibrationPoints],useUsRealTIme[calibrationPoints],numComp,fitEval,fitFile,fitFileOut,0)
                    #USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
                              PredictRandom<-predict(ModelRandom$Fit,data.matrix(useUsfp),ncomp=numComp,type=c("response"))
                          
                  RNSE[i,awesome]<-NSE(sim=ModelRandom$ObservedAndPredicted[,2],obs=ModelRandom$ObservedAndPredicted[,1])
                  #eventCalibPoints
                  #ENSE[i,awesome]<-NSE(sim=ModelRandom$ObservedAndPredicted[,2],obs=ModelRandom$ObservedAndPredicted[,1])
                  
                  AP_R<-approx(useUs,PredictRandom,n=totmin+1) #create a time series of flow at 1 minute intervals for the whole dataset
                          AP_R$y[(AP_R$y<0)]<-0
                          AP_E<-AP_R$y[AP_R$x>eventStart & AP_R$x<eventStop]
              
                  load_R[i,awesome]<-sum(AP_R$y*flowApprox$y)
                  load_E[i,awesome]<-sum(AP_E*eventFlowFlow)
                  print(i)
                  print(awesome)
                  print(sum(AP_R$y*flowApprox$y))
            }

}
   write.table(load_R,file=paste(chemicals[chemical],"_A.txt",sep=""),
               sep="\t",append = FALSE,
               row.names=FALSE,
               )
  write.table(load_E,file=paste(chemicals[chemical],"_E.txt",sep=""),
            sep="\t",append = FALSE,
            row.names=FALSE,
)
write.table(RNSE,file=paste(chemicals[chemical],"_annualNSE.txt",sep=""),
            sep="\t",append = FALSE,
            row.names=FALSE,
)
#write.table(ENSE,file=paste(chemicals[chemical],"_eventNSE.txt",sep=""),
#            sep="\t",append = FALSE,
#            row.names=FALSE,
#)

}#end for each chemical


#MES_R<-load_R
#MES_E<-load_E
#MES_RNSE<-RNSE