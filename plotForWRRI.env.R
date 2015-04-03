library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package
library(xts)
library(hydroTSM)
library(signal)

#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Specify file paths and names********************************************************#
#*******************************************************************************************************************/  
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

#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Variable Names for Use in CalibFiles************************************************#
#*******************************************************************************************************************/  
#*******************************************************************************************************************/
#column names in calibration fingerprint files
Chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
        "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity")
#constituent names
chemicals<-c("DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL","NO2","NTOT","NTOTFILT","NH4")#,"SILICA")
#and their position in the calibration file vector
chemN<-c(8,7,12,5,10,11,6,17,1,2,14,15,13,16)#number in chemistry vector


#the type of file to use (1-original 2-1st derivative 3-turb comp 4-1st der turb comp)
#mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)#number corresponding to element of ModelFilename to use
#number of components to include in the model
#NUMCOMP<-c(9,5,5,5,10,10,9,9,8,9,9,8,9)#,8)#component number to use
#changed ptot and ppo43 to 10 components
#a vector of the number of calibration points to use in the PLSR for an iterative model
#calibPoints<-c(13,17,24,36,40,46,52,60,73,91,104)#,122,183,245,365) #NUMBER OF CALIBRATION POINTS


counter<-1
j<-1 #the number of times to run through each model style


#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Output Container Variable Names ****************************************************#
#*******************************************************************************************************************/  
#******************************************************************************************************************/

iD<-matrix(nrow=6000,ncol=30)
colnames(iD)<-c("chem","fileType","ncomp","subRatio",
                "M_n","M_NSE",
                "M_VG","M_G","M_A","M_B",
                "V_VG","V_G","V_A","V_B",
                "C_n","C_r2","C_RMSE","C_NRMSE","C_slope",
                "V_n","V_r2","V_RMSE","V_NRMSE","V_slope",
                "F_n","F_r2","F_RMSE","F_NRMSE","F_slope",
                "subR")

#nash-sutcliffe efficiency output for the iterative model
# RNSE<-matrix(nrow=j,ncol=length(calibPoints))
# ENSE<-matrix(nrow=j,ncol=length(calibPoints))

#cumulative chemical load for the iterative model
# AnnualLoad<-matrix(nrow=j,ncol=length(calibPoints))
# EventLoad<-matrix(nrow=j,ncol=length(calibPoints))

#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Load Specific Files Used for This analysis**************************************#
#*******************************************************************************************************************/  
#*******************************************************************************************************************/

#excludeRows<-c(29,153,155,471,482,486,514,525,530,548.789,816,844)
excludeRows<-c(29,153,155,471,482,486,514,525,530,548,615,616,617,652,658,660,668,671,673,735,790,817,831,845,902,911)
if(!exists("britFuncLoaded",1)){
  source("C:/Users/FBlab/Desktop/Brittany/BrittanyFunctions.R", echo=FALSE)
}
if(!exists("flow",1)){
  source("C:/Users/FBlab/Desktop/Brittany/loadData.R", echo=FALSE)
}




#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Specify Parameters for a Single Run**************************************#
#*******************************************************************************************************************/  
#*******************************************************************************************************************/
startDates<-c(as.POSIXct(paste("2010-10-01 00:00:00",sep=""),tz="UTC"),
              as.POSIXct(paste("2011-10-01 00:00:00",sep=""),tz="UTC"),
              as.POSIXct(paste("2012-10-01 00:00:00",sep=""),tz="UTC"),
              as.POSIXct(paste("2010-10-01 00:00:00",sep=""),tz="UTC")
)


#jan 1 2010
stopDates<-c(as.POSIXct(paste("2011-07-30 00:00:00",sep=""),tz="UTC"),
             as.POSIXct(paste("2012-07-30 00:00:00",sep=""),tz="UTC"),
             as.POSIXct(paste("2013-07-30 00:00:00",sep=""),tz="UTC"),
             as.POSIXct(paste("2014-07-30 00:00:00",sep=""),tz="UTC")
)


#june 15 2011  




#eventStart<-as.POSIXct("2010-12-01 07:13:38",tz="UTC")
#eventStop<-as.POSIXct("2010-12-15 04:33:31",tz="UTC")



#chemical<-2   #chemical chemicals<-c("DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL",                                    #"NO2","NTOT","NTOTFILT","NH4")#,"SILICA")

Types<-c("original","prunedO",
         "1stDer","pruned1D",
         "turbComp","prunedTC",
         "1stDerTurbComp","prunedTC1D")  
#come back and change




makePlotsForWWRI<-function(year=1,FileType=1,keepSkippedDates=0,chemical=4,numComp=3,Flow=Flow,EXP=2){
  k<-FileType
  y<-year
  ksd<-keepSkippedDates #keep skipped dates
  #chemical<-14 #select Chem
  #numComp<-4 #number of comps
  
  
  
  
  if (y==1){
    startDates<-c(as.POSIXct(paste("2010-11-18 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2011-04-25 00:00:00",sep=""),tz="UTC")
    )
    stopDates<-c(as.POSIXct(paste("2011-03-01 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2011-05-20 00:00:00",sep=""),tz="UTC")
    )
  }
  if (y==2){
    startDates<-c(as.POSIXct(paste("2011-10-18 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2012-3-2 00:00:00",sep=""),tz="UTC")
    )
    stopDates<-c(as.POSIXct(paste("2012-2-8 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2012-06-01 00:00:00",sep=""),tz="UTC")
    )
  }
  
  if (y==3){
    startDates<-c(as.POSIXct(paste("2012-10-18 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2013-3-2 00:00:00",sep=""),tz="UTC")
    )
    stopDates<-c(as.POSIXct(paste("2013-3-2 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2013-07-01 00:00:00",sep=""),tz="UTC")
    )
    
  }
  if (y==4){
    startDates<-c(as.POSIXct(paste("2010-11-18 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2011-04-25 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2011-10-18 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2012-3-2 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2012-10-18 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2013-3-2 00:00:00",sep=""),tz="UTC")           
    )
    stopDates<-c(as.POSIXct(paste("2011-03-01 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2011-05-20 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2012-2-8 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2012-06-01 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2013-3-2 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2013-07-01 00:00:00",sep=""),tz="UTC")
                 
    ) 
  }
  startDate<-min(startDates)
  stopDate<-max(stopDates)
  dates<-cbind(startDates,stopDates)
  
  if(FileType%%2){
    dates<-cbind(min(startDates),max(stopDates))
  }
  
  
  
  
  flow<-subsetAll(flow=Flow,dateWindows=cbind(startDates,stopDates))$flow
  fileType<-Types[k]   #specify file TYPE 1-REGFP 2-1STDER 3-TURBCOMP 4-1STDERTURBCOMP 
  subsetRate=0
  
  
  #change here to send vector of startDates and stopDate
  calibration<-selectSpecData(fileType,"calibration")
  specDataToModel<-selectSpecData(fileType,"fingerPrints")
  
  all<-subsetAll(calibData=calibration,modelData=specDataToModel,
                 chem=chemN[chemical],
                 dateWindows=dates,
                 keepSpecData=ksd)
  
  #subsetAll(calibData=0,modelData=0,flow=0,chem,dateWindows,keepSpecData=0)
  
  
  #returns a list
  #   output$modelQuality a named vector of summary characteristics for the model and its verification, 
  #     including NSE, r2 and fitEval for the model and, if applicable validation datasets
  #   output$PredictedConcentrations
  fitFile<-paste(fitPath,  
                 Chem[chemN[chemical]],
                 "_",
                 fileType,
                 "_",
                 numComp,
                 "_",
                 paste((as.POSIXlt(startDate))$year+1900,"-",(as.POSIXlt(stopDate))$year+1900,sep=""),
                 "PLSR.in",sep=""
  )
  fitFileOut<-paste(fitPath, 
                    Chem[chemN[chemical]],
                    "_",
                    fileType,
                    "_",
                    numComp,
                    "_",
                    paste((as.POSIXlt(startDate))$year+1900,"-",(as.POSIXlt(stopDate))$year+1900,sep=""),
                    "PLSR_out.txt",sep=""
  )
  modelOutput<-modelExecution(numComp=numComp,subsampleRate=subsetRate,
                              calibData=all$calibData,
                              dataToModel=all$specData,
                              fitEval,fitFile,fitFileOut)
  
  #make plots
  #iD[counter,30]<-Chem[chemN[chemical]]
  png(file = paste(OutputPath,
                   Chem[chemN[chemical]],
                   "_",
                   fileType,
                   "_",
                   numComp,
                   "_",
                   (as.POSIXlt(startDate))$year+1900,
                   "_1.png",sep=""),
      width=1000,height=500
      
  )
  
  plot(as.POSIXct(flow$realTime,origin="1970-01-01 00:00:00"),
       flow$flow,col="black",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1])),
  )
  
  par(new=TRUE)
  plotMax<-max(modelOutput$ObservedandPredicted[,1])*1.5
  plotMin<-min((min(modelOutput$ObservedandPredicted[,1])-min(modelOutput$ObservedandPredicted[,1])*0.1),0)
  
  plot(as.POSIXct(modelOutput$ObservedandPredicted[,3],origin="1970-01-01 00:00:00"),
       modelOutput$ObservedandPredicted[,1],
       main=chemicals[chemical],
       ylim=c(plotMin,plotMax),
       xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE)),
       xlab="date",
       ylab=paste(chemicals[chemical],"concentration",sep=" "),
       cex.lab=1.5,
       col="green",cex=EXP,pch=18)
  
  mtext(paste(fileType,numComp,startDate,sep="   "))
  abline(h=0,col="grey")
  dev.off()
  
  
  png(file = paste(OutputPath,
                   Chem[chemN[chemical]],
                   "_",
                   fileType,
                   "_",
                   numComp,
                   "_",
                   (as.POSIXlt(startDate))$year+1900,
                   "_2.png",sep=""),
      width=1000,height=500
      
  )
  
  
  plot(as.POSIXct(flow$realTime,origin="1970-01-01 00:00:00"),
       flow$flow,col="black",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1])),
  )
  
  par(new=TRUE)
  plotMax<-max(modelOutput$ObservedandPredicted[,1])*1.5
  plotMin<-min((min(modelOutput$ObservedandPredicted[,1])-min(modelOutput$ObservedandPredicted[,1])*0.1),0)
  
  plot(as.POSIXct(modelOutput$ObservedandPredicted[,3],origin="1970-01-01 00:00:00"),
       modelOutput$ObservedandPredicted[,1],
       main=chemicals[chemical],
       ylim=c(plotMin,plotMax),
       xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE)),
       xlab="date",
       ylab=paste(chemicals[chemical],"concentration",sep=" "),cex.lab=1.5,
       col="green",cex=EXP,pch=18)
  
  mtext(paste(fileType,numComp,startDate,sep="   "))
  
  
  points(modelOutput$ObservedandPredicted[,3],
         modelOutput$ObservedandPredicted[,2],
         col="blue",cex=EXP,pch=20)
  abline(h=0,col="grey")
  
  
  
  
  dev.off()
  
  
  
  
  
  png(file = paste(OutputPath,
                   Chem[chemN[chemical]],
                   "_",
                   fileType,
                   "_",
                   numComp,
                   "_",
                   (as.POSIXlt(startDate))$year+1900,
                   "_3.png",sep=""),
      width=1000,height=500
      
  )
  
  
  
  plot(as.POSIXct(flow$realTime,origin="1970-01-01 00:00:00"),
       flow$flow,col="black",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1])),
  )
  
  par(new=TRUE)
  plotMax<-max(modelOutput$ObservedandPredicted[,1])*1.5
  plotMin<-min((min(modelOutput$ObservedandPredicted[,1])-min(modelOutput$ObservedandPredicted[,1])*0.1),0)
  
  plot(as.POSIXct(modelOutput$ObservedandPredicted[,3],origin="1970-01-01 00:00:00"),
       modelOutput$ObservedandPredicted[,1],
       main=chemicals[chemical],
       ylim=c(plotMin,plotMax),
       xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE)),
       xlab="date",
       ylab=paste(chemicals[chemical],"concentration",sep=" "),
       cex.lab=1.5,
       col="green",cex=EXP,pch=18)
  
  mtext(paste(fileType,numComp,startDate,sep="   "))
  
  
  points(modelOutput$ObservedandPredicted[,3],
         modelOutput$ObservedandPredicted[,2],
         col="blue",cex=EXP,pch=20)
  abline(h=0,col="grey")
  
  points(as.POSIXct(modelOutput$PredictedConcentrations[,1],origin="1970-01-01 00:00:00"),
         modelOutput$PredictedConcentrations[,2],
         col=rainbow(10)[2],
         ylim=c(plotMin,plotMax),
         xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE)),
         type="l"
  )
  points(as.POSIXct(modelOutput$PredictedConcentrations[,1],origin="1970-01-01 00:00:00"),
         modelOutput$PredictedConcentrations[,2],
         col=rainbow(10)[2],
         ylim=c(plotMin,plotMax),
         xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE)),
         type="l"
  )
  
  dev.off()
  
}
