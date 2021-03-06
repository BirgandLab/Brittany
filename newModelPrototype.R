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
    chemN<-c(8,7,12,5,10,11,6,17,1,2,14,15,13)#,16)#number in chemistry vector


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

  iD<-matrix(nrow=2000,ncol=30)
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

Types<-c("original","1stDer","turbComp","1stDerTurbComp")  
#come back and change
#the chemical number to 1:9
#and the dates
for(y in 1:length(startDates)){
  
          startDate<-startDates[y]
          stopDate<-stopDates[y]
          flow<-subset(Flow,"flow",startDate,stopDate)
        for(chemical in 1:9){
                  for(numComp in 3:11){
                  #numComp<-4    #specify number of components to use in the model
                      for(k in 1:length(Types)){
                          fileType<-Types[k]   #specify file TYPE 1-REGFP 2-1STDER 3-TURBCOMP 4-1STDERTURBCOMP 
                          subsetRate=0
                                          #   #mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)
                          iD[counter,1]<-paste((as.POSIXlt(startDate))$year+1900,"-",(as.POSIXlt(stopDate))$year+1900,sep="")
                          iD[counter,2]<-fileType
                          iD[counter,3]<-numComp
                          iD[counter,4]<-subsetRate
                        
                        calibration<-subsetSpecData(fileType,"calibration",startDate,stopDate,chemN[chemical])
                        specDataToModel<-subsetSpecData(fileType,"fingerPrints",startDate,stopDate)
                        
                     
                                         
                            for(i in 1:1){
                                  #returns a list
                                  #   output$modelQuality a named vector of summary characteristics for the model and its verification, 
                                  #     including NSE, r2 and fitEval for the model and, if applicable validation datasets
                                  #   output$PredictedConcentrations
                                  modelOutput<-modelExecution(chemical,numComp,subsetRate,calibration,specDataToModel,fitEval,fitFile,fitFileOut)
                                  
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
                                                          "_",
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
}#end for each year


#source("C:/Users/FBlab/Desktop/Brittany/makePlots.R", echo=FALSE)






