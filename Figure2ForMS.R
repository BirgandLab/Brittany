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

  iD<-matrix(nrow=100,ncol=30)
    #colnames(iD)<-c("n","NSE","VG","G","A","B","chemical","numComp","fool","o")

  
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

source("C:/Users/FBlab/Desktop/Brittany/BrittanyFunctions.R", echo=FALSE)
source("C:/Users/FBlab/Desktop/Brittany/loadData.R", echo=FALSE)



#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Specify Parameters for a Single Run**************************************#
#*******************************************************************************************************************/  
#*******************************************************************************************************************/
startDate<-as.POSIXct(paste("2010-10-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
stopDate<-as.POSIXct(paste("2011-07-30 00:00:00",sep=""),tz="UTC")    #june 15 2011  


  chemical<-6   #chemical chemicals<-c("DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL",                                    #"NO2","NTOT","NTOTFILT","NH4")#,"SILICA")
  numComp<-4    #specify number of components to use in the model
  fileType<-"original"   #specify file TYPE 1-REGFP 2-1STDER 3-TURBCOMP 4-1STDERTURBCOMP 
  subsetRate=0.5
                  #   #mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)


#make a record of the parameters used
iD[counter,11:14]<-c(numComp,chemical,fileType,0)
calibration<-subsetSpecData(fileType,"calibration",startDate,stopDate,chemN[chemical])

specDataToModel<-subsetSpecData(fileType,"fingerPrints",startDate,stopDate)
flow<-subset(flow,"flow",startDate,stopDate)
                 
for(i in 1:1){
#GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
                          
                   ModelConcentration<-PLSRFitAndTestNoNSE(calibration$fingerPrints,
                                                   calibration$ChemData,
                                                   calibration$realTime,
                                                   numComp,fitEval,fitFile,fitFileOut,subsetRate)                  
#USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
                  PredictedConcentration<-predict(ModelConcentration$Fit,
                                         as.matrix(specDataToModel$fingerPrints),
                                         ncomp=numComp,
                                         type=c("response"))  
#**************************************************************
#**************************************************************
#have updated things up to here--now need to do just a little bit more to smooth the rest of this over, and 
    #perhaps there will be a functioning tool out of the whole thing.

                    PredictedAnnualConcentrationTS<-as.data.frame(approx(subsetRealTime,PredictedConcentration,n=totmin+1)) #create a time series of [ ] at 1 minute intervals for the whole dataset
                         # PredictedAnnualConcentrationTS$y[(PredictedAnnualConcentratinTS$y<0)]<-0
                    PredictedEventConcentrationTS<-(PredictedAnnualConcentrationTS$x>as.numeric(eventStart)&PredictedAnnualConcentrationTS$x<as.numeric(eventStop))
                        
                        if(!exists("chemoGraphicTable",1)){
                          chemoGraphicTable<-matrix(nrow=length(AnnualFlowApprox$x), ncol=length(chemicals)+2)
                          chemoGraphicTable[,1]<-AnnualFlowApprox$x  #date for flow calculation
                          chemoGraphicTable[,2]<-AnnualFlowApprox$y  #flow calculation
                          colnames(chemoGraphicTable)<-c("date","flow(l*s-1","DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL","NO2","NTOT","NTOTFILT","NH4")
                          #make base plot
                               plot(as.POSIXct(PredictedAnnualConcentrationTS$x,origin="1970-01-01 00:00:00"),
                                    chemoGraphicTable[,2],col="green",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
                                    xlim=c(min(PredictedAnnualConcentrationTS$x),max(PredictedAnnualConcentrationTS$x))
                                   )
                               axis(4)
                               par(new=TRUE)  
                        }

                  chemoGraphicTable[,chemical+2]<-PredictedAnnualConcentrationTS$y
                  
                  iD[counter,1]<-ModelConcentration$Stats[1]
                  if(subsetRate>0){
                        iD[counter,2]<-NSE(sim=ModelConcentration$OaP1[,2],obs=ModelConcentration$OaP1[,1])
                        iD[counter,3:6]<-OB(ModelConcentration$OaP1[,1],ModelConcentration$OaP1[,2],fitEval,fitFile,fitFileOut)
                        iD[counter,7:10]<-OB(ModelConcentration$OaP2[,1],ModelConcentration$OaP2[,2],fitEval,fitFile,fitFileOut)
                        }
                  if(subsetRate==0){
                        iD[counter,2]<-NSE(sim=ModelConcentration$ObservedAndPredicted[,2],obs=ModelConcentration$ObservedAndPredicted[,1])
                        iD[counter,3:6]<-OB(ModelConcentration$ObservedAndPredicted[,1],ModelConcentration$ObservedAndPredicted[,2],fitEval,fitFile,fitFileOut)
                  }
                  iD[counter,15:29]<-ModelConcentration$Stats
                  iD[counter,30]<-subsetRate
if(i==1){
plot(as.POSIXct(PredictedAnnualConcentrationTS$x,origin="1970-01-01 00:00:00"),
     PredictedAnnualConcentrationTS$y,type="l",col=rainbow(10)[i],
     main=chemicals[chemical],
     ylim=c(min(PredictedAnnualConcentrationTS$y),max(PredictedAnnualConcentrationTS$y)),
     xlim=c(min(PredictedAnnualConcentrationTS$x),max(PredictedAnnualConcentrationTS$x)),
     xlab="date",
     ylab=paste(chemicals[chemical],"concentration",sep=" ")
     )
}
else{
  points(as.POSIXct(PredictedAnnualConcentrationTS$x,origin="1970-01-01 00:00:00"),
          PredictedAnnualConcentrationTS$y,type="l",col=rainbow(10)[i]
        )
}

counter<-counter+1

}

#source("C:/Users/FBlab/Desktop/Brittany/makePlots.R", echo=FALSE)






