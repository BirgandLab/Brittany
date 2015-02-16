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

#this section is out of order and has no value to initialized these variables
  
#nash-sutcliffe efficiency output for the iterative model
  RNSE<-matrix(nrow=j,ncol=length(calibPoints))
  ENSE<-matrix(nrow=j,ncol=length(calibPoints))

#cumulative chemical load for the iterative model
  AnnualLoad<-matrix(nrow=j,ncol=length(calibPoints))
  EventLoad<-matrix(nrow=j,ncol=length(calibPoints))

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
  

  chemical<-6 #chemical chemicals<-c("DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL",
                                    #"NO2","NTOT","NTOTFILT","NH4")#,"SILICA")
  numComp<-4 #specify number of components to use in the model

  fileType<-1 #specify file TYPE 1-REGFP 2-1STDER 3-TURBCOMP 4-1STDERTURBCOMP 
                  #   #mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)
  o<-50      #points to drop from the front of a calibration
  s<-0     #points to drop from the back end
  
  #make a record of the parameters used
  iD[counter,11:14]<-c(numComp,chemical,fileType,o)

  #subset the data--perhaps rewrite this as a function that is passed a vector of dates to include
  source("C:/Users/FBlab/Desktop/Brittany/subsetData.R", echo=FALSE)
                 
  

for(i in 1:1){
#GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
#                            ModelRandom<-PLSRFitAndTestNoNSE(calibrationFingerPrints[o:(length(calibrationAnalytes)-s),],
#                                                            calibrationAnalytes[o:(length(calibrationAnalytes)-s)],
#                                                            calibrationRealTime[o:(length(calibrationAnalytes)-s)],
#                                                            numComp,fitEval,fitFile,fitFileOut,0.9)         
#                           
                            ModelRandom<-PLSRFitAndTestNoNSE(calibrationFingerPrints2[o:(length(calibrationAnalytes2)-s),],
                                                             calibrationAnalytes2[o:(length(calibrationAnalytes2)-s)],
                                                             calibrationRealTime2[o:(length(calibrationAnalytes2)-s)],
                                                             numComp,fitEval,fitFile,fitFileOut,0.9)                  
#USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
                              PredictRandom<-predict(ModelRandom$calibrationFit,
                                                     data.matrix(subsetFingerPrints),
                                                     ncomp=numComp,
                                                     type=c("response"))  
                                                         
#                              PredictRandom<-predict(ModelRandom$Fit,
#                                                     data.matrix(useUsfp),
#                                                     ncomp=numComp,
#                                                     type=c("response"))#not the subtle but meaningful difference in 
#                              #useUsFP and use useUsfp. one had422 obs the other 17247 obs


                  AP_R<-approx(subsetRealTime,PredictRandom,n=totmin+1) #create a time series of [ ] at 1 minute intervals for the whole dataset
                         # AP_R$y[(AP_R$y<0)]<-0
                  event<-(AP_R$x>as.numeric(eventStart)&AP_R$x<as.numeric(eventStop))
                        
                        if(!exists("chemoGraphicTable",1)){
                          chemoGraphicTable<-matrix(nrow=length(flowApprox$x), ncol=length(chemicals)+2)
                          chemoGraphicTable[,1]<-flowApprox$x  #date for flow calculation
                          chemoGraphicTable[,2]<-flowApprox$y  #flow calculation
                          colnames(chemoGraphicTable)<-c("date","flow(l*s-1","DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL","NO2","NTOT","NTOTFILT","NH4")
                          #make base plot
                               plot(as.POSIXct(AP_R$x,origin="1970-01-01 00:00:00"),
                                    chemoGraphicTable[,2],col="green",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
                                    xlim=c(min( AP_R$x),max( AP_R$x))
                                   )
                               axis(4)
                               par(new=TRUE)  
                        }

                  chemoGraphicTable[,chemical+2]<-AP_R$y
                  
                  iD[counter,1]<-ModelRandom$Stats[1]
                  iD[counter,2]<-NSE(sim=ModelRandom$OaP1[,2],obs=ModelRandom$OaP1[,1])
                  iD[counter,3:6]<-OB(ModelRandom$OaP1[,1],ModelRandom$OaP1[,2],fitEval,fitFile,fitFileOut)
                  iD[counter,7:10]<-OB(ModelRandom$OaP2[,1],ModelRandom$OaP2[,2],fitEval,fitFile,fitFileOut)
                  iD[counter,15:29]<-ModelRandom$Stats
if(i==1){
plot(as.POSIXct(AP_R$x,origin="1970-01-01 00:00:00"),
     AP_R$y,type="l",col=rainbow(10)[i],
     main=chemicals[chemical],
     ylim=c(min(AP_R$y),max( AP_R$y)),
     xlim=c(min( AP_R$x),max( AP_R$x)),
     xlab="date",
     ylab=paste(chemicals[chemical],"concentration",sep=" ")
)
}
else{
  points(as.POSIXct(AP_R$x,origin="1970-01-01 00:00:00"),
          AP_R$y,type="l",col=rainbow(10)[i]
)
}

                  #print(i)


counter<-counter+1

}
source("C:/Users/FBlab/Documents/GitHub/Brittany/makePlots.R", echo=FALSE)






