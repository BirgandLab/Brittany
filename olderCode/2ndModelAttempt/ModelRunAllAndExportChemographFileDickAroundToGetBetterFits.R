library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package
library(xts)
library(hydroTSM)
library(signal)

#******Specify file paths and names
#FPpath<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\FichiersFP\\" #Specify folder where data is located
#path<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\"        #upstairs box
#fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/"  #fiteval_out.txt location on upstairs box

fitPath<-"C:/Users/FBlab/Desktop/FITEVAL2_win/FITEVAL2_win/"    #location on downstairs box
path<-"C:/Users/FBlab/Desktop/work_here/Data/"   #downstairs box
FPpath<-"C:/Users/FBlab/Desktop/work_here/Data/FichiersFP/" #downstairs box

fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")


ModelFilename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")
Chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
        "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity");
chemicals<-c("DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL","NO2","NTOT","NTOTFILT","NH4")#,"SILICA")
#mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)#number corresponding to element of ModelFilename to use
chemN<-c(8,7,12,5,10,11,6,17,1,2,14,15,13)#,16)#number in chemistry vector
#NUMCOMP<-c(9,5,5,5,10,10,9,9,8,9,9,8,9)#,8)#component number to use
#changed ptot and ppo43 to 10 components
#calibPoints<-c(13,17,24,36,40,46,52,60,73,91,104)#,122,183,245,365) #NUMBER OF CALIBRATION POINTS
counter<-1
iD<-matrix(nrow=100,ncol=10)
colnames(iD)<-c("n","NSE","VG","G","A","B","chemical","numComp","fool","o")


#this section is out of order and has no value to initialized these variables
chemoGraphicTable<-matrix(nrow=length(flowApprox$x), ncol=length(chemicals)+2)
chemoGraphicTable[,1]<-flowApprox$x  #date for flow calculation
chemoGraphicTable[,2]<-flowApprox$y  #flow calculation
colnames(chemoGraphicTable)<-c("date","flow(l*s-1","DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL","NO2","NTOT","NTOTFILT","NH4")
j<-1 #the number of times to run through each model style
RNSE<-matrix(nrow=j,ncol=length(calibPoints))
ENSE<-matrix(nrow=j,ncol=length(calibPoints))
load_R<-matrix(nrow=j,ncol=length(calibPoints))
load_E<-matrix(nrow=j,ncol=length(calibPoints))



source("C:/Users/FBlab/Documents/GitHub/Brittany/BrittanyFunctions.R", echo=FALSE)
source("C:/Users/FBlab/Documents/GitHub/Brittany/loadData.R", echo=TRUE)



  # mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)#number corresponding to element of ModelFilename to use
  # chemN<-c(8,7,12,5,10,11,6,17,1,2,14,15,13)#,16)#number in chemistry vector
  #  NUMCOMP<-c(9,5,5,5,10,10,9,9,8,9,9,8,9)#,8)#component number to use
  #here specify chemical
 


  chemical<-6 #chemical chemicals<-c("DIC","DOC","MES","NO3","PPO43","PTOT","SO4","TURB","CL",
                                    #"NO2","NTOT","NTOTFILT","NH4")#,"SILICA")
  numComp<-5 #specify         
  fool<-1 #specify file TYPE 1-REGFP 2-1STDER 3-TURBCOMP 4-1STDERTURBCOMP 
                  #   #mfN<-c(1,2,4,2,2,3,4,2,2,2,3,2,3)#,2)
  o<-0       #points to drop from the front of a calibration
  s<-0     #points to drop from the back end
iD[counter,8]<-numComp
iD[counter,7]<-chemical
iD[counter,9]<-fool
iD[counter,10]<-o


source("C:/Users/FBlab/Documents/GitHub/Brittany/subsetData.R", echo=FALSE)
       # n<-length(useUsChem)
          
#GENERATE A PLSR MODEL FOR THE SELECTED CALIBRATION DATA
                            ModelRandom<-PLSRFitAndTestNoNSE(useUsFP[o:(length(useUsChem)-s),],
                                                             useUsChem[o:(length(useUsChem)-s)],
                                                             useUsRealTIme[o:(length(useUsChem)-s)],
                                                             numComp,fitEval,fitFile,fitFileOut,0)
#USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
                              PredictRandom<-predict(ModelRandom$Fit,
                                                     data.matrix(useUsfp),
                                                     ncomp=numComp,
                                                     type=c("response"))#not the subtle but meaningful difference in 
                                                       #useUsFP and use useUsfp. one had422 obs the other 17247 obs
                          
  
                  AP_R<-approx(useUs,PredictRandom,n=totmin+1) #create a time series of [ ] at 1 minute intervals for the whole dataset
                         # AP_R$y[(AP_R$y<0)]<-0
                  event<-(AP_R$x>as.numeric(eventStart)&AP_R$x<as.numeric(eventStop))
                  chemoGraphicTable[,chemical+2]<-AP_R$y
                  iD[counter,1]<-ModelRandom$Stats[1]
                  iD[counter,2]<-NSE(sim=ModelRandom$ObservedAndPredicted[,2],obs=ModelRandom$ObservedAndPredicted[,1])
                  iD[counter,3:6]<-OB(ModelRandom$ObservedAndPredicted[,1],ModelRandom$ObservedAndPredicted[,2],fitEval,fitFile,fitFileOut)
                  #print(i)

source("C:/Users/FBlab/Documents/GitHub/Brittany/makePlots.R", echo=FALSE)

counter<-counter+1





