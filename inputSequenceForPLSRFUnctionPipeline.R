#load libraries for analysis
library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package

#******Specify file paths and names
path<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\" #Specify folder where data is located
fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/"#fiteval_out.txt"
fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
filename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")

#store names for the lab analytes
chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
        "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity");


#for each chemical constituent, 
#load the appropriate file, access the right #of plsr model comopnents
#run the evaluation model


#read the data specified by the vector filename
data<-read.table(file=paste(path,filename[4],sep=""),sep=",",header=TRUE,skip=0)

#parse the date
      Date<-substr(data$DateScan, 1, 10)        #the date
         T<-substr(data$DateScan, 12, 19)       #the time
         T[T==""]="00:00:00"                    #ah, but where it should be midnight, often it was " "
         D<-paste(Date,T,sep=" ")               #put back in the fixed values
         D<-strptime(D, '%d/%m/%Y %H:%M:%S')    #convert it to a _flawless_ time opject
         data$DateScan<-D                       #store it back in the matrix

realTime<-(data$DateScan)                                       #pull out time
ChemData<<-data[,(dim(data)[2]-16):(dim(data)[2])]            #pull out the analyte data
allFingerPrints<-data[,-(dim(data)[2]-20):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
allFingerPrints<-allFingerPrints[,-1:-2]                      #remove status and datetime

remove(data,D,Date)                                           #clean up


#now all chem analyses are listed in different columns 
#1 CL     2 NO2   3 NNO2      4 NO3	    5 NNO3	  6 SO4	
#7 DOC	  8 DIC	  9 UV254	    10 PPO43	11 Ptot	  12 MES	
#13 NNH4  14 Ntot	15 NTotFilt	16 Silic	17 Turb
ChemConc<-as.matrix(ChemData[,11]) #put back on the only one we care about

Components<-numberOfComponentsToUse(allFingerPrints,ChemConc)

subSample<-0.5
cMin<-5
cMax<-Components
numReps<-1000
totalIterations<-numReps
fitQuality<-matrix(nrow=totalIterations,ncol=14)
colnames(fitQuality)<-c("calibrationVG","calibrationG","calibratinA","calibrationB","calibFitVG","calibFitG","calibFitA","calibFitB","allDataVG","allDataG","allDataA","allDataB","i","j")
Stats<-matrix(nrow=totalIterations,ncol=17)
colnames(Stats)<-c("n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","i","j")

counter<-1
#for(numComponents in seq(cMin:cMax)){
    for(i in 1:numReps){
        testOutput<-PLSRFitAndTest(allFingerPrints,ChemConc,realTime,Components,fitEval,fitFile,fitFileOut,subSample)
        fitQuality[i,]<-c(testOutput$fitQuality,numComponents,i)
        Stats[i,]<-c(testOutput$Stats,numComponents,i)
      counter<-counter+1
    }
#}

