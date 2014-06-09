library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package

#******Specify file paths and names
inPath<-"C:/Users/FBLab/Documents/GitHub/Brittany/inputFiles/" #Specify folder where data is located
outPath<-"C:/Users/FBLab/Documents/GitHub/Brittany/output/"
fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/" #fiteval_out.txt"
#source<-

fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
filename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")

#store names for the lab analytes
Chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
        "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity");
bestFitnComps<-c(16,10,12,4,5,20,5,16,5,16,5,5,6,4,3,11,10)
bestFitFiles<-c(2,2,2,2,3,4,2,1,1,3,3,3,3,2,2,2,2)


fitQuality<-matrix(nrow=3400,ncol=14)
colnames(fitQuality)<-c("calibrationVG","calibrationG","calibratinA","calibrationB","calibFitVG","calibFitG","calibFitA","calibFitB","allDataVG","allDataG","allDataA","allDataB","i","j")
Stats<-matrix(nrow=3400,ncol=17)
colnames(Stats)<-c("n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","i","j")

#read the data specified by the vector filename
counter<-1
Components<-matrix(nrow=70,ncol=4)

for (reps in 1:200)
#load data for checking number of components on
for (chem in 1:17){
  Comps<-bestFitnComps[chem]
  fn<-bestFitFiles[chem]
  
  
  #all chem analytes are listed in different columns 
  #1 CL     2 NO2   3 NNO2      4 NO3      5 NNO3    6 SO4	
  #7 DOC	  8 DIC	  9 UV254	    10 PPO43	11 Ptot	  12 MES	
  #13 NNH4  14 Ntot	15 NTotFilt	16 Silic	17 Turb
  
 # jpeg(file=paste(outPath,".",Chem[chem],"best.jpg",sep=""))
            myData<-loadDataFile(inPath,filename[fn])
          #data are returned in a list myData$fingerPrints
          #                            myData$realTime
          #                            myData$ChemData
        
        #pull out the column of chem data of interest
            ChemConc<-as.matrix(myData$ChemData[,chem])                       #take out the only one we care about
                      fp<-cbind(myData$fingerPrints,myData$realTime,as.matrix(ChemConc))      #bind the two matrices together to determine complete cases
                      fp<-fp[complete.cases(fp[,2:dim(fp)[2]]),]              #just keep the fingerprints that have analytical data
            
                      ChemConc<-as.matrix(fp[,dim(fp)[2]])                    #pull chem back out
                      fp<-fp[,-dim(fp)[2]]                                    #pull off the fingerprint
        
                      realTime<-as.matrix(fp[,dim(fp)[2]])                    #pull real time off (now the last column)
                      fp<-fp[,-dim(fp)[2]]     
        #pop it off the end of the dataframe
        Comps<-bestFitnComps[chem]
        fn<-bestFitFiles[chem]
        #Comps<-numberOfComponentsToUse(fp,ChemConc)
  #Calculate Nash-Sutcliffe based goodness of fit
  fitFile<-paste(fitPath,Chem[chem],"2","PLSR.in",sep="")
  fitFileOut<-paste(fitPath,Chem[chem],"2","PLSR_out.txt",sep="")
  Output<-PLSRFitAndTest(fp,ChemConc,realTime,Comps,fitEval,fitFile,fitFileOut,0.25)
  fitQuality[counter,]<-c(Output$fitQuality,Comps,chem)
  Stats[counter,]<-c(Output$Stats,Comps,chem)
 counter<-counter+1
 #dev.off()
} #for each chemical