library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package
source('~/Documents/GitHub/Brittany/BrittanyFunctions.R')

#******Specify file paths and names
FPpath<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\lkWhlr\\" #Specify folder where data is located
path<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\"
fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/"  #fiteval_out.txt"
fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
#filename<-c("S2010-2011TurbComp.fp" ,"S2011-2012TurbComp.fp","S2012-2013TurbComp.fp","S2013-2014TurbComp.fp")
#filename<-c("S2010-20111stDer.fp", "S2011-20121stDer.fp","S2012-20131stDer.fp","S2013-20141stDer.fp")
filename<-"LkWhlrFP.csv"
#filename<-c("S2010-20111stDerTurbComp.fp" ,"S2011-20121stDerTurbComp.fp","S2012-20131stDerTurbComp.fp","S2013-20141stDerTurbComp.fp")
#LOAD DATA FOR MODEL AND MAKE PLSR MODEL BASED ON WHAT WE KNOW ABOUT WHATEVER
            ModelFilename<-c("LkWhlrTKN_Cal.csv")
            Chem<-c("TKN");
#            bestFitnComps<-c(16,10,12,4,5,20,5,16,5,16,5,5,6,4,3,11,10)
 #           bestFitFiles<-c(2,2,2,2,3,4,2,1,1,3,3,3,3,2,2,2,2)
            
          #read the data specified by the vector filename
data<-read.table(file=paste(FPpath,ModelFilename,sep=""),sep=",",header=TRUE,skip=0)
Date<-substr(data$Date.Time, 1, 10)        #the date
T<-substr(data$Date.Time, 12, 19)       #the time
T[T==""]="00:00:00"                    #ah, but where it should be midnight, often it was " "
D<-paste(Date,T,sep=" ")               #put back in the fixed values
D<-strptime(D, '%Y.%m.%d %H:%M:%S')    #convert it to a _flawless_ time opject
data$DateScan<-D                       #store it back in the matrix

  realTime<-(data$DateScan)                                       #pull out time
            ChemData<-data[,(dim(data)[2]-2)]            #pull out the analyte data
            fingerPrints<-data[,-(dim(data)[2]-7):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
            fingerPrints<-fingerPrints[,-1:-2]  
  numComp<-numberOfComponentsToUse(fingerPrints,ChemData)
            print(numComp)
            #numComp<-25  
            numComp<-9
TKNModel<-PLSRFitAndTest(fingerPrints,ChemData,realTime,numComp,fitEval,fitFile,fitFileOut,0)
 

fpData<-read.table(file=paste(FPpath,"fps.txt",sep=""),sep="\t",header=TRUE,skip=1)
            Date<-substr(fpData$Date.Time, 1, 10)        #the date
            T<-substr(fpData$Date.Time, 12, 19)       #the time
            D<-paste(Date,T,sep=" ")               #put back in the fixed values
           
D<-strptime(D, '%Y.%m.%d %H:%M:%S')    #convert it to a _flawless_ time opject

            fpData$DateScan<-D                       #store it back in the matrix

realTimeExp<-(fpData$DateScan)   
        fingerPrints<-fpData[,-(dim(fpData)[2]-5):-(dim(fpData)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
        fingerPrints<-fingerPrints[,-1:-2]  



            Predict<-predict(TKNModel$Fit,data.matrix(fingerPrints),ncomp=numComp,type=c("response"))
            realTimeExp[realTimeExp<0]<-NA


#plot(realTimeExp[2:200],Predict[2:200])


model<-as.data.frame(matrix(0,nrow=length(TKNModel$ObservedAndPredicted[,3]),ncol=3))
model[,1]<-as.POSIXct(TKNModel$ObservedAndPredicted[,3], origin="1970-01-01",tz="")
model[,2]<-round((TKNModel$ObservedAndPredicted[,1]),5)
model[,3]<-round((TKNModel$ObservedAndPredicted[,2]),5)

plot(realTimeExp[2:100],Predict[2:100],col="blue",pch=19)
points(model[,1],model[,2],col="red")


outFile<-paste(path,"TKNObservedAndPredicted.txt",sep="")

write.table(model,file=outFile, sep="\t",append = FALSE,row.names=FALSE,col.names=c("realTime","observed","predicted"))

projection<-as.data.frame(matrix(0,ncol=2,nrow=length(realTimeExp)))
projection[,1]<-as.POSIXct(realTimeExp,tz="")
projection[,2]<-round(Predict,2)

outFile<-paste(path,"TKNPredicted.txt",sep="")
write.table(projection,file=outFile, sep="\t",append = FALSE,row.names=FALSE,col.names=c("realTime","predicted"))
