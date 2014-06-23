library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package

#******Specify file paths and names
FPpath<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\FichiersFP\\" #Specify folder where data is located
path<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\"
fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/"  #fiteval_out.txt"
fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
#filename<-c("S2010-2011TurbComp.fp" ,"S2011-2012TurbComp.fp","S2012-2013TurbComp.fp","S2013-2014TurbComp.fp")
filename<-c("S2010-20111stDer.fp", "S2011-20121stDer.fp","S2012-20131stDer.fp","S2013-20141stDer.fp")
#filename<-c("S2010-20111stDerTurbComp.fp" ,"S2011-20121stDerTurbComp.fp","S2012-20131stDerTurbComp.fp","S2013-20141stDerTurbComp.fp")
#LOAD DATA FOR MODEL AND MAKE PLSR MODEL BASED ON WHAT WE KNOW ABOUT WHATEVER
            ModelFilename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")
            Chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
                    "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity");
            bestFitnComps<-c(16,10,12,4,5,20,5,16,5,16,5,5,6,4,3,11,10)
            bestFitFiles<-c(2,2,2,2,3,4,2,1,1,2,2,3,3,2,2,2,2)
            
          #read the data specified by the vector filename
            pTot<-11
            ppo43<-10
            MES<-12
            Turbidity<-17
          startDate<-as.POSIXct("2012-08-01 00:00:00")
          stopDate<-as.POSIXct("2014-08-01 00:00:00")

            myData<-loadDataFile(path,ModelFilename[bestFitFiles[Turbidity]])

            useUsRealTIme<-(myData$realTime[(myData$realTime>startDate&myData$realTime<stopDate)])
            useUsFP<-(myData$fingerPrints[(myData$realTime>startDate&myData$realTime<stopDate),])
            useUsChem<-(myData$ChemData[(myData$realTime>startDate&myData$realTime<stopDate),])

            #numComp<-numberOfComponentsToUse(myData$fingerPrints,myData$ChemData$Turb)
            numComp<-numberOfComponentsToUse(useUsFP,useUsChem$Turb)
            #numComp<-25
            #TurbModel<-PLSRFitAndTest(myData$fingerPrints,myData$ChemData$Turb,myData$realTime,numComp,fitEval,fitFile,fitFileOut,0)
            TurbModel<-PLSRFitAndTest(useUsFP,useUsChem$Turb,useUsRealTIme,numComp,fitEval,fitFile,fitFileOut,0)


#LOAD DATA FROM BRITTANY FOR PROJECTION
            data0<-read.table(file=paste(FPpath,filename[1],sep=""),sep="\t",header=TRUE,skip=1)
            data1<-read.table(file=paste(FPpath,filename[2],sep=""),sep="\t",header=TRUE,skip=1)
            data2<-read.table(file=paste(FPpath,filename[3],sep=""),sep="\t",header=TRUE,skip=1)
            data3<-read.table(file=paste(FPpath,filename[4],sep=""),sep="\t",header=TRUE,skip=1)
      #CRAM ALL THE DATA INTO ONE DATAFRAME            
            data<-rbind(data0,data1,data2,data3)
                    #CLEAN UP A LITTLE
                          remove(data0,data1,data2,data3)

      #parse the date
            realTime<-strptime(data$Date.Time, "%Y.%m.%d  %H:%M:%S")                       #store it back in the matrix
            fingerPrints<-data[,-(dim(data)[2]-3):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
            fingerPrints<-fingerPrints[,-1:-2]                      #remove status and datetime

            useUs<-(realTime[realTime>startDate&realTime<stopDate])
            useUsfp<-fingerPrints[(realTime>startDate&realTime<stopDate) ,]
#USE THE MODEL TO PREDICT OUTPUT FOR THE WHOLE TIME PERIOD
            
            #Predict<-predict(TurbModel$Fit,data.matrix(fingerPrints),ncomp=numComp,type=c("response"))
Predict<-predict(TurbModel$Fit,data.matrix(useUsfp),ncomp=numComp,type=c("response"))
            realTime[realTime<0]<-NA
#             badPoints<-is.na(realTime)
#             goodPoints<-!badPoints
#             timeAxis<-c(realTime[1],realTime[length(realTime)])
#             ylimits<-c(-0.25,2.5)
#             plot(realTime[goodPoints],Predict[goodPoints],ylim=ylimits,type="l")

#WRITE AN OUTPUT FILE WITH REALtIME, AND PREDICTED VALUES
pTotFile<-paste(path,"Turbidity_projected2011.txt",sep="")
Predict<-as.matrix(Predict[,,1])

output<-as.data.frame(matrix(0,ncol=2,nrow=length(useUs)))
output[,1]<-as.POSIXct(useUs,tz="UTC")
output[,2]<-round(Predict,5)

write.table(output,file=pTotFile, append = FALSE,row.names=FALSE,col.names=c("realTime","Predicted"))

plot(output,type="l",col="black",main="projected Turbidity",xlab="date",ylab="concentration",ylim=c(-40,600))

#WRITE ANOTHER FILE WITH REALTIME, OBSERVED AND PREDICTED VALUES
output<-as.data.frame(matrix(0,nrow=length(TurbModel$ObservedAndPredicted[,3]),ncol=3))
output[,1]<-as.POSIXct(TurbModel$ObservedAndPredicted[,3], origin="1970-01-01",tz="")
output[,2]<-round((TurbModel$ObservedAndPredicted[,1]),5)
output[,3]<-round((TurbModel$ObservedAndPredicted[,2]),5)


pTotFile<-paste(path,"Turbidity_observedAndPredicted2011.txt",sep="")
write.table(output,file=pTotFile, sep="\t",append = FALSE,row.names=FALSE,col.names=c("realTime","observed","predicted"))
#Output2010<-output
#Output2011<-output
#Output2012<-output