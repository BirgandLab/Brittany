#here we make the output for Dani
#have colnames reflect actual contents (ie name of constituent)
#date flow predicted_concentration LabValue rowSumFouling
#load some mountain of data, and set some basic environmental parameters
source("C:/Users/FBlab/Desktop/Brittany/outputForDani.env.R", echo=FALSE)

#file types are cryptic, but the code is explained here
#1  original Types<-c("original","prunedO",
#2  "pruned original" (with bad points removed)        "1stDer","pruned1D",
#3  1st derivative       "turbComp","prunedTC",
#4  "pruned" 1st derivative       "1stDerTurbComp","prunedTC1D") 
#5  turbidity compensated
#6  "pruned" turbidity compensated
#7  first derivative turbidity compensated
#8  "pruned" first derivative turbdity compensated

#chem #is used here too, and is similary cryptic. Here is the key to the code
#not all constituents are available in all years, based on the lab data we had to analyze
#1 DIC  #2 DOC      #3 MES 
#4 NNO3 #5 PPO43    #6PTOT 
#7SO4   #8 Turbdity #9 CL 
#10NO2  #11 Ntot    #12 NTotFilt 
#13NNH4 #14 Silica
#Chem[chemN[2]] changing the integer in this statement will print the name of the chemical

#for 2010-2011 #HERE we are running these with pruned files, but dropping the window data
#ksd means we are making predictions for all of the data, even those in the "bad" windows
#by changing "windows" to "range" in rangeOrwindows, calibration points contained in the fouled windows would
#be included in the prediction. In fact, I am making predictions with "range" specified, hoping to constrain
#the wildness of predictions in that region. Perhaps aquarius would help
mp2010<-as.data.frame(matrix(ncol=8,nrow=13))
colnames(mp2010)<-c("year","chem#","nComp","ksd","fileType","RorW","subset","run")
mp2010[1,]<-c(2010,1,5,1,6,"range",0,1) #DIC
mp2010[2,]<-c(2010,2,5,1,6,"range",0,1) #DOC
mp2010[3,]<-c(2010,3,4,1,6,"range",0,1) #MES
mp2010[4,]<-c(2010,4,4,1,6,"range",0,1) #NNO3
mp2010[5,]<-c(2010,5,8,1,7,"range",0,1) #PPO43
mp2010[6,]<-c(2010,6,5,1,2,"range",0,1) #PTOT
mp2010[7,]<-c(2010,7,5,1,6,"range",0,1) #SO4
mp2010[8,]<-c(2010,8,4,1,8,"range",0,1) #Turb
mp2010[9,]<-c(2010,9,3,1,4,"range",0,1) #CL
mp2010[10,]<-c(2010,10,3,1,2,"range",0,1) #NO2 10
mp2010[11,]<-c(2010,11,3,1,2,"range",0,1) #Ntot 11
mp2010[12,]<-c(2010,12,3,1,6,"range",0,1) #NtotFilt 12
mp2010[13,]<-c(2010,13,4,1,2,"range",0,1) #NNH4 13


mp2011<-as.data.frame(matrix(ncol=8,nrow=9))
colnames(mp2011)<-c("year","chem#","nComp","ksd","fileType","RorW","subset")
mp2011[1,]<-c(2011,1,4,1,6,"range",0,1) #DIC
mp2011[2,]<-c(2011,2,5,1,6,"range",0,1) #DOC
mp2011[3,]<-c(2011,3,8,1,2,"range",0,1) #MES
mp2011[4,]<-c(2011,4,4,1,6,"range",0,1) #NNO3
mp2011[5,]<-c(2011,5,8,1,8,"range",0,1) #PPO43
mp2011[6,]<-c(2011,6,5,1,2,"range",0,0) #PTOT
mp2011[7,]<-c(2011,7,5,1,6,"range",0,1) #SO4
mp2011[8,]<-c(2011,8,4,1,8,"range",0,1) #Turb
mp2011[9,]<-c(2011,9,6,1,6,"range",0,1) #CL


mp2012<-as.data.frame(matrix(ncol=8,nrow=10))
colnames(mp2012)<-c("year","chem#","nComp","ksd","fileType","RorW","subset")
mp2012[1,]<-c(2012,1,6,1,6,"range",0,1) #DIC
mp2012[2,]<-c(2012,2,5,1,6,"range",0,1) #DOC
mp2012[3,]<-c(2012,3,4,1,6,"range",0,1) #MES
mp2012[4,]<-c(2012,4,4,1,6,"range",0,1) #NNO3
mp2012[5,]<-c(2012,5,8,1,6,"range",0,1) #PPO43
mp2012[6,]<-c(2012,6,5,1,2,"range",0,0) #PTOT
mp2012[7,]<-c(2012,7,5,1,6,"range",0,1) #SO4
mp2012[8,]<-c(2012,8,4,1,2,"range",0,1) #Turb
mp2012[9,]<-c(2012,9,4,1,6,"range",0,1) #CL
mp2012[10,]<-c(2012,14,5,1,2,"range",0,1) #Silica
#test<-predictForDani(modelParameters=mp2012 ,flow=Flow,fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut,row=3) #run just one row to check output etc

mpall<-as.data.frame(matrix(ncol=8,nrow=14))
colnames(mpall)<-c("year","chem#","nComp","ksd","fileType","RorW","subset")
mpall[1,]<-c(42,1,8,1,8,"range",0,1) #DIC
mpall[2,]<-c(42,2,5,1,6,"range",0,1) #DOC
mpall[3,]<-c(42,3,3,1,6,"range",0,1) #MES
mpall[4,]<-c(42,4,4,1,6,"range",0,1) #NNO3
mpall[5,]<-c(42,5,8,1,2,"range",0,1) #PPO43
mpall[6,]<-c(42,6,5,1,2,"range",0,0) #PTOT
mpall[7,]<-c(42,7,5,1,6,"range",0,1) #SO4
mpall[8,]<-c(42,8,4,1,2,"range",0,1) #Turb
mpall[9,]<-c(42,9,6,1,6,"range",0,1) #CL
mpall[10,]<-c(42,10,3,1,2,"range",0,1) #NO2 10
mpall[11,]<-c(42,11,3,1,2,"range",0,1) #Ntot 11
mpall[12,]<-c(42,12,3,1,6,"range",0,1) #NtotFilt 12
mpall[13,]<-c(42,13,4,1,2,"range",0,1) #NNH4 13
mpall[14,]<-c(42,14,5,1,2,"range",0,1) #Silica



#test<-predictForDani(modelParameters=mp2012[1:3,] ,flow=Flow,fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)
result2010<-predictForDani(modelParameters=mp2010,flow=Flow,fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)
            result<-as.data.frame(result2010$Predictions)
            result$realTime<-as.POSIXlt(result$realTime,origin="1970-1-1",tz="UTC")
            write.table(result, file = paste(OutputPath,"2010_Predictions.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")
            
            write.table(result2010$modelStats, file = paste(OutputPath,"2010_modelQuality.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")

result2011<-predictForDani(modelParameters=mp2011,flow=Flow,fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)
            result<-as.data.frame(result2011$Predictions)
            result$realTime<-as.POSIXlt(result$realTime,origin="1970-1-1",tz="UTC")
            write.table(result, file = paste(OutputPath,"2011_Predictions.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")
            
            write.table(result2011$modelStats, file = paste(OutputPath,"2011_modelQuality.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")




result2012<-predictForDani(modelParameters=mp2012,flow=Flow,fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)
            result<-as.data.frame(result2012$Predictions)
            result$realTime<-as.POSIXlt(result$realTime,origin="1970-1-1",tz="UTC")
            write.table(result, file = paste(OutputPath,"2012_Predictions.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")
            
            write.table(result2012$modelStats, file = paste(OutputPath,"2012_modelQuality.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")

result2010_2012<-predictForDani(modelParameters=mpall,flow=Flow,fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)
            result<-as.data.frame(result2010_2012$Predictions)
            result$realTime<-as.POSIXlt(result$realTime,origin="1970-1-1",tz="UTC")
            write.table(result, file = paste(OutputPath,"2010-2013_Predictions.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")
            
            write.table(result2010_2012$modelStats, file = paste(OutputPath,"2010-2013_modelQuality.csv"),
                        row.names=FALSE,col.names = TRUE, qmethod = c("escape", "double"),
                        fileEncoding = "",sep=",")

