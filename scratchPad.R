#2010-2011 points


#2010 breaks
startDates<-c(as.POSIXct(paste("2010-11-18 00:00:00",sep=""),tz="UTC"),
                  as.POSIXct(paste("2011-04-25 00:00:00",sep=""),tz="UTC")
)
stopDates<-c(as.POSIXct(paste("2011-03-01 00:00:00",sep=""),tz="UTC"),
                 as.POSIXct(paste("2011-05-20 00:00:00",sep=""),tz="UTC")
                 )


#2011 breaks
startDates<-c(as.POSIXct(paste("2011-10-18 00:00:00",sep=""),tz="UTC"))
stopDates<-c(as.POSIXct(paste("2011-07-01 00:00:00",sep=""),tz="UTC"))










trimCalibOrig<-subset(originalmyData,"calibration",trimStartDates,trimStopDates,chem=NULL)
trimCalib1Der<-subset(DermyData,"calibration",trimStartDates,trimStopDates,chem=NULL)
trimCalibTC<-subset(TCmyData,"calibration",trimStartDates,trimStopDates,chem=NULL)
trimCalib1DerTC<-subset(TC1DmyData,"calibration",trimStartDates,trimStopDates,chem=NULL)




good3Start<-
good3Stop<-as.POSIXct(paste("2012-04-24 00:00:00",sep=""),tz="UTC")

good4Start<-as.POSIXct(paste("2012-05-05 00:00:00",sep=""),tz="UTC")
good4Stop<-as.POSIXct(paste("2012-06-01 00:00:00",sep=""),tz="UTC")

good5Start<-as.POSIXct(paste("2012-06-28 00:00:00",sep=""),tz="UTC")
good5Stop<-as.POSIXct(paste("2012-07-10 00:00:00",sep=""),tz="UTC")

good6Start<-as.POSIXct(paste("2012-10-20 00:00:00",sep=""),tz="UTC")
good6Stop<-as.POSIXct(paste("2013-03-27 00:00:00",sep=""),tz="UTC")

good7Start<-as.POSIXct(paste("2013-10-20 00:00:00",sep=""),tz="UTC")
good7Stop<-as.POSIXct(paste("2014-04-05 00:00:00",sep=""),tz="UTC")





































#lets have a generic function for loading data
#it needs a filepath--to calibration data
#it needs a filename
#it needs a fileType
#it needs a list of columns with chem Data
#it needs a list of columns with NaN data (wavelengths that don't read)
#well, it doesn't really need the last part, as we can find those columns
#it might benefit from specification of a timezone
#OUTPUTS
#return a matrix of 
#	realTime 		-- as POSIX object with timezone specified
#	labData 		-- IF LAB DATA COLS ARE SPECIFIED
#	fingerPrints	-- with columns containing data
#     start and stop time in a vector
#	summary stats for each lab analyte  --IF LAB DATA COLS ARE SPECIFIED
#	an indicator of the type of file it was supposed to have been ie fingerprint, firstderivative, turbiditycompensated, or turbiditycompensatedfirstderivative
#
#can we use the same thing for loading fingerprint data without chemData?
#labDataCols==NULL
loadGeneric<-function(filepath=filepath,filename=filename,fileType=fileType,
				labDataCols=labDataCols,analytenames=NULL,timezone="UTC")
  {

   data<-read.table(file=paste(filepath,filename,sep=""),sep=",",header=TRUE,skip=0)
   
   realTime<-strptime(data[,1],'%d/%m/%Y %H:%M:%S',tz=timezone)
   #if there is no time for an observation, we will disregard it
   		hasTime<-complete.cases(realTime)
   		realTime<-realTime[hasTime]
		data<-data[hasTime,]
   
if(labDataCols[1]!=0){
		   labData<-data[,labDataCols]
		  #drop lab data from matrix
  		   data<-data[,-min(labDataCols):-dim(data)[2]]
 		 #calculate some stas on lab data
  		     labDataStats<-rbind(colMeans(labData,na.rm=TRUE),
			    apply(labData,2,function(x) min(x,na.rm=TRUE)),
			    apply(labData,2,function(x) max(x,na.rm=TRUE)),
			    apply(labData,2,function(x) sum(!is.na(x)))
				)  
		  	rownames(labDataStats)<-c("mean","min","max","count") 
	}

  #drop realTime, and OK columns from beginning of matrix
  	data<-data[,-1:-2]
  #find columns without any fingerprint data in them
      nanCount<-apply(data,2,function(x) sum(is.na(x))) #calculate the total number of NaNs in a column
	numberOfRows<-dim(data)[1]  			    #figure out how many empty columns there are
	NanColumns<-nanCount>=(numberOfRows*0.99)         #here is a logical vector of ~empty columns
  #drop them 
      fingerPrints<-data[,-NanColumns]    

  dateRange<-c(min(realTime),max(realTime))
  
if(labDataCols[1]!=0){
		return(list(realTime=realTime, 
			labData=labData, 
			fingerPrints=fingerPrints,
			stats=labDataStats,
            	dateRange=dateRange,
			fileType=fileType,
			dataType="calibration"
			))
		}
return(list(realTime=realTime, 
			fingerPrints=fingerPrints,
            	dateRange=dateRange,
			fileType=fileType,
			dataType="fingerprint"
			))
}


















#smaller window
#plot one week at a time
#make vertical line for noon and midnight

magWStart<-as.POSIXct("2011-4-01 00:00:00",tz="UTC") #jan 1 2010
magWStop<-as.POSIXct("2011-4-15 00:00:00",tz="UTC")    #june 15 2011
hex<-(AP_R$x>as.numeric(magWStart) & AP_R$x<as.numeric(magWStop))



plot(as.POSIXct(AP_R$x[hex],origin="1970-01-01 00:00:00"),
     AP_R$y[hex],type="l",col="gray",
     main=chemicals[chemical],
 #    ylim=c(min(AP_R$y),max( AP_R$y)),
#     xlim=c(min( AP_R$x),max( AP_R$x))
)
butt <- butter(2,0.01,type="low")
points(as.POSIXct(AP_R$x[hex],origin="1970-01-01 00:00:00"),
       filtfilt(butt,AP_R$y[hex]),type="l",col="red"
       )
abline(h=0)
for(i in magWStart:magWStop){
  if(as.numeric(i)%%86400==0) {
    abline(v=i+60*60*4) 
  }
  
}




startDate<-as.POSIXct(paste("2010-10-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
#startDate<-as.POSIXct(paste(startYear[1],"-11-01 00:00:00",sep=""),tz="UTC") #jan 1 2010
stopDate<-as.POSIXct(paste("2011-07-30 00:00:00",sep=""),tz="UTC")    #june 15 2011

goodStart<-as.POSIXct(paste("2010-11-18 00:00:00",sep=""),tz="UTC")
goodStop<-as.POSIXct(paste("2011-03-15 00:00:00",sep=""),tz="UTC")
good2Start<-as.POSIXct(paste("2011-04-17 00:00:00",sep=""),tz="UTC")
good2Stop<-as.POSIXct(paste("2011-06-15 00:00:00",sep=""),tz="UTC")

good3Start<-as.POSIXct(paste("2011-12-24 00:00:00",sep=""),tz="UTC")
good3Stop<-as.POSIXct(paste("2012-04-24 00:00:00",sep=""),tz="UTC")

good4Start<-as.POSIXct(paste("2012-05-05 00:00:00",sep=""),tz="UTC")
good4Stop<-as.POSIXct(paste("2012-06-01 00:00:00",sep=""),tz="UTC")

good5Start<-as.POSIXct(paste("2012-06-28 00:00:00",sep=""),tz="UTC")
good5Stop<-as.POSIXct(paste("2012-07-10 00:00:00",sep=""),tz="UTC")

good6Start<-as.POSIXct(paste("2012-10-20 00:00:00",sep=""),tz="UTC")
good6Stop<-as.POSIXct(paste("2013-03-27 00:00:00",sep=""),tz="UTC")

good7Start<-as.POSIXct(paste("2013-10-20 00:00:00",sep=""),tz="UTC")
good7Stop<-as.POSIXct(paste("2014-04-05 00:00:00",sep=""),tz="UTC")

foo<-realTime>startDate&realTime<stopDate
bar<-(realTime>goodStart&realTime<goodStop | 
        realTime>good2Start&realTime<good2Stop | 
        realTime>good3Start&realTime<good3Stop |
        realTime>good4Start&realTime<good4Stop  |
        realTime>good5Start&realTime<good5Stop  |
        realTime>good6Start&realTime<good6Stop  |
        realTime>good7Start&realTime<good7Stop)
w<-110
plot(realTime[foo],fingerPrints[foo,w],type="l",col="red")
points(realTime[bar],fingerPrints[bar,w],type="l",col="black")
man<-as.POSIXct(AP_R$x,origin="1970-01-01 00:00:00")
blah<-(man>goodStart&man<goodStop | 
        man>good2Start&man<good2Stop)











#looking at what makes a turbidity compensated file
realTime
spec<-2000

plot(seq(200,740,2.5)[2:210],as.numeric(FingerPrints[spec,2:210]),type="l",col="black",ylim=c(-10,600))
points(seq(200,740,2.5)[2:210],as.numeric(TCFingerPrints[spec,2:210]),type="l",col="red")
plot(seq(200,740,2.5)[2:210],as.numeric(TC1DFingerPrints[spec,2:210]),type="l",col="red")
points(seq(200,740,2.5)[2:210],as.numeric(DerFingerPrints[spec,2:210]),type="l",col="black")

plot(seq(200,740,2.5)[2:210],(as.numeric(TCFingerPrints[spec,2:210]))-as.numeric(fingerPrints[spec,2:210]),type="l",col="green")



bar<-read.table(file="C:/Users/FBlab/Desktop/work_here/Data/Brittany/flow/Débit Kervidy Naizin 2010-2012.txt",
                skip=1,sep=";",dec=",",na.strings="---")
bar<-read.table(file="C:/Users/FBlab/Desktop/work_here/Data/Brittany/flow/Débit Kervidy Naizin 2012-2014.txt",
        sep=";",skip=1, dec=".",na.strings="---")
#these 2012 through 2014 flow data are not corrected. great big mess. :(
Date<-bar$V1        #the date
T<-bar$V2     #the time
D<-paste(Date,T,sep=" ")               #put back in the fixed values
D<-strptime(D, '%d/%m/%Y %H:%M:%S',tz="UTC")    #convert it to a _flawless_ time opject
bar$realTime<-D    
bar$flow<-bar$V3
flow<-flow[,-1:-3]
flow<-flow[complete.cases(flow$flow),]
flow<-flow[complete.cases(flow$realTime),]# removes all the rows for which there is a NA--keeping time in there
remove(D,T,Date)




fitFile<-paste(fitPath,  
               Chem[chemN[chemical]],
               "_",
               "test",
               "_",
               numComp,
               "_",
               (as.POSIXlt(startDate))$year+1900,
               "PLSR.in",sep=""
)
fitFileOut<-paste(fitPath, 
                  Chem[chemN[chemical]],
                  "_",
                  "test",
                  "_",
                  numComp,
                  "_",
                  (as.POSIXlt(startDate))$year+1900,
                  "PLSR_out.txt",sep=""
)
modelOutput<-modelExecution(chemical,numComp,subsetRate,calibration,specDataToModel,fitEval,fitFile,fitFileOut)

excludeRows<-c(29,115,153,155,200,250,471,482,486,514,525,530,548,615,616,617,652,660,668,677,671,673,735,762,790,817,831,840,841,845,902,911,934,945)
prunedO<-loadDataFile(CalibrationFingerPrintsPath,ModelFilename[1],excludeRows) 
foo<-subsetSpecData("prunedO","calibration",startDates,stopDates,chemN[1])

plot(foo$realTime,foo$ChemData)











# begin<-min(dateWindows)
# end<-max(dateWindows) 
# totmin<-as.numeric(end-begin)/60
# 
# #calculate the flow values for all 1 minute intervals within the year of interest
# flowApprox<-approx(flow$realTime,flow$flow,n=totmin+1) #create a time series of flow at 1 minute 
# flowApprox$x<-align.time(flowApprox$x,60)
# a<-merge(modelData,flowApprox,by.x="realTime",by.y="x",all.x="TRUE")
# b<-merge(calibData,flowApprox,by.x="realTime",by.y="x",all.x="TRUE")










###run a model
a<-runModel(chem="DIC",
         numComp=5,
         dates=cbind(startDates,stopDates),
         calibData=,
         modelData=TurbidityCompensated,
         type="turbComp",
         flow=Flow,
         subset=0.75,
         iterations=100,
         fitEval=fitEval,fitFile=fitFile,fitFileOut=fitFileOut)










runIterativeModel(numComp=5,iterations=1,subsetRatio=0,
 chemToRun=c("DIC"),
 calibData,
 specData,
 fitFilePaths,
 timeWindow)



#run a set of models with specified parameters
runIterativeModel<-function(numComp=5,iterations=1,subsetRatio=0,
                  chemToRun=NULL,calibData=NULL,specData=NULL,
                  fitFilePaths=NULL,timeWindow=NULL)
{
  #ssDATA<-subsetAll(calibData,specData,chem=chemToRun,flow=Flow, dateWindows=timeWindow)
  a<-runModel(chem=chemToRun,
                  numComp=numComp,
                  dates=timeWindow,
                  calibData=calibData,
                  modelData=specData,
                  type=" ",
                  flow=flow,
                  subset=subsetRatio,
                  iterations=iterations,
                  fitEval=filePaths[1],fitFile=filePaths[2],fitFileOut=filePaths[3])
      
      #splat out files that might be reloaded for analytical purposes laytah
      #this is a list with several types of data
      fN[i]<-paste(chem,subsetRatio,"_",iterations,"_",as.POSIXlt(min(timeWindow),origin="1970-1-1")$year+1900,"-",as.POSIXlt(max(timeWindow),origin="1970-1-1")$year+1900,sep="")
      save(a,file=fN[i]);

      # close(fN[i])
}

###SHEER BRILLIANCE: combine all types of fingerprints and see if PLSR makes better fits this way!!!!
#monster Calib Data



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






#mcd<-cbind(originalmyData$fingerPrints,TCmyData$fingerPrints,
#           DermyData$fingerPrints,TC1DmyData$fingerPrints)
#MCD<-originalmyData
#mcd<-cbind(prunedO$fingerPrints,prunedTC$fingerPrints,
#           pruned1D$fingerPrints,prunedTC1D$fingerPrints)

#MCD<-pruned1D


MCD$fingerPrints<-mcd
rm(mcd)
mcd<-cbind(original$fingerPrints,TurbidityCompensated$fingerPrints,
           Derivative$fingerPrints,FstDerivativeTurbidityCompensated$fingerPrints)
MSD<-original
MSD$fingerPrints<-mcd
rm(mcd)







startDate<-min(startDates)
stopDate<-max(stopDates)
flow<-subsetAll(flow=Flow,dateWindows=cbind(startDates,stopDates))$flow
for(year in 1:4){
for(chemical in 1:14){
  for(numComp in 3:9){
#chemical<-1
#numComp<-6
subsetRate<-0
all<-subsetAll(calibData=MCD,modelData=MSD,
               chem=chemN[chemical],
               dateWindows=cbind(startDates,stopDates),
               keepSpecData=0)
fileType<-"monsterup"
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
                            fitEval,fitFile,fitFileOut
                            )

png(file = paste(OutputPath,
                 Chem[chemN[chemical]],
                 "_",
                 fileType,
                 "_",
                 numComp,
                 "_",
                 (as.POSIXlt(startDate))$year+1900,
                 ".png",sep=""),
    unit="in",
    width=8,height=4,res=300
    
)
plot(as.POSIXct(flow$realTime,origin="1970-01-01 00:00:00"),
     flow$flow,col="black",type="l",xaxt="n",yaxt="n",
     xlab="date",ylab="",
     xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE))
    )

par(new=TRUE)
plotMax<-max(modelOutput$ObservedandPredicted[,1])*1.5
plotMin<-min((min(modelOutput$ObservedandPredicted[,1])-min(modelOutput$ObservedandPredicted[,1])*0.1),0)


plot(as.POSIXct(modelOutput$ObservedandPredicted[,3],origin="1970-01-01 00:00:00"),
     modelOutput$ObservedandPredicted[,1],
     xlab="",
     type="p",
     col="green",
     ylim=c(plotMin,plotMax),
     ylab=paste(chemicals[chemical],sep=" "),
     cex.lab=1.3,
     pch=20,
     cex=0.7,
    xlim=c(min(modelOutput$PredictedConcentrations[,1]),max(modelOutput$PredictedConcentrations[,1],na.rm=TRUE))
     )

points(modelOutput$ObservedandPredicted[,3],
       modelOutput$ObservedandPredicted[,2],
       type="p",
       pch=19,
       cex=0.7,
       col="blue"
       )

points(modelOutput$PredictedConcentrations[,1],
       modelOutput$PredictedConcentrations[,2],
       type="l",
       col="red"
       )

legend("topright",
       c("hydrograph","lab value","pred lab","pred ts"),
       pch=c(NA,20,19,NA),
       lty=c(1,NA,NA,1),
       col=c("black","green","blue","red"),cex=0.8
       )
dev.off()
}#for the #comp loop
}#for each chem
}#for each time window