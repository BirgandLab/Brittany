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
				labDataCols=labDataCols,analytenames=NULL,timezone="UTC"){

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


