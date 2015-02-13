
#load and condition flow data

flow<-read.table(file="C:/Users/FBlab/Desktop/work_here/Data/Brittany/flow/Débit Kervidy Naizin 2010-2012.txt",
               skip=1,sep=";",dec=",",na.strings="---")
#flow<-read.table(file="C:/Users/FBlab/Desktop/work_here/Data/Brittany/flow/Débit Kervidy Naizin 2010-2012.txt",skip=1,sep=";")
#flow<-read.table(file="C:/Users/FBlab/Desktop/work_here/Data/Brittany/flow/Débit Kervidy Naizin 2012-2014.txt",sep=";",skip=1)
#these 2012 through 2014 flow data are not corrected. great big mess. :(
Date<-flow$V1        #the date
T<-flow$V2     #the time
D<-paste(Date,T,sep=" ")               #put back in the fixed values
D<-strptime(D, '%d/%m/%Y %H:%M:%S',tz="UTC")    #convert it to a _flawless_ time opject
flow$realTime<-D    
flow$flow<-flow$V3
flow<-flow[,-1:-3]
flow<-flow[complete.cases(flow$flow),]
flow<-flow[complete.cases(flow$realTime),]# removes all the rows for which there is a NA--keeping time in there
remove(D,T,Date)


filename<-matrix(nrow=4,ncol=4)
#LOAD FINGERPRINTS FROM DEPLOYMENT FOR MODELING CONCENTRATIONS
filename[1,]<-c("S2010-2011.fp" ,"S2011-2012.fp","S2012-2013.fp","S2013-2014.fp")
filename[2,]<-c("S2010-20111stDer.fp", "S2011-20121stDer.fp","S2012-20131stDer.fp","S2013-20141stDer.fp")
filename[3,]<-c("S2010-2011TurbComp.fp" ,"S2011-2012TurbComp.fp","S2012-2013TurbComp.fp","S2013-2014TurbComp.fp")
filename[4,]<-c("2010-20111stDerTurbComp.fp", "2011-20121stDerTurbComp.fp","2012-20131stDerTurbComp.fp","2013-20141stDerTurbComp.fp")


#LOAD DATA FROM BRITTANY FOR PROJECTION
#load fp
OrigmyData<-loadDataFile(path,ModelFilename[1])   #returns an object with    $realTime $ChemData and $fingerPrints 

data0<-read.table(file=paste(FPpath,filename[1,1],sep=""),sep="\t",header=TRUE,skip=1)
data1<-read.table(file=paste(FPpath,filename[1,2],sep=""),sep="\t",header=TRUE,skip=1)
data2<-read.table(file=paste(FPpath,filename[1,3],sep=""),sep="\t",header=TRUE,skip=1)
data3<-read.table(file=paste(FPpath,filename[1,4],sep=""),sep="\t",header=TRUE,skip=1)
data<-rbind(data0,data1,data2,data3)
remove(data0,data1,data2,data3)
data<-data[data$Status_0=="Ok",]

realTime<-strptime(data$Date.Time, "%Y.%m.%d  %H:%M:%S",tz="UTC")                       #store it back in the matrix
fingerPrints<-data[,-(dim(data)[2]-3):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
fingerPrints<-fingerPrints[,-1:-2]  
remove(data)
#align data to closest minute
realTime<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)
#left with fingerPrints and realTime
FingerPrints<-fingerPrints
RealTime<-realTime




#load 1stDerivative)

DermyData<-loadDataFile(path,ModelFilename[2])   #returns an object with    $realTime $ChemData and $fingerPrints 

data0<-read.table(file=paste(FPpath,filename[2,1],sep=""),sep="\t",header=TRUE,skip=1)
data1<-read.table(file=paste(FPpath,filename[2,2],sep=""),sep="\t",header=TRUE,skip=1)
data2<-read.table(file=paste(FPpath,filename[2,3],sep=""),sep="\t",header=TRUE,skip=1)
data3<-read.table(file=paste(FPpath,filename[2,4],sep=""),sep="\t",header=TRUE,skip=1)
data<-rbind(data0,data1,data2,data3)
remove(data0,data1,data2,data3)
data<-data[data$Status_2=="Ok",]

realTime<-strptime(data$Date.Time, "%Y.%m.%d  %H:%M:%S",tz="UTC")                       #store it back in the matrix
fingerPrints<-data[,-(dim(data)[2]-3):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
fingerPrints<-fingerPrints[,-1:-2]  
remove(data)
#align data to closest minute
realTime<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)
#left with fingerPrints and realTime
DerFingerPrints<-fingerPrints
DerRealTime<-realTime



#load turb comp
TCmyData<-loadDataFile(path,ModelFilename[3])   #returns an object with    $realTime $ChemData and $fingerPrints 

data0<-read.table(file=paste(FPpath,filename[3,1],sep=""),sep="\t",header=TRUE,skip=1)
data1<-read.table(file=paste(FPpath,filename[3,2],sep=""),sep="\t",header=TRUE,skip=1)
data2<-read.table(file=paste(FPpath,filename[3,3],sep=""),sep="\t",header=TRUE,skip=1)
data3<-read.table(file=paste(FPpath,filename[3,4],sep=""),sep="\t",header=TRUE,skip=1)
data<-rbind(data0,data1,data2,data3)
remove(data0,data1,data2,data3)
data<-data[data$Status_1=="Ok",]


realTime<-strptime(data$Date.Time, "%Y.%m.%d  %H:%M:%S",tz="UTC")                       #store it back in the matrix
fingerPrints<-data[,-(dim(data)[2]-3):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
fingerPrints<-fingerPrints[,-1:-2]  
remove(data)
#align data to closest minute
realTime<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)
#left with fingerPrints and realTime
TCFingerPrints<-fingerPrints
TCRealTime<-realTime

#load 1stDer Turb Comp
TC1DmyData<-loadDataFile(path,ModelFilename[4])   #returns an object with    $realTime $ChemData and $fingerPrints 

data0<-read.table(file=paste(FPpath,filename[4,1],sep=""),sep="\t",header=TRUE,skip=1)
data1<-read.table(file=paste(FPpath,filename[4,2],sep=""),sep="\t",header=TRUE,skip=1)
data2<-read.table(file=paste(FPpath,filename[4,3],sep=""),sep="\t",header=TRUE,skip=1)
data3<-read.table(file=paste(FPpath,filename[4,4],sep=""),sep="\t",header=TRUE,skip=1)
data<-rbind(data0,data1,data2,data3)
remove(data0,data1,data2,data3)
data<-data[data$Status_5=="Ok",]

realTime<-strptime(data$Date.Time, "%Y.%m.%d  %H:%M:%S",tz="UTC")                       #store it back in the matrix
fingerPrints<-data[,-(dim(data)[2]-3):-(dim(data)[2])]    #remove chemical data and the 4 columns of NaNs in 742.5-750 nM bins
fingerPrints<-fingerPrints[,-1:-2]  
remove(data)
#align data to closest minute
realTime<-align.time(realTime,60) #rounds UP to begining of next interval (specified in seconds)
#left with fingerPrints and realTime
TC1DFingerPrints<-fingerPrints
TC1DRealTime<-realTime