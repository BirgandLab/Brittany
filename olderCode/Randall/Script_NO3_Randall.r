pathD<-"C:/Users/birgand/Documents/Publications/2014/PLSR_Brittany/" #Specify folder where calibration/lab data is located
filename<-"Fingerprint_R_1D.csv"  #Specify file where calibration/lab data is located
data<-read.table(file=paste(pathD,filename,sep=""),sep=",",skip=1)  #Import calibration/lab data
#data<-data1[-nrow(data1),]
attach(data)

colnames(data) <- c( 'Date/Time', 'Status') #Add column names
colnames(data)[3:(length(data)-2)]<-c(seq(200,750,2.5)) #Add column names

FP<-data[,-1:-10] #Remove NAs at low wavelengths
FP<-FP[,-(dim(FP)[2]-13):-dim(FP)[2]] #Remove NO3-N and NAs at high wavelengths
FP<-data.matrix(FP) #Convert to data matrix
NO3<-data.matrix(data[dim(data)[2]-1])  #Make matrix of NO3-N values
Sal<-data.matrix(data[dim(data)[2]])  #Make matrix of salinity values if included in spreadsheet

library(pls) #Load the pls package

fit<-plsr(NO3~data.matrix(FP),ncomp=14)  #Create PLSR with the number of components from cross validation
summary(fit)  #Summary of model

setwd("C:/Users/birgand/Documents/Publications/2014/PLSR_Brittany/") #Set the working directory where the output will be saved
pathD<-"C:/Users/birgand/Documents/Publications/2014/PLSR_Brittany/"  #Specify the folder where spectrometer output is located
filename<-"Sandy.fp"  #Specify the file where the spectrometer output is located
data1<-read.table(file=paste(pathD,filename,sep=""),sep="",skip=2) #Import the spectrometer output

colnames(data1) <- c( 'Date', 'Time', 'Status') #Add column names
colnames(data1)[4:length(data1)]<-c(seq(200,750,2.5)) #Add column names

Date_<-substr(data1$Date, 1, 10)  #Extract date from spectrometer data
T<-substr(data1$Time, 1, 8) #Extract time from spectrometer data
D<-paste(Date_,T,sep=" ") #Combine the date and time values
Dat<-strptime(D, "%Y.%m.%d %H:%M:%S") #Create record of date and time

FP5<-data1[,-1:-11] #Remove NAs at low wavelengths
FP5<-FP5[,-(dim(FP5)[2]-11):-dim(FP5)[2]] #Remove NO3-N and NAs at high wavelengths
FP5<-data.matrix(FP5) #Convert spectrometer output to matrix

NO3P<-as.data.frame(matrix(0,dim(data1)[1],2))  #Create data frame for date/time and predicted NO3-N values
NO3P[,1]<-as.character(Dat, "%m/%d/%Y %H:%M:%S")  #Add date/time to data frame

Pfit<-predict(fit,FP5,ncomp=14,type=c("response"))  #Predict NO3-N values using PLSR model and the spectrometer output
NO3P[,2]<-as.data.frame(Pfit[1:dim(data1)[1]])  #Add predicted NO3-N values to data frame

#NO3P[is.na(NO3P$V2),"V2"]=0
NO3P[NO3P$V2<0.1,"V2"]=0.05 #If any of the NO3-N values are lower than the detection limit change to one-half of the detection limit.

write.csv(NO3P,file="NO3_DN.csv",row.names=FALSE) #Export results to .csv file

plot(Dat,NO3P[,2],type="l")  #Plot chemograph
