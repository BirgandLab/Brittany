
# #IDENTIFTY TIME WINDOW FOR THE "ANNUAL MODEL"
# startDate<-c(as.POSIXct("2010-10-01 00:00:00",tz="UTC"),
#              as.POSIXct("2012-01-01 00:00:00",tz="UTC"),
#              as.POSIXct("2013-01-01 00:00:00",tz="UTC")
#              )
#              
# stopDate<-c(as.POSIXct("2011-01-1 00:00:00",tz="UTC"),
#             as.POSIXct("2012-03-1 00:00:00",tz="UTC"),            
#             as.POSIXct("2013-03-1 00:00:00",tz="UTC")
#           )
# a<-subset(originalmyData,"calibration",startDate,stopDate)
# b<-subset(original,"fingerPrints",startDate,stopDate)
# c<-subset(flow,"flow",startDate,stopDate)


subset<-function(data,type,startDates,stopDates){
  
  logicalIndexofInclusion<-matrix(nrow=length(data$realTime),ncol=length(startDates))
  
  for(i in 1:length(startDates)){
    logicalIndexofInclusion[,i]<-(data$realTime>startDates[i]&data$realTime<stopDates[i])
  }
    
  keep<-as.logical(rowSums(logicalIndexofInclusion,na.rm=TRUE))
  

  if (type=="calibration"){
    realTime<-data$realTime[keep]
    fingerprint<-data$fingerPrint[keep,]
    ChemData<-data$ChemData[keep,]
    return(list(realTime=realTime,fingerPrints=fingerprint,ChemData=ChemData))
    #returnStatement
  }
  
  if (type=="fingerPrints"){
    realTime<-data$realTime[keep]
    fingerprint<-data$fingerPrint[keep,]
    return(list(realTime=realTime,fingerPrints=fingerprint))
    
  }
  
  if(type=="flow"){
    realTime<-data$realTime[keep]
    flow<-data$flow[keep]
    return(list(realTime=realTime,flow=flow))
  }
    
}