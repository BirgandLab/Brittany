#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Pre-Load Data for Flow *************************************************************#
#*******************************************************************************************************************/  
#*******************************************************************************************************************/
#load and condition flow data
message("load flow data")
      flow1<-read.table(file=paste(flowPath,"Débit Kervidy Naizin 2010-2012.txt",sep=""),
                     skip=1,sep=";",dec=",",na.strings="---")
      
      flow2<-read.table(file=paste(flowPath,"Débit Kervidy Naizin 2012-2014.txt",sep=""),
                    skip=1, sep=";",na.strings="---")

      #put the two matrices together
      flow<-rbind(flow1,flow2)
        #clean up
        remove(flow1,flow2)


      #calculate a time object from the date and time
      flow$realTime<-strptime(paste(flow$V1,flow$V2,sep=" "),'%d/%m/%Y %H:%M:%S',tz="UTC")               #put back in the fixed values

      #rename columns and remove the unnecessary ones
      colnames(flow)[3]<-"flow"
      flow<-flow[,-1:-2]
  
      #remove any missing data
      Flow<-flow[complete.cases(flow$flow),]
      Flow<-flow[complete.cases(flow$realTime),]# removes all the rows for which there is a NA--keeping time in there

      #plot(flow$realTime,flow$flow,type="l")

#*******************************************************************************************************************/
#*******************************************************************************************************************/
  


#*******************************************************************************************************************/
#*******************************************************************************************************************/
#*******************************Pre-Load FingerPrint Data for Predictions*******************************************#
#*******************************************************************************************************************/  
#*******************************************************************************************************************/

#definitions of filenames
      filename<-matrix(nrow=4,ncol=4)
      #the unmodified fingerprints
        filename[1,]<-c("S2010-2011.fp" ,"S2011-2012.fp","S2012-2013.fp","S2013-2014.fp")
      #first derivative fingerprints
        filename[2,]<-c("S2010-20111stDer.fp", "S2011-20121stDer.fp","S2012-20131stDer.fp","S2013-20141stDer.fp")
      #Turbidity Compensated Finger Prints
        filename[3,]<-c("S2010-2011TurbComp.fp" ,"S2011-2012TurbComp.fp","S2012-2013TurbComp.fp","S2013-2014TurbComp.fp")
      #1st Derivative Turbidity Compensated Finger Prints
        filename[4,]<-c("2010-20111stDerTurbComp.fp", "2011-20121stDerTurbComp.fp","2012-20131stDerTurbComp.fp","2013-20141stDerTurbComp.fp")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Load the unmodified finger print records %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#load the calibration fingerprints with their associated lab data
       message("load calibration data")
        originalmyData<-loadDataFile(CalibrationFingerPrintsPath,ModelFilename[1])   
                    #returns a list with objects  $realTime $ChemData and $fingerPrints 
#load the fingerprint data file for each year
       message("load fingerprints for prediction")
        original<-loadFingerPrints(FingerPrintPath,filename,"original")
                    #returns a list with $realTime and $fingerprints

      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#load 1stDerivative calibration data and observation fingerprints
message("load calibration data first derivative")
        DermyData<-loadDataFile(CalibrationFingerPrintsPath,ModelFilename[2])
                    #returns an object with    $realTime $ChemData and $fingerPrints
message("load fingerprint for predictions")
        Derivative<-loadFingerPrints(FingerPrintPath,filename,"1stDerivative")
                    #returns a list with $realTime and $fingerprints

      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#load turbidity compensated calibration data and fingerprints
message("load calibration data turbidity compensated")
        TCmyData<-loadDataFile(CalibrationFingerPrintsPath,ModelFilename[3])   
                    #returns an object with $realTime $ChemData and $fingerPrints 
message("load fingerprints for prediction")
        TurbidityCompensated<-loadFingerPrints(FingerPrintPath,filename,"TurbidityCompensated")

      
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#load 1stDer Turb Comp
message("load calibration data turbidity compensated first derivative")
        TC1DmyData<-loadDataFile(CalibrationFingerPrintsPath,ModelFilename[4])   
                    #returns an object with    $realTime $ChemData and $fingerPrints 
message("load fingerprints for prediction")
        FstDerivativeTurbidityCompensated<-loadFingerPrints(FingerPrintPath,filename,"1stDerivativeTurbidityCompensated")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
