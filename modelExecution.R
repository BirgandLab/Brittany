

modelExecution<-function(chemical,numComp,subsampleRate,calibData,dataToModel,fitEval,fitFile,fitFileOut){
  modelQuality<-matrix(ncol=26)
  colnames(modelQuality)<-c("M_n","M_NSE",
                            "M_VG","M_G","M_A","M_B",
                            "V_VG","V_G","V_A","V_B",
                            "C_n","C_r2","C_RMSE","C_NRMSE","C_slope",
                            "V_n","V_r2","V_RMSE","V_NRMSE","V_slope",
                            "F_n","F_r2","F_RMSE","F_NRMSE","F_slope",
                            "subR")
  
  #for any model--first subset the data, and pass subsetted data to the execution function
  if(length(calibData$ChemData)>numComp*1.2){
  ModelConcentration<-PLSRFitAndTestNoNSE(calibData$fingerPrints,
                                          calibData$ChemData,
                                          calibData$realTime,
                                          numComp,fitEval,fitFile,fitFileOut,subsampleRate)                  
  #USE THE MODEL TO PREDICT OUTPUT FOR THE TARGET TIME WINDOW        
  PredictedConcentration<-predict(ModelConcentration$Fit,
                                  as.matrix(dataToModel$fingerPrints),
                                  ncomp=numComp,
                                  type=c("response"))  
  PredictedConcentration<-cbind(as.numeric(dataToModel$realTime),PredictedConcentration)

  #calculate stats for goodness of fit etc

  
  modelQuality[1]<-ModelConcentration$Stats[1]
  if(subsampleRate>0){
    modelQuality[2]<-NSE(sim=ModelConcentration$OaP1[,2],obs=ModelConcentration$OaP1[,1])
    modelQuality[3:6]<-OB(ModelConcentration$OaP1[,1],ModelConcentration$OaP1[,2],fitEval,fitFile,fitFileOut)
    modelQuality[7:10]<-OB(ModelConcentration$OaP2[,1],ModelConcentration$OaP2[,2],fitEval,fitFile,fitFileOut)
  }
  if(subsampleRate==0){
    modelQuality[2]<-NSE(sim=ModelConcentration$ObservedAndPredicted[,2],obs=ModelConcentration$ObservedAndPredicted[,1])
    modelQuality[3:6]<-OB(ModelConcentration$ObservedAndPredicted[,1],ModelConcentration$ObservedAndPredicted[,2],fitEval,fitFile,fitFileOut)
  }
  modelQuality[11:25]<-ModelConcentration$Stats
  modelQuality[26]<-subsampleRate
  #return summary statistics, and chemographic prediction vector (in case someone wants to keep it)
  #return(PredictedConcentration)
  #this can be repeated 1x or
  return(list(modelQuality=modelQuality,PredictedConcentrations=PredictedConcentration,ObservedandPredicted=ModelConcentration$ObservedAndPredicted))
  
}
else return(list(modelQuality=modelQuality,PredictedConcentrations=c(NaN,NaN),ObservedandPredicted=c(NaN,NaN,NaN)))

}


#foobar<-modelExecution(chemical,numComp,subsetRate,calibration,specDataToModel,fitEval,fitFile,fitFileOut)
