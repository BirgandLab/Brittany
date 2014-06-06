library(oce)
library(hydroGOF)
library(pls)  #Load the pls package
Components<-4
iterations<-5
#subsetting<-1 #for random picks
#subsetting<-2 #for contiguous windows
#for (subsetting in c(2)){
  #******LOAD DATA SPECIFYING PATH AND FILE NAME
  path<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\" #Specify folder where data is located
  chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
          "NNH4",  "Ntot",  "NTotFilt",  "Silica",	"Turbidity");
  filename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")
  
  despikingRounds<-1 #best results when set to 1
  counter<-0
  Stats<-matrix(nrow=1000,ncol=33)
  colnames(Stats)<-c("dfSS","r2SS","RMSESS","NRMSESS","MSS","CISS","maxSS","minSS","meanSS","df","r2","RMSE","NRMSE","M","CI","max1","min1","mean1","dfw","r2w","RMSEw","NRMSEw","Mw","CIw","max1w","min1w","mean1w","max","min","chem","file","sampleSize","iteration")
  badDates<-matrix(nrow=20,ncol=300)
  differences<-matrix(nrow=20,ncol=300)
  #for (i in 224:241) { #loop through the chemical constituents
  for (i in 240:240){ #pick a chemical constituent
    #now all chem analyses are listed in different columns 
    #224CL  225NO2  2276NNO2	227NO3	228NNO3	229SO4	
    #230DOC	231DIC	232UV254	233PPO43	234Ptot	
    #235MES	236NNH4	237Ntot	238NTotFilt	239Silica	240Turb
    for (j in 1:1){ #loop through the different spectra
      
      data<-read.table(file=paste(path,filename[j],sep=""),sep=",",header=TRUE,skip=0)  #Import data as .csv file
      
      #******PULL OUT DATE AND CONVERT IT INTO A REAL DATE, THEN PUT IT BACK INTO THE MATRIX
      Date<-substr(data$DateScan, 1, 10) 
      T<-substr(data$DateScan, 12, 19)  
      T[T==""]="00:00:00"               
      D<-paste(Date,T,sep=" ")        
      D<-strptime(D, '%d/%m/%Y %H:%M:%S')
      data$DateScan<-D
      #*****SELECT A WINDOW OF DATA TO USE FOR ANALYSIS--I HAVE NO IDEA WHY THIS HAS BEEN DONE
      date1<-strptime("14/4/2011",'%d/%m/%Y') #beginning of data window (Why these?)
      date2<-strptime("30/1/2014",'%d/%m/%Y') #end of data (Why these dates?)
      
      allData<-data[,-(dim(data)[2]-20):-(dim(data)[2])]
      allData<-cbind(allData,data[,i])
      allData<-allData[complete.cases(allData[,2:dim(allData)[2]]),] # removes all the rows for which there is a NA--keeping time in there
      ChemConc<-data[complete.cases(data[i]),i]

      #lets go through and make progressively larger data sets
      #OR pick 100 point data sets and have them start at differt points?
      values<-c(25,50,100,200,300,400,500,600,700,800)
     # for (k in 800){ #number of values(*100) to use in modeling
        
     for (l in 1:iterations){#run through a bunch of times for stats fun
          counter<-counter+1
          allData<-data[,-(dim(data)[2]-20):-(dim(data)[2])]
          allData<-cbind(allData,data[,i])
          allData<-allData[complete.cases(allData[,2:dim(allData)[2]]),] # removes all the rows for which there is a NA--keeping time in there
          
          subdata<-allData
#SORT DATA FOR STRATIFICATION
FP<-subdata  #FP has data, OK, spectra and chem for a complete dataset
FP<-FP[order(FP[,dim(FP)[2]]),] #SORTED IN ASCENDING ORDER
ChemConc<-FP[,dim(FP)[2]]
minimum<-min(ChemConc)
maximum<-max(ChemConc)
#ChemPerc<-ChemConc/maximum
divisions<-200
#for each division, count the values, and pick some...
increment<-round((maximum-minimum)/divisions)
steps<-seq(minimum,maximum,increment) 
#OK simple approach isn't warranted
#because not uniformly distributed, break into %of maximum? and subset by %
#How about using native HIST function with Sturgess, scott and DF methods for breaking.
#Then, if possible take 1 point from each bin
histKyle<-hist(ChemConc,breaks=steps, plot=FALSE)
#histSturges<-hist(ChemConc,breaks="Sturges")
#histScott<-hist(ChemConc,breaks="Scott",plot=FALSE)
#histFD<-hist(ChemConc,breaks="FD",plot=FALSE)

Histogram<-histKyle

sampleRatio<-1/10
count<-0
for (q in seq(1,(length(Histogram$breaks)-1)*0.50,1)){
#for (q in seq((length(Histogram$breaks)-1)*0.25,(length(Histogram$breaks)-1),1)){
  inTheBin<-FP[(ChemConc>=Histogram$breaks[q])&(ChemConc<=Histogram$breaks[q+1]),]
  sampleNum<-sum(ceiling(Histogram$counts[q]*sampleRatio))
  someNumber<-inTheBin[sample(1:dim(inTheBin)[1],sampleNum,replace=FALSE),]
  if(count==0){
    calibration<-someNumber
  }
  if(count>0){
    calibration<-rbind(calibration,someNumber)
  }
  count<-count+1
}
remove(someNumber)


# calibration<-as.data.frame(matrix(0,nrow=length(Histogram$breaks),ncol=length(FP)))
# #calibration<-rep(0,(length(Histogram$breaks)-1),dim(FP)[2])
# for (q in seq(1,(length(Histogram$breaks)-1),1)){
#   inTheBin<-FP[(ChemConc>=Histogram$breaks[q])&(ChemConc<=Histogram$breaks[q+1]),]
#    if(dim(inTheBin)[1]>0){
#     calibration[q,]<-inTheBin[sample(1:dim(inTheBin)[1],1,replace=FALSE),]
#   }
# }

  


        
          FP<-calibration[calibration[,1]>0,]
        
          colnames(FP)[dim(FP)[2]]<-"Conc" 
          colnames(allData)[dim(allData)[2]]<-"Conc"

          parameter<-chem[i-223]
          
          FPConc<-FP[,dim(FP)[2]]
          FPDate<-FP[,1]
#plot(FPConc)          

          
          #*******STORE DATE AND CONCENTRATION VECTOR, STRIP DATE FROM THE MATRIX
          
          realtimeSubSample<-FP$DateScan
          realTimeAllData<-allData$DateScan #complete cases--have actual chem data
          FP<-FP[,-1:-2] #Remove 1st two columns of data (date-time and OK) 
          allData<-allData[,-1:-2]
          
          ConcSubSample<-as.matrix(FP$Conc) # creates a vector of CHEMICAL CONCENTRATION DATA ONLY
          ConcAllData<-allData$Conc # creates a vector of CHEMICAL CONCENTRATION DATA ONLY
          FP<-FP[,-(dim(FP)[2])] #drop the concentration data from the end??
          allData<-allData[,-(dim(allData)[2])] #drop the concentration data from the end??
          
          xlim<-c(min(ConcAllData),max(ConcAllData))
          ylim<-c(min(ConcAllData),max(ConcAllData))
          
          FP<-data.matrix(FP) #Convert to data matrix
          allData<-data.matrix(allData)
          
    #RUN PLSR MODEL FOR SUBSAMPLE
          fitSubSample<-plsr(ConcSubSample~data.matrix(FP),ncomp=Components,validation="CV")  #PLSR model to predict NO3-N with cross validation
    #RUN PLSR MODEL FOR WHOLE DATSET
          fitWholeModel<-plsr(ConcAllData~allData,ncomp=15,validation="CV")
    #predict values for data with PLSR model
          SubSampleFitWithSubSample<-predict(fitSubSample,FP,ncomp=Components,type=c("response")) #Predict NO3-N concentrations based on PLSR model

          ConcP.SubSample.SubSample<-as.data.frame(matrix(0,1,dim(FP)[1])) #Create data frame for predicted NO3-N values
          ConcP.SubSample.SubSample<-SubSampleFitWithSubSample[1:dim(FP)[1]]  #Insert predicted NO3-N values into data frame
          #print(summary(fit))  #See summary of PLSR model to choose number of components
    
  #predict values for whole dataset from subsample PLSR model       
          allDataFitWithSubSample<-predict(fitSubSample,allData,ncomp=Components,type=c("response")) #USE MODEL TO FIT ALL DATA (SOME NOT INCLUDED IN MODE L)
          allDataFitWithAllData<-predict(fitWholeModel,allData,ncomp=Components,type=c("response"))

          ConcP.allData.SubSample<-as.data.frame(matrix(0,1,dim(allData)[1])) #Create data frame for predicted NO3-N values
          ConcP.allData.SubSample<-allDataFitWithSubSample[1:dim(allData)[1]]  #Insert predicted NO3-N values into data frame
        
          ConcP.allData.allData<-as.data.frame(matrix(0,1,dim(allData)[1])) #Create data frame for predicted NO3-N values
          ConcP.allData.allData<-allDataFitWithAllData[1:dim(allData)[1]]  #Insert predicted NO3-N values into data frame



      #Calculate statistics for how well PLSR model predicts lab Values
      #SUBSAMPLE PLSR predicts SUBSAMPLE lab data
          SubSample.SubSample.FitStats<-lm(as.matrix(ConcP.SubSample.SubSample)~ConcSubSample) #Linear regression of predicted and lab NO3-N values
          subsetFit<-summary(SubSample.SubSample.FitStats)
          confInt<-predict(lm(ConcSubSample~as.matrix(ConcP.SubSample.SubSample)),interval='prediction',level=0.95)

          
      #SUBSAMPLE PLSR predicts ALLDATA lab data    
          AllData.SubSample.FitStats<-lm(as.matrix(ConcP.allData.SubSample)~ConcAllData) #Linear regression of predicted and lab NO3-N values
          subsetFit2<-summary(AllData.SubSample.FitStats)
          confIntR<-predict(lm(ConcAllData~as.matrix(ConcP.allData.SubSample)),interval='prediction',level=0.95)
          
      #Whole model predicts all data
          AllData.AllData.FitStats<-lm(ConcAllData~as.matrix(ConcP.allData.allData)) #Linear regression of predicted and lab NO3-N values
          subsetFit3<-summary(AllData.AllData.FitStats)
          confInt3<-predict(lm(ConcAllData~as.matrix(ConcP.allData.allData)),interval='prediction',level=0.95)

# 
# 
predicteds<-ConcP.SubSample.SubSample
observeds<-FPConc

predicted<-ConcP.allData.SubSample
observed<-ConcAllData

plot(observed, predicted,col='black',main="comparison of models")
abline(AllData.SubSample.FitStats,col='black')


points(ConcAllData,ConcP.allData.allData,col='gray')
abline(AllData.AllData.FitStats,col='gray')

points(observeds,predicteds,col='green')
abline(SubSample.SubSample.FitStats,col='green')

readline()

#STATS FOR FITTING SUBSET TO SUBSET
          Stats[counter,1]<-subsetFit$df[2]
          Stats[counter,2]<-subsetFit$r.squared
          Stats[counter,3]<-rmse(obs=ConcSubSample,sim=as.matrix(ConcP.SubSample.SubSample))
          Stats[counter,4]<-nrmse(obs=ConcSubSample,sim=as.matrix(ConcP.SubSample.SubSample))
          Stats[counter,5]<-subsetFit$coefficients[2,1] #slope
          Stats[counter,6]<-mean(confInt[,3]-confInt[,1])
          Stats[counter,7]<-max(ConcSubSample)
          Stats[counter,8]<-min(ConcSubSample)
          Stats[counter,9]<-mean(ConcSubSample)
          
#STATS FOR FITTING ALL DATA TO SUBSET MODEL
          Stats[counter,10]<-subsetFit2$df[2]
          Stats[counter,11]<-subsetFit2$r.squared
          Stats[counter,12]<-rmse(obs=ConcAllData,sim=(ConcP.allData.SubSample))
          Stats[counter,13]<-nrmse(obs=ConcAllData,sim=(ConcP.allData.SubSample))
          Stats[counter,14]<-subsetFit2$coefficients[2,1] #slope
          Stats[counter,15]<-mean(confIntR[,3]-confIntR[,1])
          Stats[counter,16]<-max(ConcP.allData.SubSample)
          Stats[counter,17]<-min(ConcP.allData.SubSample)
          Stats[counter,18]<-mean(ConcP.allData.SubSample)

#STATS FOR FITTING ALL DATA TO ALL DATA MODEL
          Stats[counter,19]<-subsetFit3$df[2]
          Stats[counter,20]<-subsetFit3$r.squared
          Stats[counter,21]<-rmse(obs=ConcAllData,sim=(ConcP.allData.allData))
          Stats[counter,22]<-nrmse(obs=ConcAllData,sim=(ConcP.allData.allData))
          Stats[counter,23]<-subsetFit3$coefficients[2,1] #slope
          Stats[counter,24]<-mean(confInt[,3]-confInt[,1])
          Stats[counter,25]<-max(ConcP.allData.allData)
          Stats[counter,26]<-min(ConcP.allData.allData)
          Stats[counter,27]<-mean(ConcP.allData.allData)

          Stats[counter,28]<-max(ConcAllData)
          Stats[counter,29]<-min(ConcAllData)

          Stats[counter,30]<-i
          Stats[counter,31]<-j
          #Stats[counter,32]<-k
          Stats[counter,33]<-l
    
          


        } #for each of the different spectra (make sure file read is inside here)
        
      #}#number of values(*100) to use in modeling
    }#run through a bunch of times for stats fun
    
  }#for each constituent

#}#for each subset approach

