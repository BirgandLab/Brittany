library(oce) #used for despiking
library(hydroGOF) #forgot what I used this for
library(pls)  #Load the pls package

#******Specify file paths and names
path<-"C:\\Users\\FBLab\\Desktop\\workHere\\Data\\" #Specify folder where data is located
fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/"#fiteval_out.txt"
fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
filename<-c("OriginalBrittany.csv" ,"Brittany1stDerative.csv","TubidityCompensatedBrittany.csv","TurbidityCompensated1stDerivativeBrittany.csv")

#store names for the lab analytes
Chem<-c("CL", "NO2", "NNO2","NO3","NNO3","SO4","DOC","DIC","UV254", "PPO43","Ptot", "MES",
        "NNH4",  "Ntot",  "NTotFilt",  "Silica",  "Turbidity");

#read the data specified by the vector filename
counter<-1
Components<-matrix(nrow=70,ncol=4)


#load data for checking number of components on
for (chem in 1:17){
  #all chem analytes are listed in different columns 
  #1 CL     2 NO2   3 NNO2      4 NO3      5 NNO3    6 SO4	
  #7 DOC	  8 DIC	  9 UV254	    10 PPO43	11 Ptot	  12 MES	
  #13 NNH4  14 Ntot	15 NTotFilt	16 Silic	17 Turb
  
  #open a file to make a 4 panel plot (one for each of the 4 fingerprint files)
  jpeg(file=paste(path,"\\images\\",".",Chem[chem],"1.jpg",sep=""))
  par(mfrow=c(2,2))
        for (fn in 1:4){#for each filename
        #load the data
            myData<-loadDataFile(path,filename[1])
          #data are returned in a list myData$fingerPrints
          #                            myData$realTime
          #                            myData$ChemData
        
        #pull out the column of chem data of interest
            ChemConc<-as.matrix(myData$ChemData[,chem])                       #take out the only one we care about
                      fp<-cbind(myData$fingerPrints,myData$realTime,as.matrix(ChemConc))      #bind the two matrices together to determine complete cases
                      fp<-fp[complete.cases(fp[,2:dim(fp)[2]]),]              #just keep the fingerprints that have analytical data
            
                      ChemConc<-as.matrix(fp[,dim(fp)[2]])                    #pull chem back out
                      fp<-fp[,-dim(fp)[2]]                                    #pull off the fingerprint
        
                      realTime<-as.matrix(fp[,dim(fp)[2]])                    #pull real time off (now the last column)
                      fp<-fp[,-dim(fp)[2]]                                    #pop it off the end of the dataframe
                      Comps<-30    
        
        #calculate the PLSR model for the available data
        
                #subset if you like (comment out if not subset)
                      subset<-densityDependentSubset(ChemConc,realTime,fp,0.5,TRUE)
                      if(length(subset$ChemConc)<=39){
                            Comps<-round(length(subset$ChemConc)*0.5)
                          }  
                      modelRMSEP<-RMSEP(plsr(subset$ChemConc~data.matrix(subset$fingerPrint),ncomp=Comps,validation="CV"))
                          
                     
       
                    #  modelRMSEP<-RMSEP(plsr(ChemConc~data.matrix(fp),ncomp=30,validation="CV"))
        #and pull out the RMSEP values for each component
#*****Thishas a probelm the comlet cases might get the x axis off by removing a value, but not keeping the row count correct.
#could be fixed by 1.)avoiding nans in the RMSEP
#or adding a column of index variables so the complete cases keeps the index even though it throws out a row of bad data

                      rmsepIndex<-complete.cases(as.matrix(modelRMSEP$val[2,1,]))
                      rmsepValues<-modelRMSEP$val[2,1,rmsepIndex]
                      nComps<-min(which(rmsepValues==min(rmsepValues)))-1
                    
            
                        plot(0:(length(rmsepValues)-1),rmsepValues,xlab=c("number of comoponents"),ylab=c("RMSEP"),main=paste(filename[fn], Chem[chem],sep=" " ))
                        points(nComps,rmsepValues[nComps+1],col="green")
                        if(nComps>15){
                        nComps2<-which(abs(diff(rmsepValues))==min(abs(diff(rmsepValues))))
                        points(nComps2,rmsepValues[nComps2],col="red")   
                        }
                
            
            #find a place where the diff is minimized try that as a logical breakpoint
                      #}
            
            
            Components[counter,1]<-nComps
            Components[counter,2]<-fn
            Components[counter,3]<-chem
           if(nComps>15){
            Components[counter,4]<-nComps2
           }
          counter<-counter+1
        #}
        
        } #for each file
  dev.off()
} #for each chemical