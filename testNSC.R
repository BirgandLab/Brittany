fitPath<-"C:/Users/FBLab/Downloads/FITEVAL2_win/FITEVAL2_win/"#fiteval_out.txt"
fitEval<-paste(fitPath,"fiteval",sep="")
fitFile<-paste(fitPath,"PLSR.in",sep="")
fitFileOut<-paste(fitPath,"PLSR_out.txt",sep="")
totalIterations<-20
fitQuality<-matrix(nrow=totalIterations,ncol=14)
colnames(fitQuality)<-c("calibrationVG","calibrationG","calibratinA","calibrationB","calibFitVG","calibFitG","calibFitA","calibFitB","allDataVG","allDataG","allDataA","allDataB","i","j")
Stats<-matrix(nrow=totalIterations,ncol=17)
colnames(Stats)<-c("n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","n","r2","rmse","nrmse","slope","i","j")


oandp<-cbind(seq(1:500),seq(1:500)+rnorm(500,mean=0,sd=10))

fitQuality<-matrix(ncol=8,nrow=20)
#values<-c(1,250,500,1)
values<-c(1,100,200,300,400,500)
#values<-c(1,50,100,150,200,250,300,350,400,450,500)
for (i in 1:(length(values)-1)){
  subset<-values[i]:values[i+1]
fitQuality[i,1:4]<-OB(oandp[subset,1],oandp[subset,2],fitEval,fitFile,fitFileOut)
fitQuality[i,5:7]<-c(length(subset),values[i],values[i+1])
}
  
