pathD<-"C:/Users/birgand/Documents/Publications/2014/PLSR_Brittany/" #Specify folder where data is located
filename<-"Fingerprint_R_1D.csv"  #Specify file where data is located
data<-read.table(file=paste(pathD,filename,sep=""),sep=",",skip=1)  #Import data as .csv file
#data<-data1[-nrow(data1),]
attach(data)

colnames(data) <- c( 'Date/Time', 'Status') #Add column names
colnames(data)[3:(length(data)-2)]<-c(seq(200,750,2.5)) #Add column names

FP<-data[,-1:-10] #Remove NAs at low wavelengths
FP<-FP[,-(dim(FP)[2]-13):-dim(FP)[2]] #Remove NO3-N values and NAs at high wavelengths
FP<-data.matrix(FP) #Convert to data matrix
NO3<-data.matrix(data[dim(data)[2]-1]) #Make matrix of NO3-N values
Sal<-data.matrix(data[dim(data)[2]]) #Make matrix of salinity values if included in spreadsheet

library(pls)  #Load the pls package

fit<-plsr(NO3~data.matrix(FP),ncomp=20,validation="CV")  #PLSR model to predict NO3-N with cross validation
summary(fit)  #See summary of PLSR model to choose number of components


Pfit<-predict(fit,FP,ncomp=14,type=c("response")) #Predict NO3-N concentrations based on PLSR model

NO3P<-as.data.frame(matrix(0,1,dim(data)[1])) #Create data frame for predicted NO3-N values
NO3P<-as.data.frame(Pfit[1:108])  #Insert predicted NO3-N values into data frame

plot(NO3,as.matrix(NO3P)) #Compare predicted and lab values of NO3-N

fit2<-lm(NO3~as.matrix(NO3P)) #Linear regression of predicted and lab NO3-N values
summary(fit2) #Summary of linear regression

library(lars) #Load the lars package for Lasso and stepwise regression testing
mod1<-lars(FP,NO3,normalize=FALSE,type=c("lasso"),intercept=TRUE) #Lasso model to predict NO3-N values
mod1cv<-cv.lars(FP,NO3,K=10,normalize=FALSE,type=c("lasso"),intercept=TRUE) #Cross validation of Lasso model
sAtBest<-mod1cv$index[which.min(mod1cv$cv)] #Choose Lasso model based on cross validation
Pfit1<-as.data.frame(predict(mod1,FP,s=sAtBest,mode=c("fraction"))) #Predict NO3-N values based on chosen Lasso model

fit3<-lm(NO3~as.matrix(Pfit1[4])) #Linear regression lab and Lasso model predicted NO3-N values
summary(fit3) #Summary of linear regression
sqrt(min(mod1cv$cv))  #Calculate the Root Mean Square Error of Prediction (RMSEP) of the Lasso Model
coef1<-as.data.frame(coef(mod1,s=sAtBest,mode=c("fraction"))) #Determine the number of wavelengths included in the Lasso model
sum(coef1!=0) #Count the number of wavelengths included in the Lasso Model

mod2<-lars(FP,NO3,normalize=FALSE,type=c("stepwise"),intercept=TRUE)  #Stepwise regression model to predict NO3-N values
mod2cv<-cv.lars(FP,NO3,K=10,normalize=FALSE,type=c("stepwise"),intercept=TRUE)  #Cross validation of stepwise regression model
sAtBest2<-mod2cv$index[which.min(mod2cv$cv)]  #Choose stepwise regression model based on chosen Lasso model
Pfit2<-as.data.frame(predict(mod2,FP,s=sAtBest2,mode=c("step")))  #Predict NO3-N values based on chosen stepwise regression model

fit4<-lm(NO3~as.matrix(Pfit2[4])) #Linear regression of lab and stepwise regression model prediction NO3-N values
summary(fit4) #Summary of linear regression
sqrt(mod2cv$cv[sAtBest2]) #Calculate the RMSEP of the stepwise regression model
coef2<-as.data.frame(coef(mod2,s=sAtBest2,mode=c("step")))  #Determine the number of wavelengths included in the stepwise regression model
sum(coef2!=0) #Count the number of wavelengths included in the stepwise regression model

library(HH)
drop1<-matrix(0,dim(coef1)[1],1)
data1<-FP

for (i in 1:dim(coef1)[1]){
  if (coef1[i,]==0){
    drop1[i,]<-i
  }
}

row_sub = apply(drop1, 1, function(row) all(row !=0 ))
drop1<-drop1[row_sub,]
data1<-subset(data1,select=-c(drop1))
test1<-lm(NO3~data1)
vif(test1)

drop2<-matrix(0,dim(coef2)[1],1)
data2<-FP

for (i in 1:dim(coef2)[1]){
  if (coef2[i,]==0){
    drop2[i,]<-i
  }
}

row_sub = apply(drop2, 1, function(row) all(row !=0 ))
drop2<-drop2[row_sub,]
data2<-subset(data2,select=-c(drop2))
test2<-lm(NO3~data2)
vif(test2)

cex.lab=3
cex.axis=3
windowsFonts(times=windowsFont("Times"))
x11(width = 12, height = 10)
par(mar=c(6,6.5,1.5,1.5),family="times")
plot(as.matrix(NO3P),NO3,
     type="p",
     pch=16,
     cex=2,
     #col=tim.colors(175)[round(CDOM)],  #the whole trick is here
     lwd=2,                              #I chose tim.colors as a 
     yaxt="n",                           #nice color theme
     xaxt="n",
     xlim=c(0,5),
     ylim=c(0,5),
     yaxs="i",
     xaxs="i",        
     bty="n",        
     lty=1,
     xlab="",
     ylab=expression(Lab~Measured~NO[3]-N~(mg~L^-1)),
     cex.lab=cex.lab,
     cex.axis=cex.axis)
par(new=TRUE)
axis(2,cex.axis=cex.lab)
par(mgp=c(3,1.5,0),new=TRUE)
axis(1,cex.axis=cex.lab)
par(new=TRUE)
abline(0,1,lwd=2)
text(0.75,4.5,paste("R�=",signif(summary(fit2)$r.squared,digits=3),sep=""),cex=3)
title(xlab=expression(Predicted~NO[3]-N~(mg~L^-1)),cex.lab=3,mgp=c(5,1,0))