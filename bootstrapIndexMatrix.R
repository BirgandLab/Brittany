#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#                                       TSBOOTSTRAP MATRIX
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#uses
#library(tseries)

bsIndMat<-function(realTime,numReps){
#bsIndMat(time,reps)
#here is a dummy function to return, but not print the result of the bootstrap 
#because we are just sending POSIX integer dates, the output is the date of the 
#data that should be included in the bootstrap model
    
    mystat<-function(y){
      invisible(y)  
    }


#this function needs a list of times, and number of bootstrap iterations
#and will return a matrix of block subsamples from that list of times
    #realTime<-originalmyData$realTime

#make sure there are no NANs
    realTime<-realTime[!is.na(realTime)]

#convert them to a numeric vector
    bsTime<-(as.numeric(realTime))

#run the bootstrapping procedure--it returns "$statistic", which is a matrix of each BS
#of dates, rows contain dates from each BS, columns individual dates
      #turns out, it would have worked about right if i didn't do anything...(ie stat=NULL sigh)
    totalBS<-tsbootstrap(bsTime,nb=numReps,statistic=mystat)

    indices<-totalBS$statistic

}
# 
# set.seed(1)
# 
# a<-bsIndMat(originalmyData$realTime[1:100],1000)
# plot(a[1,])
# for(i in 2:100){
#   points(a[i,])
# }
