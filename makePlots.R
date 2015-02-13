#PLOT ANNUAL SCALE
plot(as.POSIXct(AP_R$x,origin="1970-01-01 00:00:00"),
     AP_R$y,type="l",col="gray",
     main=chemicals[chemical],
     ylim=c(min(AP_R$y),max( AP_R$y)),
     xlim=c(min( AP_R$x),max( AP_R$x)),
     xlab="date",
     ylab=paste(chemicals[chemical],"concentration",sep=" ")
)
butt <- butter(2,0.005,type="low")
points(as.POSIXct(AP_R$x,origin="1970-01-01 00:00:00"),
       filtfilt(butt,chemoGraphicTable[,chemical+2]),type="l",col="red",
       ylim=c(min(AP_R$y),max( AP_R$y)),
       xlim=c(min( AP_R$x),max( AP_R$x))
)
abline(h=0)
points(as.POSIXct(ModelRandom$ObservedAndPredicted[,3],origin="1970-01-01 00:00:00"),
       ModelRandom$ObservedAndPredicted[,2],
       ylim=c(min(AP_R$y),max( AP_R$y)),
       xlim=c(min(AP_R$x),max( AP_R$x)),
       col="gray"
)
points(as.POSIXct(ModelRandom$ObservedAndPredicted[,3],origin="1970-01-01 00:00:00"),
       ModelRandom$ObservedAndPredicted[,1],
       ylim=c(min(AP_R$y),max( AP_R$y)),
       xlim=c(min(AP_R$x),max( AP_R$x)),
       col="blue"
)
par(new=TRUE)
plot(as.POSIXct(AP_R$x,origin="1970-01-01 00:00:00"),
     chemoGraphicTable[,2],col="green",type="l",xaxt="n",yaxt="n",xlab="",ylab="",
      xlim=c(min( AP_R$x),max( AP_R$x))
      )
mtext(paste(iD[counter,9],iD[counter,8],iD[counter,3],iD[counter,4],iD[counter,5],iD[counter,6],sep=","))
axis(4)
#eventScale
plot(as.POSIXct(AP_R$x[event],origin="1970-01-01 00:00:00"),
     AP_R$y[event],type="l",col="gray",
     main=chemicals[chemical],
     ylim=c(min(AP_R$y),max( AP_R$y)),
     xlim=c(min( AP_R$x[event]),max( AP_R$x[event])),
     xlab="date",
     ylab=paste(chemicals[chemical],"concentration",sep=" ")
)

points(as.POSIXct(AP_R$x[event],origin="1970-01-01 00:00:00"),
       filtfilt(butt,chemoGraphicTable[event,chemical+2]),type="l",col="red",
       ylim=c(min(AP_R$y[event]),max( AP_R$y[event])),
       xlim=c(min( AP_R$x[event]),max(AP_R$x[event]))
)
points(useUsRealTIme,useUsChem,
       ylim=c(min(AP_R$y[event]),max( AP_R$y[event])),
       xlim=c(min(AP_R$x[event]),max( AP_R$x[event]))
)
abline(h=0)
par(new=TRUE)
plot(as.POSIXct(AP_R$x[event],origin="1970-01-01 00:00:00"),
     chemoGraphicTable[event,2],col="green",
     type="l",xaxt="n", yaxt="n",xlab="",ylab="",
     xlim=c(min(AP_R$x[event]),max( AP_R$x[event]))
)

mtext(paste(iD[counter,9],iD[counter,8],iD[counter,3],iD[counter,4],iD[counter,5],iD[counter,6],sep=","))
axis(4)