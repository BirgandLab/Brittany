library(hydroTSM)
library(signal)
U1<-read.table("G:/claridge/09-04-2014/Upstream/SONTEKUP09-04-2014.LOG",skip=2)
U2<-read.table("G:/claridge/08-20-2014/Upstream/SONTEKClrUp2014-08-20.LOG",skip=2)
U3<-read.table("G:/claridge/2014-08-07/Upstream/SONTEKUp2014-08-07.LOG",skip=2)
U4<-read.table("G:/claridge/07-24-2014/upstream/SONTEK.LOG",skip=2)
U5<-read.table("G:/claridge/07-14-2014/UP/SONTEKUS07-14-2014.LOG",skip=2)
U6<-read.table("G:/claridge/07-09and10-2014/UP/SONTEKUS07-09-2014.LOG",skip=2)
Usontek<-rbind(U6,U5,U4,U3,U2,U1)
colnames(Usontek)<-c("SN","sample","YYYY","MM","DD","HH","mm","ss","flow","stage","MVel","TVol","zP","IndxV", 
                     "XsecA","Temp","SStat","CXV","CZV","LXV","RXV","Bat","Pitch","Roll","Submerged","zAcoustic",
                     "z3","Totv+","TotV-","ECell","SNRB1","SNRB2","SNRB3","SNRB4")
success<-paste(Usontek$YYYY,"-",Usontek$MM,"-",Usontek$DD," ",
               Usontek$HH,":",Usontek$mm,":",Usontek$ss,sep="")
Usontek$DateTime<-as.POSIXct(tz="EST",success)
remove(success,U1,U2,U3,U4,U5,U6)

M1<-read.table("G:/claridge/09-04-2014/Middlestream/SONTEKMD09-04-2014.LOG",skip=3)
M2<-read.table("G:/claridge/08-20-2014/Middlestream/SONTEKClrMd08-20-2014.LOG",skip=2)
M3<-read.table("G:/claridge/2014-08-07/Middlestream/SONTEKMD08-07-2014.LOG",skip=3)
M4<-read.table("G:/claridge/07-24-2014/middle/SONTEK.LOG",skip=2)
M5<-read.table("G:/claridge/07-14-2014/MD/SONTEKMD07-14-2014.LOG",skip=2)
M6<-read.table("G:/claridge/07-09and10-2014/MD/SONTEKMD07-10-2014.LOG",skip=2)
Msontek<-rbind(M6,M5,M4,M3,M2,M1)

colnames(Msontek)<-c("SN","sample","YYYY","MM","DD","HH","mm","ss","flow","stage","MVel","TVol","zP","IndxV",
                     "XsecA","Temp","SStat","CXV","CZV","LXV","RXV","Bat","Pitch","Roll","Submerged","zAcoustic",
                     "z3","Totv+","TotV-","ECell","SNRB1","SNRB2","SNRB3","SNRB4")

success<-paste(Msontek$YYYY,"-",Msontek$MM,"-",Msontek$DD," ",
               Msontek$HH,":",Msontek$mm,":",00,sep="")
Msontek$DateTime<-strptime(tz="EST",success,format="%Y-%m-%d %H:%M:%S")

remove(success,M6,M5,M4,M3,M2,M1)


D1<-read.table("G:/claridge/09-04-2014/Downstream/SONTEKDn2014-09-04.LOG",skip=2)
D2<-read.table("G:/claridge/08-20-2014/Downstream/SONTEKDn2014-08-20.LOG",skip=2)
D3<-read.table("G:/claridge/2014-08-07/Downstream/SONTEKDn2014-08-07.LOG",skip=2)
D4<-read.table("G:/claridge/07-24-2014/Downstream/SONTEK.LOG",skip=2)
D5<-read.table("G:/claridge/07-14-2014/DN/SONTEKDn2014-07-14.LOG",skip=2)
D6<-read.table("G:/claridge/07-09and10-2014/Dn/SONTEKDn2014-07-10.LOG",skip=2)
Dsontek<-rbind(D6,D5,D4,D3,D2,D1)


colnames(Dsontek)<-c("SN","sample","YYYY","MM","DD","HH","mm","ss","flow","stage","MVel","TVol","zP","IndxV",
                     "XsecA","Temp","SStat","CXV","CZV","LXV","RXV","Bat","Pitch","Roll","Submerged","zAcoustic",
                     "z3","Totv+","TotV-","ECell","SNRB1","SNRB2","SNRB3","SNRB4")

success<-paste(Dsontek$YYYY,"-",Dsontek$MM,"-",Dsontek$DD," ",
               Dsontek$HH,":",Dsontek$mm,":",Dsontek$ss,sep="")
Dsontek$DateTime<-as.POSIXct(tz="EST",success)
remove(success,D6,D5,D4,D3,D2,D1)

 




#clean up bits of data
#where stage is 0 replace with NAN
#Usontek$stage[Usontek$stage==0]<-NaN
#Msontek$stage[Msontek$stage==0]<-NaN
#Dsontek$stage[Dsontek$stage==0]<-NaN

#Usontek$MVel[Usontek$MVel==0]<-NaN
#Msontek$MVel[Msontek$MVel==0]<-NaN
#Dsontek$MVel[Dsontek$MVel==0]<-NaN

plot(Dsontek$DateTime,Dsontek$stage,type="l",col="blue",ylab="stage [m]",xlab="date")
points(Usontek$DateTime,Usontek$stage,type="l",col="red")
points(Msontek$DateTime,Msontek$stage,type="l")


plot(Usontek$DateTime,Usontek$flow,type="l",col="red",ylab="flow [m3 s-1]",xlab="date")
points(Msontek$DateTime,Msontek$flow,type="l")
points(Dsontek$DateTime,Dsontek$flow,type="l",col="blue")

window<-4

plot(Usontek$stage,ma(abs(Usontek$MVel),window),type="l",col="red",ylab="flow [m s-1]",xlab="stage [m]")
points(Msontek$stage,ma(abs(Msontek$MVel),window),type="l")
points(Dsontek$stage,ma(abs(Dsontek$MVel),window),type="l",col="blue")

plot(Usontek$MVel,type="l")
points(Usontek$stage,type="l",col="green")
points(despike(abs(Usontek$MVel),reference="median",n=4),type="l",col="red")


butt <- butter(2,0.2,type="low")
b <- filtfilt(butt,Usontek$MVel)

plot(Usontek$MVel,type="l",col="black")
points(b,type="l",col="red")
points(Usontek$stage,type="l",col="green")

