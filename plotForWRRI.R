source("C:/Users/FBlab/Desktop/Brittany/plotForWRRI.env.R", echo=FALSE)
#y<-3
#
#Types 
# 1="original"
# 2="prunedO",
# 3="1stDer"
# 4="pruned1D",
# 5="turbComp"
# 6="prunedTC",
# 7="1stDerTurbComp"
# 8= "prunedTC1D"  
#chemical
#  1 = "DIC"  2 = "DOC"   3 = "MES"
#  4 = "NNO3" 5 = "PPO43" 6 = "Ptot"
#  7 = "SO4"  8 = "Turb"  9 = "CL"
#  10 = "NO2" 11 = "Ntot" 12 = "NTotFilt"
#  13 = "NH4" 14 = "Silica"


#NO3
makePlotsForWWRI(year=1,FileType=1,keepSkippedDates=1,chemical=4,numComp=3,Flow=Flow,EXP=1.8)
makePlotsForWWRI(year=1,FileType=2,keepSkippedDates=1,chemical=4,numComp=3,Flow=Flow,EXP=1.8)

#DOC
makePlotsForWWRI(year=1,FileType=1,keepSkippedDates=1,chemical=2,numComp=3,Flow=Flow,EXP=1.8)
makePlotsForWWRI(year=1,FileType=2,keepSkippedDates=1,chemical=2,numComp=3,Flow=Flow,EXP=1.8)

#Silica
makePlotsForWWRI(year=3,FileType=1,keepSkippedDates=1,chemical=14,numComp=3,Flow=Flow,EXP=1.8)





