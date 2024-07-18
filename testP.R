rm(list=ls())
time1<-Sys.time()
#source("powerN.R")
source("powerD.R")
source("DC.R")
source("select_h.R")

model=1
n=100
a=0
sd0=sqrt(0.25)

#result <- POWER_CALCULATIONn(model, n, sd0, a)
result <- POWER_CALCULATIONd(model, n, sd0, a)
print(result)
time2<-Sys.time()
print(time2 - time1)