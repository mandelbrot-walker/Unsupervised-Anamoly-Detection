
library(kernlab)

ncentrality2<-read.csv("NCentrality_p2p_Gnutella04.csv")
ncentrality2<-within(ncentrality2, rm(X)) 


kpc <- kpca(~.,data=ncentrality2,kernel="rbfdot",
            kpar=list(sigma=0.2),features=2)

pcv<-pcv(kpc)

plot(rotated(kpc),col=rainbow(8),
     xlab="1st Principal Component",ylab="2nd Principal Component")