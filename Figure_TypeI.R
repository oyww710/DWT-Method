#setwd ("G:/prostate_cancer_microarray_data")
setwd ("D:/DWT");

# x is a sclar
# y is a vector or matrix
# n = length(y)
##############################################################################
xgrt.rate<-function(a, y, n){(length(which(y<a)))/n;} #faster than mean(a>y)

minp<-function(p1, p2)
{
	pmn<-.5*((p1+p2)-abs(p1-p2));
	return (pmn);
	#much faster than the inner pmin(); 
}

xgrt2.rate<-function(a, y1, y2, n)
{
	abf=1-sqrt(1-a); #Bonferroni 
	pwr=length(which( minp(y1, y2)<abf ))/n;
	return(pwr);
}  

plot.PWR<-function(T, methods, colv, lwdv, ltyv, x.mk, y.mk, main.txt)
{
	nrep=dim(T)[1]; #number of replicates (rows of T)
	nmethods = length(methods);
	npoints = 500;

	a<- (seq(0, 0.05, length.out=npoints));
	d= 1.96*sqrt(a*(1-a)/nrep);
	Uc<-(a+d);
	Lc<-(a-d);

	P<-T[, c(1, 2, 3, 4, 14, 15,9)];
	
	y1<-T[,1];
	y2<-T[,6];
	
	ePWR<-matrix(0, nrow=npoints, ncol=nmethods);
	for(i in c(1:(nmethods-1)))
	{
		ePWR[,i] <- sapply(a, xgrt.rate, P[,i], nrep);
	}
	ePWR[,nmethods]<-sapply(a, xgrt2.rate, y1, y2, nrep);
	colnames(ePWR)<-methods;


	plot(a, ePWR[,1], type="l", pch=21, lty=ltyv[1], lwd=lwdv[1], col=colv[1], 
		xlim=c(0, 0.05), ylim=c(0,0.05),axes=FALSE, ann=FALSE);

	points(a, ePWR[,2], type="l", pch=21, lty=ltyv[2], lwd=lwdv[2], col=colv[2])
	points(a, ePWR[,3], type="l", pch=22, lty=ltyv[3], lwd=lwdv[3], col=colv[3])
	points(a, ePWR[,4], type="l", pch=23, lty=ltyv[4], lwd=lwdv[4], col=colv[4])
	points(a, ePWR[,5], type="l", pch=24, lty=ltyv[5], lwd=lwdv[5], col=colv[5])
	points(a, ePWR[,6], type="l", pch=25, lty=ltyv[6], lwd=lwdv[6], col=colv[6])
	points(a, ePWR[,7], type="l", pch=19, lty=ltyv[7], lwd=lwdv[7], col=colv[7])
	points(a, ePWR[,8], type="l", pch=18, lty=ltyv[8], lwd=lwdv[8], col=colv[8])

	polygon(c(a, sort(a, decreasing = TRUE)), c(Lc, sort(Uc, decreasing = TRUE)), 
		col = colv[1+nmethods], border = colv[1+nmethods]);

	points(a, ePWR[,1], type="l", pch=21, lty=ltyv[1], lwd=lwdv[1], col=colv[1]);
	points(a, ePWR[,2], type="l", pch=21, lty=ltyv[2], lwd=lwdv[2], col=colv[2])
	points(a, ePWR[,3], type="l", pch=22, lty=ltyv[3], lwd=lwdv[3], col=colv[3])
	points(a, ePWR[,4], type="l", pch=23, lty=ltyv[4], lwd=lwdv[4], col=colv[4])
	points(a, ePWR[,5], type="l", pch=24, lty=ltyv[5], lwd=lwdv[5], col=colv[5])
	points(a, ePWR[,6], type="l", pch=25, lty=ltyv[6], lwd=lwdv[6], col=colv[6])
	points(a, ePWR[,7], type="l", pch=19, lty=ltyv[7], lwd=lwdv[7], col=colv[7])
	points(a, ePWR[,8], type="l", pch=18, lty=ltyv[8], lwd=lwdv[8], col=colv[8])

	axis(1, las = 1, at = x.mk, lab = sprintf("%g", x.mk), font.axis = 2);
	axis(2, las = 1, at = y.mk, lab = sprintf("%g", y.mk), font.axis = 2);

	legend("bottomright", methods, col=colv, lty=ltyv, lwd=lwdv, ncol=2, bg="white", 
		text.font=2, xjust = 0.5, cex=1, title.adj=0.2);

	title( main=main.txt, font.main=2,	cex.lab = 1, srt=1);
	title(xlab="Significance Level",font.lab=2);
	title(ylab="False Positive Rate",font.lab=2);
}



tiff('Figure2.tiff', 
	width = 720, height = 640, 
	units = "px", pointsize = 12, bg = "white")

par(mfrow = c(2, 2), mgp=c(2.5, .6, 0), pty = "m", bty='n', 
	pin=c(0.5, 0.5), mai=c(0.85, 0.85, .4, .1), 
	omi=c(.12, .12, .12, .12), 
	font.axis=2, font.lab=2, cex=1.025)


methods<-c("WT","MT","MWT", "STSD","FWT","IMVT","DWT","SMWT");
colv<-c("black", "black", "red", "red", "blue", "blue", "orange", "black", "gray75");
lwdv<-c(1, 1, 2, 2, 3, 3, 4, 4);
ltyv<-c(1, 2, 3, 1, 2, 3, 2, 3);

x.mk<-y.mk<-seq(0, 0.05, 0.01);

###(a) Normal: n1=n2=5
main.txt<-"(a) n1=n2=5"; #expression("(a)"~n[1]~"="~n[2]~"="~5);
T<- read.csv(file="TypeI_Normal_5.csv",head=TRUE,sep=",")
plot.PWR(T, methods, colv, lwdv, ltyv, x.mk, y.mk, main.txt);

###(b) Normal: n1=n2=10
main.txt<-"(b) n1=n2=10"; #expression("(b)"~n[1]~"="~n[2]~"="~10);
T<- read.csv(file="Type_Norm_10.csv",head=TRUE,sep=",")
plot.PWR(T, methods, colv, lwdv, ltyv, x.mk, y.mk, main.txt);

###(c) Normal: n1=n2=20
main.txt<-"(c) n1=n2=20"; #expression("(c)"~n[1]~"="~n[2]~"="~20);
T<- read.csv(file="TypeI_Normal_20.csv",head=TRUE,sep=",")
plot.PWR(T, methods, colv, lwdv, ltyv, x.mk, y.mk, main.txt);

###(d) Normal: n1=n2=40
main.txt<-"(d) n1=n2=40"; #expression("(d)"~n[1]~"="~n[2]~"="~40);
T<- read.csv(file="Power_Type_10_10_0_1_0_0.csv",head=TRUE,sep=",")
plot.PWR(T, methods, colv, lwdv, ltyv, x.mk, y.mk, main.txt);




dev.off()


