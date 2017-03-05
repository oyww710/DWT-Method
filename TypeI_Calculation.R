###Simulation#############
##LF###
library(multcomp)
library(mwt)
library(limma)
rm(list=ls())
setwd ("C:/Users/oyww710/Desktop/DWT")

Levene.test1 <- function (y, group) 
{
    group <- as.factor(group) 
    means <- tapply(y, group, mean, na.rm = TRUE)
    resp <- abs(y - means[group])
    return(c(anova(lm(resp ~ group))[,4][1],anova(lm(resp ~ group))[,5][1]))
} 


##BF###
Levene.test2 <- function (y, group) 
{
    group <- as.factor(group) 
    means <- tapply(y, group, median, na.rm = TRUE)
    resp <- abs(y - means[group])
    return(c(anova(lm(resp ~ group))[,4][1],anova(lm(resp ~ group))[,5][1]))
} 


Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p),lower.tail=FALSE)
  return(c(Xsq = Xsq, p.value = p.val))
}

###0908####


tF.normal.P<-function(m,n,a,b,t,h,N)
{

G<-c(rep(0,m),rep(1,n))
M1 <- matrix(rnorm(N*m,a,b),N,m)  #N*m#
M2<-matrix(rnorm(N*n,a+t,b+h),N,n) #N*n#
M<-cbind(M1,M2)
fit <- lmFit(M,design=G)
fit1<-ebayes(fit)
t_m<-as.numeric(fit1$p.value)
t_mw<- mwt(M,G)$pvalue

LRT_p<-t_n<-t<-F<-L_F<-L_F1<-DWT<-DWT_Adj<-tF<-tLF<-tBF<-SMT<-SMT1<-SMT2<-array()
for (i in 1:N)
{

s1i<-as.numeric(M1[i,]);
s2i<-as.numeric(M2[i,]);
	m1i<-mean(s1i);
	m2i<-mean(s2i);
ss1i1<-(sqrt((m+1)/(m-1))*(s1i-m1i))^2
ss1i<-(s1i-m1i)^2
ss2i<-(s2i-m1i)^2
S1i<-s1i^2
S2i<-s2i^2
	v1i<-var(s1i); 
	v2i<-var(s2i); 
     v<-var(c(s1i,s2i))
s1i_star<-s1i/sd(s1i);
s2i_star<-s2i/sd(s2i);
#####LRT############################################################
LRT<-(m+n)*log(v*((m+n-1)/(m+n)))-m*log(v1i*((m-1)/m))-n*log(v2i*((n-1)/n))
#LRT<--2*log(max1/max2)#
LRT_p[i]<-pchisq(LRT,2,lower.tail=FALSE) 

#########t on normalized data##########################
t_n[i]<-t.test(s1i_star,s2i_star,var.equal=T)$p.value
############################################################
dwt1<-t.test(ss1i,ss2i)$p.value;  #DWT#
dwt2<-t.test(ss1i1,ss2i)$p.value;  #DWT#
G1<-c(rep(0,length(s1i)),rep(1,length(s2i)))
Y<-c(s1i,s2i)

t[i]<-t.test(s1i,s2i)$p.value;  #Welch t test#
t1<-t.test(s1i,s2i)$p.value
F[i]<-var.test(s1i,s2i)$p.value; #F test#
F1<-var.test(s1i,s2i)$p.value;
L_F[i]<-Levene.test1(Y,G1)[2]   #LF#
L_F_1<-Levene.test1(Y,G1)[2]
L_F1[i]<-Levene.test2(Y,G1)[2]  #BF#
L_F_11<-Levene.test2(Y,G1)[2]
SMT[i]<-t.test(S1i,S2i)$p.value;
SMT1[i]<-dwt1
SMT2[i]<-dwt2
tF[i]<-Fisher.test(p = c(t1,F1))[[2]] ##tF## 
tLF[i]<-Fisher.test(p=c(t1,L_F_1))[[2]] ##tLF##
tBF[i]<-Fisher.test(p=c(t1,L_F_11))[[2]] ##tBF##
DWT[i]<-Fisher.test(p=c(t1,dwt1))[[2]] #
DWT_Adj[i]<-Fisher.test(p=c(t1,dwt2))[[2]] 

}


	tF=cbind(t,t_m,t_mw,t_n,F,L_F,L_F1,LRT_p,DWT,DWT_Adj,SMT,SMT1,SMT2,tF,tLF,tBF);

	return(tF);
}

PValue<-tF.normal.P(10,10,0,1,0,0,100000)

write.csv(PValue,"Power_Type_10_10_0_1_0_0.csv")


