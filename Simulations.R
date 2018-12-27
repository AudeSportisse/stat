######
#Univariate
######

source('Univariate.R') #Paralell computation with 20 clusters
set.seed(4)
n=100
p=4
r=1
a=3
b=0
nbsim=50
bruit=0.8
Ns=500
dat=rnorm(n)
Xtrue=cbind.data.frame(dat,dat*1.1,dat*1.2,dat*1.3)
resultProbitUnivariate=ComparMNAR_Univariate(Xtrue,a,b,r,bruit,Ns,"probit","MNAR",nbsim)
resultUnivariate=ComparMNAR_Univariate(Xtrue,a,b,r,bruit,Ns,"logit","MNAR",nbsim)
resultMARUnivariate=ComparMNAR_Univariate(Xtrue,a,b,r,bruit,Ns,"logit","MAR",nbsim)


######
#Bivariate
######

source('Bivariate.R') #Paralell computation with 20 clusters
set.seed(4)
n=100
p=50
r=4
a=3
b=0
a2=2
b2=1
nbsim=50
bruit=0.8
Ns=500
res=simu(n,p,r,sqrt(bruit))
Xtrue=res$mu
resultMARBivariate=ComparMNAR_Bivariate(Xtrue,a,b,a2,b2,r,bruit,Ns,colbis=4,m1=10,m2=20,nbsim) #MAR mechanism
resultMNARBivariate=ComparMNAR_Bivariate(Xtrue,a,b,a2,b2,r,bruit,Ns,colbis=4,m1=1,m2=4,nbsim) #MNAR mechnanism

#######
#Multivariate
#######

source('Multivariate.R') #Paralell computation with 20 clusters
source('General_tools.R')
set.seed(4)
n=100
p=20
r=4
a=3
b=0
nbsim=50
bruit=0.5 
Ns=500
Xtrue=simu(n,p,r,sqrt(bruit))$mu
result2=ComparMNAR_Multivariate(Xtrue,a,b,r,bruit,Ns,10,nbsim)
