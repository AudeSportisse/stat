source('PlotSimulation.R')

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
noise=0.8
Ns=500
dat=rnorm(n)
Xtrue=cbind.data.frame(dat,dat*1.1,dat*1.2,dat*1.3)

##Figure 2 (logit, MNAR)
resultUnivariate=ComparMNAR_Univariate(Xtrue,a,b,r,noise,Ns,"logit","MNAR",nbsim)
Plot(resultUnivariate)
##Figure 5 (logit, MAR)
resultMARUnivariate=ComparMNAR_Univariate(Xtrue,a,b,r,noise,Ns,"logit","MAR",nbsim) 
Plot(resultMARUnivariate,type="WithoutMask")
#Figure 7 (probit, MNAR)
resultProbitUnivariate=ComparMNAR_Univariate(Xtrue,a,b,r,noise,Ns,"probit","MNAR",nbsim) 
Plot(resultProbitUnivariate,type="WithoutMask")

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
noise=0.8
Ns=500
res=simu(n,p,r,sqrt(noise))
Xtrue=res$mu

##Figure 3 (MNAR)
resultMNARBivariate=ComparMNAR_Bivariate(Xtrue,a,b,a2,b2,r,noise,Ns,colbis=4,m1=1,m2=4,nbsim)
Plot(resultMNARBivariate)
##Figure 6 (MAR)
resultMARBivariate=ComparMNAR_Bivariate(Xtrue,a,b,a2,b2,r,noise,Ns,colbis=4,m1=10,m2=20,nbsim) 
Plot(resultMARBivariate,type="WithoutMask")

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
noise=0.5 
Ns=500
Xtrue=simu(n,p,r,sqrt(noise))$mu

##Figure 4 (for different noise values)
resultMNARMultivariate=ComparMNAR_Multivariate(Xtrue,a,b,r,noise,Ns,10,nbsim)
Plot(resultMNARMultivariate)
