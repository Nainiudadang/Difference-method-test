
########## Difference-based test methods (DBT_Gau, DBT_Gen, DBT_Perm)

diff.fn = function(n,lmd,y,method=1){
  # This is the test function using the difference based method
  # n is sample size; lmd is the parameter to control the bandwidth; 
  # y is the response
  # two methods for the bandwidth selection
  if(method==1){m =round((1+lmd - sqrt(lmd^2+2*lmd))*sqrt(14*n))};
  if(method==2){r <- optimize(beta.v,interval=c(0.7,0.85),n=n,y=y)$minimum;
  m=round(n^r)}
  N=n*m - m*(m+1)/2;
  wk = seq(n-1,n-m)/N;   
  dk = (seq(1,m)/n)^2;
  dwb=sum(wk*dk);
  skn = vector();
  for(k in 1:m)
  { 
    skn[k] = sum( (y[(1+k):n] - y[1:(n-k)])^2); 
  }
  sk = skn/(2*(n-seq(1:m)));
  beta.h =sum(wk*sk*(dk - dwb))/sum(wk*(dk-dwb)^2); ## estimated beta hat
  sg2.h = sum(wk*sk) - beta.h*dwb;  ## estimated sigma^2
  test.Gau= beta.h/sqrt((15/28*n^4/m+(45/4)*(n^5)/(m^3))*sg2.h^2/(n*m-m*(m+1)/2)^2) ## Using Gaussian Estimator
  e4.Evan = m4.evan.f(n,y);
  test.Gen= beta.h/sqrt((15/56*(e4.Evan - sg2.h^2)*n^4/m+45/4*(n^5)/(m^3)*sg2.h^2)/(n*m-m*(m+1)/2)^2) ## Using the General Estimator
  return(list(test.Gau=test.Gau, test.Gen=test.Gen, beta.h=beta.h, sg2.h=sg2.h,m=m))
}

beta.v <- function(r,n,y){
  m <- n^r
  N <- n*m - m*(m+1)/2;
  wk = seq(n-1,n-m)/N;   
  dk = (seq(1,m)/n)^2;
  dwb=sum(wk*dk);
  skn = vector();
  for(k in 1:m)
  { 
    skn[k] = sum( (y[(1+k):n] - y[1:(n-k)])^2); 
  }
  sk = skn/(2*(n-seq(1:m)));
  beta.h =sum(wk*sk*(dk - dwb))/sum(wk*(dk-dwb)^2); ## estimated beta hat
  sg2.h = sum(wk*sk) - beta.h*dwb;
  MSE <- (15/28*n^4/m+(45/4)*(n^5)/(m^3))*sg2.h^2/(n*m-m*(m+1)/2)^2 + (beta.h*(m+1)/(2*n-m-1))^2
  return(MSE)
}

m4.evan.f = function(n,y){
  ## This function is used to estimate the fourth moment using the method proposed by Evans and Jones (2008).
  ga = NULL
  i=1;ga[1] = (y[i]- y[i+1])*(y[i]- y[i+2])*(y[i]- y[i+3])*(y[i]- y[i+4])
  i=2;ga[2] = (y[i]- y[i-1])*(y[i]- y[i+1])*(y[i]- y[i+2])*(y[i]- y[i+3])
  for(i in 3:(n-2)){
    ga[i] = (y[i]- y[i-2])*(y[i]- y[i-1])*(y[i]- y[i+1])*(y[i]- y[i+2])
  }
  i=n-1;ga[n-1] = (y[i]- y[i-3])*(y[i]- y[i-2])*(y[i]- y[i-1])*(y[i]- y[i+1])
  i=n;ga[n] = (y[i]- y[i-4])*(y[i]- y[i-3])*(y[i]- y[i-2])*(y[i]- y[i-1])
  gam.evan = sum(ga)/n
  return(gam.evan)
}

diff.fn.b = function(n,m,y){
  ## This function calculates beta only
  #m=round(n^r);
  N=n*m - m*(m+1)/2;
  wk = seq(n-1,n-m)/N;   
  dk = (seq(1,m)/n)^2;
  dwb=sum(wk*dk);
  skn = vector();
  for(k in 1:m){ skn[k] = sum( (y[(1+k):n] - y[1:(n-k)])^2)}
  sk = skn/(2*(n-seq(1:m)));
  beta.h =sum(wk*sk*(dk - dwb))/sum(wk*(dk-dwb)^2);
  return(beta.h)
}

diff.test = function(n,sg,ns,a,lmd=0,f,method){
  # This is the simulation function of the difference based test.
  # n=sample size, sg=sigma, ns=number of simulation, lmd is to control the bandwidth,
  # a is the coefficient, f is the function number
  k=k2=matrix(0,ns,length(a));
  x<<- 1:n/n;
  for(j in 1:length(a)){
    beta=vector() ; sg2.h = vector();
    for(i in 1:ns){
      if(f==1){y=1+a[j]*5*(x^2-x)+rnorm(n,0,sg)}      # function g_1
      if(f==2){y=1+a[j]*sin(4*pi*x)+rnorm(n,0,sg)}    # function g_2
      if(f==3){y=1+a[j]*5*(x^2-x)+rt(n,df=3)/10*sqrt(3)} # case of non-Gaussian random errors.
      result = diff.fn(n,lmd,y,method);
      if(result$test.Gau > qnorm(0.95) ) k[i,j]=k[i,j]+1   # DBT-Gau
      if(result$test.Gen > qnorm(0.95) ) k2[i,j]=k2[i,j]+1 # DBT-Gen
      beta = c(beta,result$beta.h);sg2.h = c(sg2.h,result$sg2.h)
    }
  }
  DBT.Gau=colMeans(k)     ## Test result using T test statistics
  DBT.Gen=colMeans(k2)    ## Test result using G test statistics
  rate=rbind(DBT.Gau,DBT.Gen)
  return(rate)
}

diff.Perm.test = function(n,sg,ns,a,lmd=0,f,method){
  # This is the simulation function of the difference based test.
  # n=sample size, sg=sigma, ns=number of simulation, lmd is to control the bandwidth,
  # a is the coefficient, f is the function number
  DBT.Perm=NULL; x<<- 1:n/n;
  for(j in 1:length(a)){
    beta=vector() ; sg2.h = vector(); pv <- NULL
    for(i in 1:ns){
      if(f==1){y=1+a[j]*5*(x^2-x)+rnorm(n,0,sg)}      # function g_1
      if(f==2){y=1+a[j]*sin(4*pi*x)+rnorm(n,0,sg)}    # function g_2
      if(f==3){y=1+a[j]*5*(x^2-x)+rt(n,df=3)/10*sqrt(3)} # case of non-Gaussian random errors.
      result = diff.fn(n,lmd,y,method);
      beta = c(beta,result$beta.h);sg2.h = c(sg2.h,result$sg2.h); m<-result$m;
      perm=NULL;
      for (l in 1:200) perm=c(perm, diff.fn.b(n,m,y=sample(y))) # We set the number of permutation as 200
      pv=c(pv, mean(perm>result$beta.h))
    }
    DBT.Perm=c(DBT.Perm,mean(pv<.05)) ## Test result using permutation method
  }
  DBT.Perm <- matrix(DBT.Perm,nrow=1)
  rownames(DBT.Perm)<- "DBT.Perm"
  return(DBT.Perm)
}

a=c(0,0.2,0.5,0.7,1);
## result for the first function g_1
for(n in c(30,50,100) ){
startTime <- Sys.time()
set.seed(1)
print(round(diff.test(n,0.3,1000,a,0,1,2),digits=3))
endTime <- Sys.time()
print(endTime - startTime)

startTime <- Sys.time()
set.seed(1)
print(round(diff.Perm.test(n,0.3,1000,a,0,1,2),digits=3))
endTime <- Sys.time()
print(endTime - startTime)
}

## result for the second function g_2

set.seed(4)
round(diff.test(30,0.5,1000,a,0.2,2,1),digits=3)
set.seed(4)
round(diff.Perm.test(30,0.5,1000,a,0.2,2,1),digits=3)

set.seed(4)
round(diff.test(50,0.5,1000,a,0.05,2,1),digits=3)
set.seed(4)
round(diff.Perm.test(50,0.5,1000,a,0.05,2,1),digits=3)

set.seed(4)
round(diff.test(100,0.5,1000,a,0,2,1),digits=3)
set.seed(4)
round(diff.Perm.test(100,0.5,1000,a,0,2,1),digits=3)

## result for the third function g_3

for(n in c(30,50,100) ){
  set.seed(1)
  print(round(diff.test(n,0.3,1000,a,0,3,2),digits=3))
  set.seed(1)
  print(round(diff.Perm.test(n,0.3,1000,a,0,3,2),digits=3))
}
 
########## F-Permutation method ##########

library(assist)
fstat=function(x,y,n){
  fit=ssr(y~1, linear(x))
  S1=sum((fit$fit-mean(y))**2)
  S3=sum((y-mean(y))**2)
  H=hat.ssr(fit)
  g1=sum(diag(H%*%H))
  return((n-g1)*S1/((g1-1)*(S3-S1)))
}

perm=function(n,sg,nsim,a,f){
  power=NULL
  x<<-1:n/n; nperm=200;
  for(k in 1:length(a)){
    pv <- NULL
    for(i in 1:nsim){
      ep=rnorm(n,0,sg)
      if(f==1){y=1+a[k]*5*(x^2-x)+ep}      # function g_1
      if(f==2){y=1+a[k]*sin(4*pi*x)+ep}    # function g_2
      if(f==3){y=1+a[k]*5*(x^2-x)+rt(n,df=3)/10*sqrt(3)} 
      # case of non-Gaussian random errors.
      f1=fstat(x=x,y=y,n=n)
      f2=NULL
      for (j in 1:nperm) f2=c(f2, fstat(x=x,y=sample(y),n))
      pv=c(pv, mean(f2>f1))
    }
    power=c(power,mean(pv<.05))
  }
  print(power)
}

a=c(0,0.2,0.5,0.7,1)

for(n in c(30,50,100)){
startTime <- Sys.time()
set.seed(1)
perm(n,0.3,1000,a,f=1)
endTime <- Sys.time()
print(endTime - startTime)
}

for(n in c(30,50,100)){
  set.seed(4)
  perm(n,0.5,1000,a,f=2)
}

for(n in c(30,50,100)){
  set.seed(1)
  perm(n,0.3,1000,a,f=3)
}


########## LMP test ##########

install.packages(assist)
library(assist) 

LMP.test <- function(n,sg,a,f,ns=1000,simusize=1000){
  pv <- matrix(NA,ns,length(a)) 
  for (j in 1:length(a)){
    for (i in 1:ns) {
      x<<-1:n/n
      if(f==1){y=1+a[j]*5*(x^2-x)+rnorm(n,0,sg)}
      if(f==2){y=1+a[j]*sin(4*pi*x)+rnorm(n,0,sg)}
      if(f==3){y=1+a[j]*5*(x^2-x)+rt(n,df=3)/10*sqrt(3)}
      fit <- ssr(y~1,linear(x),spar="m")
      tmp <- anova(fit, simu.size=simusize)
      pv[i,j] <- mean(tmp$lmp.test$simu>tmp$lmp.test$value)
    }}
  LMP.result <- apply(pv<.05,2,mean)
  print(LMP.result)
}

a=c(0,0.2,0.5,0.7,1);
for(n in c(30,50,100)){
startTime <- Sys.time()
set.seed(1)
LMP.test(n=n,sg=0.3,a,f=1,ns=1000,simusize=1000)
endTime <- Sys.time()
print(endTime - startTime)}

for(n in c(30,50,100)){
  set.seed(4)
  LMP.test(n=n,sg=0.5,a,f=2,ns=1000,simusize=1000)
}

for(n in c(30,50,100)){
  set.seed(1)
  LMP.test(n=n,sg=0.3,a,f=3,ns=1000,simusize=1000)
}

########## Yatchew's specification method ##########

install.packages("VarED")
require(VarED)
comb = function(n,x){
  factorial(n)/factorial(n-x)/factorial(x)
}

Yat.ord.fn = function(n,sigma,a,r,f){  ## ordinary difference sequence
  x=1:n/n
  if(f==1){y=1+a*5*(x^2-x)+sigma*rnorm(n)}
  if(f==2){y=1+a*sin(4*pi*x)+sigma*rnorm(n)}
  if(f==3){y=1+a*5*(x^2-x)+rt(n,df=3)/10*sqrt(3)}
  s.res = sum((y - mean(y))^2)/(n-1)
  d=vector();D=matrix(0,n,n);d.m = matrix(0,r+1,r+1)
  for(j in 1:(r+1)) d[j] = (-1)^(j-1) *comb(2*r,r)^(-1/2)*comb(r,j-1)    #vector of differencing weights
  for(i in 1:(n-r)) for(j in 1:(r+1)) D[i,j+i-1] = d[j]                 # differencing matrix
  for(i in 1:(r+1)) for(j in 1:(r+1)) if(i+j <= r+1) d.m[i,j] = d[i+j]  # d.m is a matrix used to calculate delta
  delta =sum((t(d) %*%d.m)^2 )
  s.ord = sum((D%*%y)^2)/(n-r)                #ordinary differencing sequence
  test.Y = sqrt(n/(4*delta))*(s.res - s.ord)/s.res   # test statistic
  return(test.Y)
}

Yat.opt.fn = function(n,sigma,a,r,f){  ## ordinary difference sequence
  x=1:n/n
  if(f==1){y=1+a*5*(x^2-x)+sigma*rnorm(n)}
  if(f==2){y=1+a*sin(4*pi*x)+sigma*rnorm(n)}
  if(f==3){y=1+a*5*(x^2-x)+rt(n,df=3)/10*sqrt(3)}
  s.res = sum((y - mean(y))^2)/(n-1)
  D=matrix(0,n,n);d.m = matrix(0,r+1,r+1)
  d=optseq(r)            #vector of differencing weights
  for(i in 1:(n-r)) for(j in 1:(r+1)) D[i,j+i-1] = d[j]                 # differencing matrix
  s.opt = sum((D%*%y)^2)/(n-r)                #ordinary differencing sequence
  test.Y = sqrt(n*r)*(s.res - s.opt)/s.res  # test statistic
  return(test.Y)
}

Yat.test = function(n,sg,ns,r,a,f){  
  #n=sample size, sg=sd of error, r=seq order, a=coeff of mean function, f=# of function
  k=k2=matrix(0,ns,length(a))
  ci = qnorm(0.975)
  for(j in 1:length(a)){
    for(i in 1:ns){
      test.ord = Yat.ord.fn(n,sg,a[j],r,f)
      test.opt = Yat.opt.fn(n,sg,a[j],r,f)
      if(test.ord > ci | test.ord < -ci) k[i,j]=k[i,j]+1
      if(test.opt > ci | test.opt < -ci) k2[i,j]=k2[i,j]+1
    }
  }
  k.c=colSums(k)
  rate.ord = k.c/ns
  k2.c=colSums(k2)
  rate.opt = k2.c/ns
  rate=rbind(rate.ord,rate.opt )
  return(rate)
}

a=c(0,0.2,0.5,0.7,1)
for(n in c(30,50,100)){
  set.seed(1)
  startTime <- Sys.time()
  print(Yat.test(n,0.3,1000,r=2,a,f=1))
  endTime <- Sys.time()
  print(endTime - startTime)
}

for(n in c(30,50,100)){
  set.seed(4)
  print(Yat.test(n,0.5,1000,r=2,a,f=2))
}

for(n in c(30,50,100)){
  set.seed(1)
  print(Yat.test(n,0.3,1000,r=2,a,f=3))
}

### Kolmogorov-Smirnov Type Statistic ###
install.packages("EnvStats")
library(EnvStats)
kernel.quad.n=function(u) { 1/sqrt(2*pi)*exp(-u^2/2) }

kernel.mean=function(y,x,bdw,K=kernel.quad.n,dx=1,dy=1,cv=F,ind.print=F)  
  ## kernel estimator E(Y|X=xi) for all i
  ## x is nxp, y is nxd, x0 is px1, product kernel is used with same kernel function
{
  if(dy==1) n=length(y) else n=nrow(y);  ## y is a matrix if dy>1
  ones=rep(1,n); KD=matrix(1,n,n);  ## initial weight set to 1
  if(dx==1) x=as.matrix(x,ncol=1)
  for(r in 1:dx)## take the r-th column and find the wights K((xir-xjr)/bdw_r)
  {  s=x[,r];   
  s.mtx=s%*%t(ones); D.s=s.mtx-t(s.mtx)  ## matrix with entry (i,j)=si-sj
  KD=KD*K(D.s/bdw[r]);   ## update the weights for the r-th dimension in x
  }
  if(cv) diag(KD)=0  ## if cross validation is used, set diag=0
  Ky=t(KD)%*%y;  ## ith row is \sum_i K((x_i-x_j)/bdw)*yi
  Ks=colSums(KD);  ## Ks is nx1, jth col of KD is for E(X|beta'X=beta'Xj)
  
  if(cv) {index=(Ks==0); Ks[index]=1; Ey.x=Ky/Ks*(1-index);}  ## no other points in the interval
  ##(Ks-K(0)^dx)/K(0)^dx<1e-4 ); 
  ## if Ks=0 or if Ks=K(0)^dx, meaning the row in KD is (0,...,c,...0), 
  ## which makes x.b=x for that 
  else {if(any(Ks==0)) {print("zero denominator"); 
  return("a")} else Ey.x=Ky/Ks; }
  #if(sum(abs(Ey.x))<1e-8) print(cbind(Ky,Ks))
  if(dy==1) Ey.x=as.vector(Ey.x)    ## if dy=1, make it a vector (since y was a vector)
  return(Ey.x);   ## xK/Ks is nxp where division is by column
}

Keil.test <- function(x,y,n){
  h0=1.06*sd(x)*n^(-1/6);
  F.d <- numeric()
  m.hat <- kernel.mean(y,x,h0,K = kernel.quad.n)
  sg2.hat <- mean((y-m.hat)^2)
  e.hat <- (y - m.hat)/sqrt(sg2.hat)
  e.hat0 <- (y - mean(y))/sqrt(sg2.hat)
  e.s <- sort(e.hat)
  for(i in 1:length(e.s)) {
  F.d[i] <- abs(sum(e.hat <= e.s[i])/n - sum(e.hat0 <= e.s[i])/n )} 
  # F.d is the difference of the Cdf
  Tks <- sqrt(n)*max(F.d)
  
  e.sd <-(e.hat - mean(e.hat))/sd(e.hat)
  B = 1000
  D.Tks <- numeric();
  for(b in 1:B){
    e.star <- remp(n,e.sd)
    h1 <- 1.06*sd(e.star)*n^(-1/5);
    f.e <- numeric()
    for(i in 1:n) {
    f.e[i] <- mean(kernel.quad.n((e.star[i] - e.star)/h1))/h1 }
    D.Tks[b] <- max(abs(f.e))* abs(rnorm(1,0,1))
  }
  Tks.t <- Tks > quantile(D.Tks,0.95)
  return(list(Tks.t=Tks.t))
}

K.Test <- function(nv,sg,ns,a,f){
  p.Tks <- matrix(NA,length(nv),length(a))
  for(k in 1:length(nv)){
    n <- nv[k]
    x <- seq(1:n)/n
    for(j in 1:length(a)){
      Tks.t.sim <- numeric(ns);
      for(i in 1:ns){
        if(f==1){y=1+a[j]*5*(x^2-x)+rnorm(n,0,sg)}
        if(f==2){y=1+a[j]*sin(4*pi*x)+rnorm(n,0,sg)}
        if(f==3){y=1+a[j]*5*(x^2-x)+rt(n,df=3)/10*sqrt(3)}
        Keil.t <- Keil.test(x,y,n)
        Tks.t.sim[i] <- Keil.t$Tks.t
      }
      p.Tks[k,j] <- mean(Tks.t.sim);
    }
  }
  return(p.Tks)
}

for(n in c(30,50,100)){
startTime <- Sys.time()
set.seed(1)
print(K.Test(nv=n,sg=0.3,ns=100,a=c(0,0.2,0.5,0.7,1),f=1))
endTime <- Sys.time()
print(endTime - startTime)
}

for(n in c(30,50,100)){
  set.seed(4)
  print(K.Test(nv=n,sg=0.5,ns=1000,a=c(0,0.2,0.5,0.7,1),f=2))
}

for(n in c(30,50,100)){
  set.seed(1)
  print(K.Test(nv=n,sg=0.3,ns=1000,a=c(0,0.2,0.5,0.7,1),f=3))
}