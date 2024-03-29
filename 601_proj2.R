### STAT-601 Spring-2023-Project-Part-2 by SAMIPAN MAJUMDER
### Date: 22-April-2023.

### To DO: Parametric Bootstrap for Confidence Intervals in the case of Question-2


setwd("D:/Dell_Laptop/Desktop/ISU US 19aug_or_ 25Aug2021and09Oct2021 onward/STAT 601 Advanced Statistical Methods S2022/Spring2023/Project-Part-2/")
getwd()

## Load Dataset

s601.proj2.data <- read.table(file = "projdat.txt",header = T,fill = T)
str(s601.proj2.data)
print(s601.proj2.data)
nrow(s601.proj2.data)
plot(x = s601.proj2.data$x,y = s601.proj2.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")

## Loading XYdata - Part-2 Project

XYdata.2 <- cbind(s601.proj2.data$x,s601.proj2.data$y)
#XYdata.2
colnames(XYdata.2) <- c("X","Y")
XYdata.2

## Functions/Methods from Part-1 Project


## Non-Linear Algorithm by Prof. Kaiser in STAT601 course

"gauss.newton"<-
  function(xmat, y, ps, fctn, ders, wts)
  {
    # fctn is a function that computes the current 
    # expectations, ders is a function that computes
    # the current matrix of derivatives of the 
    # expectation function (n by p matrix), wts is a function
    # that computes weights as 1/g^2 where g is
    # the variance model
    # xmat are covariates that can be in any
    # form (vector, matrix, list, etc.) that is
    # expected by fctn, ders, and wts
    cps <- ps
    cnt <- 0
    crit <- 1e-08
    repeat {
      cnt <- cnt + 1
      cps1 <- cps
      ee <- fctn(xmat, cps)
      # print(ee)
      V <- ders(xmat, cps)
      #  print(V)
      W <- wts(xmat, cps)
      # print(W)
      d <- solve(t(V) %*% W %*% V)
      #  cat("d:",d,fill=T)
      d <- d %*% t(V) %*% W %*% as.matrix(y - ee)
      cps <- cps1 + d
      cat("New Estimates at Iteration",cnt,":", fill = T)
      cat(cps, fill = T)
      dist <- (sum((cps - cps1)^2))^0.5
      if(dist < crit) break
    }
    cat("Convergence Criterion of ", crit, " met", fill = T)
    cat("Final Estimates: ", cps, fill = T)
    cps
  }

"nonlin"<-
  function(xmat, ys, ps, fctn, ders, wts)
  {
    # external function definitions for fctn,
    # ders and wts are as in gauss.newton
    # output is list containing betahat(bs),
    # sigma squared hat (sshat), the covariance
    # matrix for betahat (cov), fitted values
    # (yh), studentized residuals (bi), 
    # absolute studentized residuals to
    # (2/3) power (abi),
    # and the matrix of derivatives of the
    # expectation function (fb)
    N <- length(ys)
    P <- length(ps)
    bs <- gauss.newton(xmat, ys, ps, fctn, ders, wts)
    yh <- fctn(xmat, bs)
    r <- ys - yh
    w <- wts(xmat, bs)
    g2 <- 1/(diag(w))
    ri <- r/(g2^0.5)
    fb <- ders(xmat, bs)
    G <- matrix(sqrt(g2), N, P, byrow = F)
    xx <- fb/G
    H <- xx %*% (solve(t(xx) %*% xx)) %*% t(xx)
    h <- diag(H)
    sshat <- (1/(N - P)) * sum(ri^2)
    cov <- matrix(0, P, P)
    cnt <- 0
    repeat {
      cnt <- cnt + 1
      tfb <- as.matrix(fb[cnt,  ])
      tfbfb <- tfb %*% t(tfb)
      tel <- tfbfb/g2[cnt]
      cov <- cov + tel
      if(cnt == N)
        break
    }
    cov <- sshat * solve(cov)
    bi <- ri/((sshat * (1 - h))^0.5)
    abi <- (abs(bi))^(2/3)
    result <- list(bs=bs, sshat=sshat, covb=cov, yhat=yh, stdres=bi, absres=abi, fb=fb)
    result
  }

mu.i.func <- function(X,par){
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  mu.i <- alpha*X^(beta*X^(-gamma))
  return  (mu.i) 
}

der.func <- function(X,par){
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  mu.i <- alpha*X^(beta*(X^(-gamma)))
  dmu_i.dalpha <- X^(beta*X^(-gamma))
  dmu_i.dbeta <-  mu.i*(X^(-gamma))*log(X,base = exp(1))
  #dmu_i.dgamma <- mu.i*log(mu.i/alpha,base=exp(1))*log(X^(-1),base=exp(1))
  #dmu_i.dgamma <- mu.i*log(X,base=exp(1))*log(X^(-1),base=exp(1))*X^(-gamma)
  dmu_i.dgamma <- beta*mu.i*log(X,base=exp(1))*log(X^(-1),base=exp(1))*X^(-gamma)
  
  V.j <- matrix(data = c(dmu_i.dalpha,dmu_i.dbeta,dmu_i.dgamma),nrow =length(X),ncol = length(par))
  return (V.j)
}

wts.func <- function(X,par){
  return (diag(length(X)))
}



# GLS algorithm from STAT520

Q4.GLS.1 <- function(par, XYdata,iter=1000){
  cnt <- 0
  X <- XYdata[,1]
  Y <- XYdata[,2]
  repeat{
    cnt <- cnt+1   
    X <- XYdata[,1]
    Y <- XYdata[,2]
    n <- nrow(XYdata)
    alpha <- par[1]
    beta <- par[2]
    gamma <- par[3]
    mu.i <- alpha*X^(beta*X^(-gamma))
    g2.phi <- 1
    W <- diag(n)
    dmu_i.dalpha <- X^(beta*X^(-gamma))
    dmu_i.dbeta <-  mu.i*(X^(-gamma))*log(X,base = exp(1))
    # dmu_i.dgamma <- mu.i*log(mu.i/alpha,base=exp(1))*log(X^(-1),base=exp(1))
    dmu_i.dgamma <- beta*mu.i*log(X,base=exp(1))*log(X^(-1),base=exp(1))*X^(-gamma)
    V.j <- matrix(data = c(dmu_i.dalpha,dmu_i.dbeta,dmu_i.dgamma),nrow =n,ncol = 3)
    Y.j <- Y-mu.i
    #W.mat <- diag(W)
    del <- solve(t(V.j) %*% W %*% V.j) %*% t(V.j) %*% W %*% as.matrix(Y.j)
    par.old <- par
    #cat(par,"\n")
    #cat(par.old,"\n")
    #par <- par+diag(del)
    par <- par+(del)
    #cat(par,"\n")
    #if (all(abs(par-par.old)<1e-10)){
    if (((sum((par - par.old)^2))^0.5)<1e-10){
      break
    } else if (cnt>iter){
      #cat("Convergence Not Obtained Within ",cnt," Steps","\n")
      break
    }
    cat("Estimates at Step",cnt,par,fill=T)
  }
  
  alpha.est <- par[1]
  beta.est <- par[2]
  gamma.est <- par[3]
  p <- 3 # Number of parameters in the Systematic Component
  n <- nrow(XYdata)
  mu.i.est <- alpha.est*X^(beta.est*X^(-gamma.est))
  g2.phi.est <- 1
  sigma.sq.1 <- (1/(n-p))*sum(((Y-mu.i.est)/g2.phi.est)^2)
  
  W.est <- diag(n)
  dmu_i.dalpha.est <- X^(beta.est*X^(-gamma.est))
  dmu_i.dbeta.est <-  mu.i.est*(X^(-gamma.est))*log(X,base = exp(1))
  #dmu_i.dgamma.est <- mu.i.est*log(mu.i.est/alpha.est,base=exp(1))*log(X^(-1),base=exp(1))
  dmu_i.dgamma.est <- beta.est*mu.i.est*log(X,base=exp(1))*log(X^(-1),base=exp(1))*X^(-gamma.est)
  V.j.est <- matrix(data = c(dmu_i.dalpha,dmu_i.dbeta,dmu_i.dgamma),nrow = n,ncol = 3)
  
  sigma.mat <- solve(t(V.j.est)%*%W.est%*%V.j.est)
  # sigma.mat <- (1/n)*(t(V.j.est)%*%V.j.est)*W.est
  
  
  # vcov.vector <- matrix(data = 0,nrow = 3,ncol = 3,byrow = F)
  # vcov.sum <- vcov.vector
  # for (i in 1:n){
  #   vcov.vector <- (V.j.est[i,])%*%(t(V.j.est[i,]))*(W.est[i,i])
  #   vcov.sum = vcov.sum + vcov.vector
  # }
  # #sigma.mat.1 <- sum(vcov.vector)/n
  # sigma.mat.1 <- vcov.sum*(1/n)
  # sigma.mat.1 <- solve(as.matrix(sigma.mat.1))
  
  
  # return (par)
  #return (list(pars=par,sigma.2=sigma.sq.1,sigma.Matrix=sigma.mat,vcov.Mat=sigma.mat.1))
  return (list(pars=par,sigma.2=sigma.sq.1,sigma.Matrix=sigma.mat,vcov.Mat=sigma.sq.1*sigma.mat,Derivatives.Mat=V.j.est))
}


# Bootstrap Sample Function

boot.sample <- function(M=2500,alpha=0.05,theta,var,n,type){
  theta.star <- c()
  sigma.sq.star <- c()  
  
  for (i in 1:M){
    w.i.boot <- c()
    e.i.boot <- c()
    w.0 <- rnorm(n = 1,mean = 0,sd = sqrt(var/(1-theta^2)))
    e.i.boot[1] <- rnorm(n = 1,mean = 0,sd = sqrt(var))
    w.i.boot[1] <- theta*w.0 + e.i.boot[1]
    for (j in 2:n){
      e.i.boot[j] <- rnorm(n = 1,mean = 0,sd = sqrt(var))
      w.i.boot[j] <- theta*w.i.boot[j-1] + e.i.boot[j]
    }
    
    if (type=="naive"){
      theta.star[i] <- sum(w.i.boot)/sum(w.i.boot[1:(n-1)])
      sigma.sq.star[i] <- (sum((w.i.boot[2:n]-theta.star[i]*w.i.boot[1:(n-1)])^2))/n
    } else if (type=="optimal"){
      theta.star[i] <- sum(w.i.boot[2:n]*w.i.boot[1:(n-1)])/sum(w.i.boot[1:(n-1)]^2)
      sigma.sq.star[i] <- (sum((w.i.boot[2:n]-theta.star[i]*w.i.boot[1:(n-1)])^2))/n
    }
    
    
  }
  bias.theta <- sum(theta.star-theta)/M
  bias.theta
  
  bias.sigma.sq <- sum(sigma.sq.star-var)/M
  bias.sigma.sq
  
  theta.percentile.lower <- quantile(x = theta.star,probs = alpha/2)
  theta.percentile.upper <- quantile(x = theta.star,probs = 1-alpha/2)
  
  sigma.sq.percentile.lower <- quantile(x = sigma.sq.star,probs = alpha/2)
  sigma.sq.percentile.upper <- quantile(x = sigma.sq.star,probs = 1-alpha/2)
  
  theta.basic.lower <- 2*theta-quantile(x = theta.star,probs = 1-alpha/2)
  theta.basic.upper <- 2*theta-quantile(x = theta.star,probs = alpha/2)
  
  sigma.sq.basic.lower <- 2*var-quantile(x = sigma.sq.star,probs = 1-alpha/2)
  sigma.sq.basic.upper <- 2*var-quantile(x = sigma.sq.star,probs = alpha/2)
  
  #Result <- cbind(bias=as.numeric(c(bias.theta,bias.sigma.sq)),Basic.Interval=as.numeric(c(c(theta.basic.lower,theta.basic.upper),c(sigma.sq.basic.lower,sigma.sq.basic.upper))),Percentile.Interval=as.numeric(c(c(theta.percentile.lower,theta.percentile.upper),c(sigma.sq.percentile.lower,sigma.sq.percentile.upper))))
  # row.names(Result) <- c("Theta","Sigma.sq")
  
  Result1 <- cbind(c(bias=as.numeric(bias.theta),basic.interval=as.numeric(c(theta.basic.lower,theta.basic.upper)),Percentile.Interval=as.numeric(c(theta.percentile.lower,theta.percentile.upper))))
  colnames(Result1) <- "theta"
  Result2 <- cbind(c(bias=as.numeric(bias.sigma.sq),basic.interval=as.numeric(c(sigma.sq.basic.lower,sigma.sq.basic.upper)),Percentile.Interval=as.numeric(c(sigma.sq.percentile.lower,sigma.sq.percentile.upper))))
  colnames(Result2) <- "Sigma.sq"
  
  Result3 <- cbind(Result1,Result2)
  
  return (list(theta = theta.star,var=sigma.sq.star,type=type,result=Result3,res.rounded=round(Result3,4)))
}


# Initial values : Case-4 (From Part-1 Project) Initial values for alpha, beta, gamma.

alpha.init <- seq(from=0,to=2,length.out=10)
beta.init <- seq(from=0,to=2,length.out=10)
gamma.init <- seq(from=0,to=2,length.out=10)
length(alpha.init)
length(beta.init)
length(gamma.init)

initvals <- list()
cnt <- 0
for (i in 1:length(alpha.init)){
  for (j in 1:length(beta.init)){
    for (k in 1:length(gamma.init)){
      pars.4 <- c(alpha.init[i],beta.init[j],gamma.init[k])
      tryCatch(
        expr = {
          mod.GLS <- Q4.GLS.1(par = pars.4,XYdata = XYdata.2)
          mod.nonlin <- nonlin(xmat=XYdata.2[,1], ys=XYdata.2[,2], ps=pars.4, fctn=mu.i.func, ders=der.func, wts=wts.func)
          message("\nGLS and Non-Linear functions have completed successfully\n")
          cat("Init values Correct:","alpha = ",pars.4[1],"beta= ",pars.4[2],"gamma=",pars.4[3],fill=T,"\n")
          cnt <- cnt+1
          initvals[[cnt]] <- c(pars.4[1],pars.4[2],pars.4[3])
        }, error = function(e){
          message(e)
          cat("\nInit values wrong:","alpha = ",pars.4[1],"beta= ",pars.4[2],"gamma=",pars.4[3],fill=T,"\n")
        }, warning = function(w){
          message(w)
          cat("\nInit values warning:","alpha = ",pars.4[1],"beta= ",pars.4[2],"gamma=",pars.4[3],fill=T,"\n")
        },finally = {
          message("\nDONE with try catch\n")
          #cat("Init values Correct:","alpha = ",pars.4[1],"beta= ",pars.4[2],"gamma=",pars.4[3],fill=T,"\n")
          #cnt <- cnt+1
          #initvals[[cnt]] <- c(pars.4[1],pars.4[2],pars.4[3])
        }
      )
    }
  }
}
initvals
length(initvals)
initvals[[1]]
# alpha= 0.2222222,beta= 2.0000000,gamma= 0.2222222 Initial values.


### Q.1. Maximum Likelihood Estimation (MLE)


# 1.a) ## log-likelihood MLE

"Pro2.Lhood.MLE"<- function(X,Y,pars){
  alpha <- pars[1]
  beta <- pars[2]
  gamma <- pars[3]
  theta <- pars[4]
  sigma.sq <- pars[5]
  n <- length(Y)
  mu.i <- alpha*X^(beta*(X^((-1)*gamma)))
  
  # cat("Mu.i",mu.i,fill=T)
  # cat("Y[-1]",Y[-1],fill=T)
  # cat("Mean of Y_i 2 to n",mu.i[-1]+theta*(Y[-n]-mu.i[-n]),fill=T)
  # cat("Y[1]",Y[1],fill=T)
  # cat("mu.i.1",mu.i[1],fill=T)
  # cat("sigma.sq",sigma.sq,"sd",sqrt(sigma.sq),fill=T)
  # cat("dnorm.1",dnorm(x = Y[1],mean = mu.i[1],sd = sqrt(sigma.sq)),fill=T)
  # cat("dnorm.2",dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq)),fill=T)
  # cat("prod.dnorm.2",prod(dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq))),fill=T)
  
  Lhood <- dnorm(x = Y[1],mean = mu.i[1],sd = sqrt(sigma.sq))*prod(dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq)))
  logLhood <- log(x = Lhood,base = exp(1))
  
  logLhood.1 <- ((-1)*(n/2))*log(2*pi*sigma.sq,base = exp(1))-((Y[1]-mu.i[1])^2/(2*sigma.sq))-sum((Y[-1]-(mu.i[-1]+theta*(Y[-n]-mu.i[-n])))^2/(2*sigma.sq))
  
  logLhood.2 <- log(dnorm(x = Y[1],mean = mu.i[1],sd = sqrt(sigma.sq))) + sum(log(dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq))))
  
  sum.2.n <- 0
  for (i in 2:n){
    sum.2.n <- sum.2.n + (Y[i]-(mu.i[i]+theta*(Y[i-1]-mu.i[i-1])))^2
  }
  
  logLhood.3 <- (-n/2)*log(2*pi*sigma.sq)-(Y[1]-mu.i[1])^2/(2*sigma.sq) - sum.2.n/(2*sigma.sq)
  
  return (list(Lhood=Lhood,logLhood=logLhood,logLhood.1=logLhood.1,logLhood.2=logLhood.2,logLhood.3=logLhood.3))
}

"Pro2.negLhood.MLE"<-function(X,Y,pars){
  
  #neg.Lhood <- (-1)*Pro2.Lhood.MLE(X,Y,pars)$logLhood
  neg.Lhood <- (-1)*Pro2.Lhood.MLE(X,Y,pars)$logLhood.1
  #neg.Lhood <- (-1)*Pro2.Lhood.MLE(X,Y,pars)$logLhood.2
  #neg.Lhood <- (-1)*Pro2.Lhood.MLE(X,Y,pars)$logLhood.3
  
  return (neg.Lhood)
}





# 1.b.) ## MLE Estimation Procedure

## MLE Initial values - Case-1  *** # alpha, beta, gamma (First set from the list of initial values from Part-1), theta and sigma.sqaure (Randomly Chosen).

pro2.pars.init.1 <- c(0.2222222, 2.0000000, 0.2222222,0,0.5)  # alpha, beta, gamma, theta, sigma.sq 
pro2.pars.init.1

optim(par = pro2.pars.init.1,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead") # Error: function cannot be evaluated at initial parameters

optim(par = pro2.pars.init.1,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS") # Error: initial value in 'vmmin' is not finite

nlm(f = Pro2.negLhood.MLE,p = pro2.pars.init.1,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2]) # Initial Value is final Estimate. Hessian is ZERO Matrix, Gradient is zERO.

Pro2.Lhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.pars.init.1)

Pro2.negLhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.pars.init.1)


## MLE Initial values - Case-2 *** # alpha, beta, gamma (First set from the list of initial values from Part-1), theta and sigma.sqaure (Randomly Chosen).

pro2.pars.init.2 <- c(0.2222222, 2.0000000, 0.2222222,0.1,0.5)  # alpha, beta, gamma, theta, sigma.sq 
pro2.pars.init.2

optim(par = pro2.pars.init.2,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead") # Error: function cannot be evaluated at initial parameters

optim(par = pro2.pars.init.2,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

nlm(f = Pro2.negLhood.MLE,p = pro2.pars.init.2,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2]) # Initial Value is final Estimate. Hessian is ZERO Matrix, Gradient is zERO.


## MLE Initial values - Case-3 *** # alpha, beta, gamma (First set from the list of initial values from Part-1), theta and sigma.sqaure (Randomly Chosen).

pro2.pars.init.3 <- c(0.22, 2.00, 0.22,0.1,0.5)  # alpha, beta, gamma, theta, sigma.sq 
pro2.pars.init.3

optim(par = pro2.pars.init.3,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead") # Error: function cannot be evaluated at initial parameters

optim(par = pro2.pars.init.3,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS") # Error: initial value in 'vmmin' is not finite

nlm(f = Pro2.negLhood.MLE,p = pro2.pars.init.3,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2]) # Initial Value is final Estimate. Hessian is ZERO Matrix, Gradient is zERO.

Pro2.Lhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.pars.init.3)

Pro2.negLhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.pars.init.3)


## MLE Initial Values : Case-4



initvals
initvals[[1]]
# alpha= 0.2222222,beta= 2.0000000,gamma= 0.2222222 Initial values.

theta.init.1 <- seq(from=-0.9999,to=0.9999,length.out=10)
theta.init.1
length(theta.init.1)
theta.init.1 <- c(theta.init.1,0)
theta.init.1 <- sort(theta.init.1)
theta.init.1
length(theta.init.1)


sigma.sq.init.1 <- seq(from = 0.0001, to = 2, length.out = 11)
sigma.sq.init.1

pro2.pars.init.4 <- c() # alpha, beta, gamma
pro2.pars.init.4
#pro2.pars.init.4[4]

df0.1 <- data.frame()
df0.2 <- data.frame()
df0.3 <- data.frame()

"pro2.case4" <- function(initvals,theta.init.1,sigma.sq.init.1,data=XYdata.2){
  df0.11 <- data.frame()
  df0.21 <- data.frame()
  df0.31 <- data.frame()
for (k1 in 1:length(initvals)){
  pro2.pars.init.4 <- c()
  pro2.pars.init.4 <- c(initvals[[k1]],NA,NA)
for (i1 in 1:length(sigma.sq.init.1)){
  pro2.pars.init.4[5] <- sigma.sq.init.1[i1]
  for (j1 in 1:length(theta.init.1)){
    pro2.pars.init.4[4] <- theta.init.1[j1]
    cat("Initial values: ",pro2.pars.init.4,fill=T,"\n")
    #print(Pro2.Lhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.pars.init.4))
 # print(nlm(f = Pro2.negLhood.MLE,p = pro2.pars.init.4,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])$estimate)
  tryCatch(expr = {
    mod.Lhood.0 <- Pro2.Lhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.pars.init.4)
    mod.MLE.0.1 <- optim(par = pro2.pars.init.4,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")
    mod.MLE.0.2 <- optim(par = pro2.pars.init.4,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")
    mod.MLE.0.3 <- nlm(f = Pro2.negLhood.MLE,p = pro2.pars.init.4,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
    res1 <- c(k1,pro2.pars.init.4[4],pro2.pars.init.4[5],mod.Lhood.0$logLhood.1,mod.Lhood.0$logLhood.2,mod.MLE.0.1$par)
    df0.11 <- rbind(df0.11,res1)
    
    res2 <- c(k1,pro2.pars.init.4[4],pro2.pars.init.4[5],mod.Lhood.0$logLhood.1,mod.Lhood.0$logLhood.2,mod.MLE.0.2$par)
    df0.21 <- rbind(df0.21,res2)
    
    res3 <- c(k1,pro2.pars.init.4[4],pro2.pars.init.4[5],mod.Lhood.0$logLhood.1,mod.Lhood.0$logLhood.2,mod.MLE.0.3$estimate)
    df0.31 <- rbind(df0.31,res3)
  },error=function(e){
    message(e)
  },warning=function(w){
    message(w)
  },finally={
    message("\nDONE with try catch\n")
  })
  

  }
}
}
  return (list(df1=df0.11,df2=df0.21,df3=df0.31))
}

#pro2.case4.init <- pro2.case4(initvals = initvals,theta.init.1 = theta.init.1,sigma.sq.init.1 = sigma.sq.init.1,data = XYdata.2)

if (exists("pro2.case4.init")){
df0.1 <- pro2.case4.init$df1
df0.2 <- pro2.case4.init$df2
df0.3 <- pro2.case4.init$df3
}

if (nrow(df0.1)>0){
colnames(df0.1) <- c("Initval.Index","theta.init","var.init","logLhood1","logLhood2","alpha","beta","gamma","theta","var")}
df0.1
if (nrow(df0.2)>0){
colnames(df0.2) <- c("Initval.Index","theta.init","var.init","logLhood1","logLhood2","alpha","beta","gamma","theta","var")}
df0.2
if (nrow(df0.3)>0){
colnames(df0.3) <- c("Initval.Index","theta.init","var.init","logLhood1","logLhood2","alpha","beta","gamma","theta","var")}
df0.3



## MLE Initial Values : Case-5


pro2.alpha.init <- seq(from=0,to=2,length.out=6)
pro2.beta.init <- seq(from=0,to=2,length.out=6)
pro2.gamma.init <- seq(from=0,to=2,length.out=6)
pro2.theta.init <- seq(from=-0.99,to=0.99,length.out=6)
# pro2.theta.init <- c(0)
#  if(0 %in% seq(from=-0.99,to=0.99,length.out=5),seq(from=-0.99,to=0.99,length.out=6),sort(c(seq(from=-0.99,to=0.99,length.out=5),0)))
pro2.var.init <- seq(from=0.001,to=2,length.out=6)
length(pro2.alpha.init)
length(pro2.beta.init)
length(pro2.gamma.init)
pro2.alpha.init
pro2.beta.init
pro2.gamma.init
pro2.theta.init
pro2.var.init
length(pro2.var.init)
length(pro2.theta.init)

pro2.intval <- list()
"pro2.brute.force" <- function(){
cnt1 <- 0
stop1 = F
for (i2 in 1:length(pro2.alpha.init)){
  for (j2 in 1:length(pro2.beta.init)){
    for (k2 in 1:length(pro2.gamma.init)){
      for (l2 in 1:length(pro2.theta.init)){
        for (m2 in 1:length(pro2.var.init)){
      pars.5 <- c(pro2.alpha.init[i2],pro2.beta.init[j2],pro2.gamma.init[k2],pro2.theta.init[l2],pro2.var.init[m2])
      tryCatch(
        expr = {
          mod.MLE.1 <- optim(par = pars.5,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")
          mod.MLE.2 <- optim(par = pars.5,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")
          mod.MLE.3 <- nlm(f = Pro2.negLhood.MLE,p = pars.5,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
          
          message("\nGLS and Non-Linear functions have completed successfully\n")
          cat("Init values Correct:","alpha = ",pars.5[1],"beta= ",pars.5[2],"gamma=",pars.5[3],"theta = ",pars.5[4],"var= ",pars.5[5],fill=T,"\n")
          cnt1 <- cnt1+1
          pro2.intval[[cnt1]] <- c(pars.5[1],pars.5[2],pars.5[3],pars.5[4],pars.5[5])
          if (cnt1==1000){
            stop1=T
            break
          }
        }, error = function(e){
          message(e)
          cat("\nInit values wrong:","alpha = ",pars.5[1],"beta= ",pars.5[2],"gamma=",pars.5[3],"theta = ",pars.5[4],"var= ",pars.5[5],fill=T,"\n")
        }, warning = function(w){
          message(w)
          cat("\nInit values warning:","alpha = ",pars.5[1],"beta= ",pars.5[2],"gamma=",pars.5[3],"theta = ",pars.5[4],"var= ",pars.5[5],fill=T,"\n")
        },finally = {
          message("\nDONE with try catch\n")
          #cat("Init values Correct:","alpha = ",pars.4[1],"beta= ",pars.4[2],"gamma=",pars.4[3],fill=T,"\n")
          #cnt <- cnt+1
          #initvals[[cnt]] <- c(pars.4[1],pars.4[2],pars.4[3])
        }
      )
        }
        if(stop1==T){break}
      } 
      if(stop1==T){break}
    } 
    if(stop1==T){break}
  } 
  if(stop1==T){break}
}
return (pro2.intval)
}

#pro2.initvals <- pro2.brute.force()

if (exists(x = "pro2.initvals")){
pro2.initvals
length(pro2.initvals)


if (length(pro2.initvals)>0){

  #pro2.initvals[[1]]
#pro2.initvals[[which()]]

# optim(par = pro2.initvals[[1]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")
# optim(par = pro2.initvals[[1]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")
# nlm(f = Pro2.negLhood.MLE,p = pro2.initvals[[1]],hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])

# Pro2.Lhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.initvals[[1]])
# Pro2.Lhood.MLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.initvals[[100]])


df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()

for (i3 in 1:length(pro2.initvals)){
  
  cat("Initial values:","alpha = ",pro2.initvals[[i3]][1],"beta= ",pro2.initvals[[i3]][2],"gamma=",pro2.initvals[[i3]][3],"theta = ",pro2.initvals[[i3]][4],"var= ",pro2.initvals[[i3]][5],fill=T,"\n")
  
  mod1.res <- optim(par = pro2.initvals[[i3]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")
  if (mod1.res$convergence==0){
  res1 <- c(i3,mod1.res$par)
  df1 <- rbind(df1,res1)
  }
  
  mod2.res2 <- optim(par = pro2.initvals[[i3]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")
  if (mod2.res2$convergence==0){
    res2 <- c(i3,mod2.res2$par)
    df2 <- rbind(df2,res2)
  }
  
  
  mod3.res3 <- nlm(f = Pro2.negLhood.MLE,p = pro2.initvals[[i3]],hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
  if (mod3.res3$code==1){
    res3 <- c(i3,mod3.res3$estimate)
    df3 <- rbind(df3,res3)
  }
  
  
}

#df1
colnames(df1) <- c("Init.index","alpha","beta","gamma","theta","var")
df1
colnames(df2) <- c("Init.index","alpha","beta","gamma","theta","var")
df2
colnames(df3) <- c("Init.index","alpha","beta","gamma","theta","var")
df3

df1[df1$theta<1,]
df2[df2$theta<1,]
df3[df3$theta<1,]

nlm(f = Pro2.negLhood.MLE,p = pro2.initvals[[12]],hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
optim(par = pro2.initvals[[12]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS") #absurd values
optim(par = pro2.initvals[[12]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")

optim(par = pro2.initvals[[50]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")
nlm(f = Pro2.negLhood.MLE,p = pro2.initvals[[50]],hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
optim(par = pro2.initvals[[50]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")


pro2.initvals[[12]] #df3 nlm
pro2.initvals[[50]] #df2 bfgs optim
pro2.initvals[[51]] #df1 nelder optim and below
pro2.initvals[[52]] 
pro2.initvals[[106]]
pro2.initvals[[140]] # 0.8000 1.2000 1.6000 0.9900 0.8006 # estimates close to GLs.
pro2.initvals[[141]]
pro2.initvals[[158]]

optim(par = pro2.initvals[[140]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")
optim(par = pro2.initvals[[140]],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS") # Absurd Values
nlm(f = Pro2.negLhood.MLE,p = pro2.initvals[[140]],hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2]) # Absurd Values

}
}


## MLE Initial Values : Case-6


### Final Estimates from Part-1

#theta.naive ˆtheta.optimal
# 0.0589        0.5422

# var.GLS  var.naive  var.optimal
# 1.0149     0.8378      0.629

# alpha   beta   gamma
# 0.5063 1.2220 0.1614

##### Case:6.1
nlm(f = Pro2.negLhood.MLE,p = c(0.5063, 1.2220, 0.1614,0.0589,1.0149),hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])

optim(par = c(0.5063, 1.2220, 0.1614,0.0589,1.0149),fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

optim(par = c(0.5063, 1.2220, 0.1614,0.0589,1.0149),fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")

##### Case:6.2
nlm(f = Pro2.negLhood.MLE,p = c(0.5063, 1.2220, 0.1614,0.5422,0.629),hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])

optim(par = c(0.5063, 1.2220, 0.1614,0.5422,0.629),fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

optim(par = c(0.5063, 1.2220, 0.1614,0.5422,0.629),fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")


###### Case: 6.3  # Same as Case # 6.2
pars <- c(0.5063, 1.2220, 0.1614,0.5422,0.629)
optim(par = pars,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")
optim(par = pars,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")
nlm(f = Pro2.negLhood.MLE,p = pars,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])


##### Case:6.4
nlm(f = Pro2.negLhood.MLE,p = c(0.5063, 1.2220, 0.1614,0.0589,0.8378),hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])

optim(par = c(0.5063, 1.2220, 0.1614,0.0589,0.8378),fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

optim(par = c(0.5063, 1.2220, 0.1614,0.0589,0.8378),fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")


#### Case-7 # Same Initial values for alpha, beta, gamma, only thing varies are theta and variance


# Initial values: alpha= 0.2222222,beta= 2.0000000,gamma= 0.2222222 Initial values, from Part-1 GLS .

theta.init.2 <- seq(from=0,to=0.99,by=0.01)
theta.init.2
#length(theta.init.1)
#theta.init.1 <- c(theta.init.1,0)
#theta.init.1 <- sort(theta.init.1)
#theta.init.1
length(theta.init.2)


sigma.sq.init.2 <- seq(from=0,to=0.99,by=0.01)
sigma.sq.init.2
length(sigma.sq.init.2)

pro2.pars.init.5 <- c(0.2222222,2.0000000,0.2222222,NA,NA) # alpha, beta, gamma
pro2.pars.init.5

pro2.pars.init.6 <- c(0.22,2.00,0.22,NA,NA) # alpha, beta, gamma
pro2.pars.init.6

df7.1.init <- data.frame()
df7.1.final.NM <- data.frame()
df7.1.final.BFGS <- data.frame()
df7.1.final.nlm <- data.frame()
#df7.2 <- data.frame()

for (i4 in 1:length(theta.init.2)){
  pro2.pars.init.5[4] <- theta.init.2[i4]
  for (j4 in 1:length(sigma.sq.init.2)){
    pro2.pars.init.5[5] <- sigma.sq.init.2[j4]
    
    tryCatch(expr={
      mod1.res7.1 <- optim(par = pro2.pars.init.5,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")
      mod2.res7.2 <- optim(par = pro2.pars.init.5,fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")
      mod3.res7.3 <- nlm(f = Pro2.negLhood.MLE,p = pro2.pars.init.5,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
      
      if (mod1.res7.1$convergence==0 && mod2.res7.2$convergence==0 && mod3.res7.3$code==1){
        df7.1.init <- rbind(df7.1.init,pro2.pars.init.5)
        df7.1.final.NM <- rbind(df7.1.final.NM,mod1.res7.1$par)
        df7.1.final.BFGS <- rbind(df7.1.final.BFGS, mod2.res7.2$par)
        df7.1.final.nlm <- rbind(df7.1.final.nlm, mod3.res7.3$estimate)
      }
    },error=function(e){
      message(e)
    },warning=function(w){
      message(w)
    },finally={
      message("\nDONE with try catch\n")
    })
    
    
    
    
    
    
  }
}

df7.1.init
colnames(df7.1.init) <- c("alpha.init","beta.init","gamma.init","theta.init","var.init")
df7.1.init
df7.1.final.BFGS
colnames(df7.1.final.BFGS) <- c("alpha","beta","gamma","theta","var")
df7.1.final.BFGS
df7.1.final.nlm
colnames(df7.1.final.nlm) <- c("alpha","beta","gamma","theta","var")
df7.1.final.nlm
df7.1.final.NM
colnames(df7.1.final.NM) <- c("alpha","beta","gamma","theta","var")
df7.1.final.NM

optim(par = df7.1.init[1,],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")

optim(par = df7.1.init[1,],fn = Pro2.negLhood.MLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

nlm(f = Pro2.negLhood.MLE,p = as.numeric(df7.1.init[1,]),hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])


### Final Full MLE Estimates of Part-2 Q.1.

# We are choosing Full MLE estimates provided by nlm() function, because of the consistency of the Maximized Log-likelihood values, obtained through various combinations of initial values.

# Possible Candidates : 1

# (0.2222222,	2,	0.2222222,	0.99,	0.61) # Initial Values

# (0.7460794,	1.0215273,0.1568744,0.5597127,0.6249797) # Final MLE estimates. # nlm()

# Maximized Log-likelihood value: 29.59801


# Possible Candidates : 2

# (0.2222222,	2.0000000,	0.2222222,	0.1,	0.5) # Initial Values

# (0.7467031,1.0209172,0.1568288,0.5597117,0.6249792) # Final MLE estimates. # nlm()

# Maximized Log-likelihood value: 29.59801





# 1.c) ## MLE Estimates, 95% Wald Intervals, Maximized log-likelihood value.

### Case-1 : Question-1


pro2.final.init.q1.1 <- c(0.2222222,	2,	0.2222222,	0.99,	0.61) # Initial values
pro2.final.init.q1.1
round(pro2.final.init.q1.1,4)


pro2.final.mod.q1.1 <- nlm(f = Pro2.negLhood.MLE,p = pro2.final.init.q1.1,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
pro2.final.mod.q1.1
pro2.final.est.q1.1 <- pro2.final.mod.q1.1$estimate
pro2.final.est.q1.1
pro2.final.hessian.q1.1 <- pro2.final.mod.q1.1$hessian
pro2.final.vcov.q1.1 <- solve(pro2.final.hessian.q1.1)
pro2.final.vcov.q1.1
diag(pro2.final.vcov.q1.1)
pro2.final.se.q1.1 <- sqrt(diag(pro2.final.vcov.q1.1))
pro2.final.se.q1.1
pro2.final.maxLL.q1.1 <- pro2.final.mod.q1.1$minimum
pro2.final.maxLL.q1.1 # MAximum Log-likelihood value
round(pro2.final.maxLL.q1.1,4) # MAx LL
pro2.final.gradient.q1.1 <- pro2.final.mod.q1.1$gradient
pro2.final.gradient.q1.1 # Gradients
round(pro2.final.gradient.q1.1,4)

## For alpha # Q1- case-1

alpha.est.q1.1 <- pro2.final.est.q1.1[1]
alpha.LowerCI.q1.1 <- alpha.est.q1.1 - 1.96*pro2.final.se.q1.1[1]
alpha.UpperCI.q1.1 <- alpha.est.q1.1 + 1.96*pro2.final.se.q1.1[1]
c(alpha.est.q1.1,alpha.LowerCI.q1.1,alpha.UpperCI.q1.1,alpha.UpperCI.q1.1-alpha.LowerCI.q1.1)

## For beta # Q1- case-1

beta.est.q1.1 <- pro2.final.est.q1.1[2]
beta.LowerCI.q1.1 <- beta.est.q1.1 - 1.96*pro2.final.se.q1.1[2]
beta.UpperCI.q1.1 <- beta.est.q1.1 + 1.96*pro2.final.se.q1.1[2]
c(beta.est.q1.1,beta.LowerCI.q1.1,beta.UpperCI.q1.1,beta.UpperCI.q1.1-beta.LowerCI.q1.1)


## For gamma # Q1- case-1

gamma.est.q1.1 <- pro2.final.est.q1.1[3]
gamma.LowerCI.q1.1 <- gamma.est.q1.1 - 1.96*pro2.final.se.q1.1[3]
gamma.UpperCI.q1.1 <- gamma.est.q1.1 + 1.96*pro2.final.se.q1.1[3]
c(gamma.est.q1.1,gamma.LowerCI.q1.1,gamma.UpperCI.q1.1,gamma.UpperCI.q1.1-gamma.LowerCI.q1.1)


## For theta # # Q1- case-1

theta.est.q1.1 <- pro2.final.est.q1.1[4]
theta.LowerCI.q1.1 <- theta.est.q1.1 - 1.96*pro2.final.se.q1.1[4]
theta.UpperCI.q1.1 <- theta.est.q1.1 + 1.96*pro2.final.se.q1.1[4]
c(theta.est.q1.1,theta.LowerCI.q1.1,theta.UpperCI.q1.1,theta.UpperCI.q1.1-theta.LowerCI.q1.1)


## For variance ## # Q1- case-1

var.est.q1.1 <- pro2.final.est.q1.1[5]
var.LowerCI.q1.1 <- var.est.q1.1 - 1.96*pro2.final.se.q1.1[5]
var.UpperCI.q1.1 <- var.est.q1.1 + 1.96*pro2.final.se.q1.1[5]
c(var.est.q1.1,var.LowerCI.q1.1,var.UpperCI.q1.1,var.UpperCI.q1.1-var.LowerCI.q1.1)




CI.Int.q1.1 <- c()
CI.Int.q1.1 <- cbind(Alpha=as.numeric(c(alpha.est.q1.1,alpha.LowerCI.q1.1,alpha.UpperCI.q1.1,alpha.UpperCI.q1.1-alpha.LowerCI.q1.1)),Beta=as.numeric(c(beta.est.q1.1,beta.LowerCI.q1.1,beta.UpperCI.q1.1,beta.UpperCI.q1.1-beta.LowerCI.q1.1)),Gamma=as.numeric(c(gamma.est.q1.1,gamma.LowerCI.q1.1,gamma.UpperCI.q1.1,gamma.UpperCI.q1.1-gamma.LowerCI.q1.1)),Theta=as.numeric(c(theta.est.q1.1,theta.LowerCI.q1.1,theta.UpperCI.q1.1,theta.UpperCI.q1.1-theta.LowerCI.q1.1)),Var=as.numeric(c(var.est.q1.1,var.LowerCI.q1.1,var.UpperCI.q1.1,var.UpperCI.q1.1-var.LowerCI.q1.1)))
row.names(CI.Int.q1.1) <- c("Point Estimate","Lower_95%_CI","Upper_95% CI","Width of CI")
CI.Int.q1.1
round(CI.Int.q1.1,4)
#t(round(CI.Int.q1.1,4))


### Case-2 : Question-1


pro2.final.init.q1.2 <- c(0.2222222,	2.0000000,	0.2222222,	0.1,	0.5) # Initial values
pro2.final.init.q1.2
round(pro2.final.init.q1.2,4)

pro2.final.mod.q1.2 <- nlm(f = Pro2.negLhood.MLE,p = pro2.final.init.q1.2,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
pro2.final.mod.q1.2
pro2.final.est.q1.2 <- pro2.final.mod.q1.2$estimate
pro2.final.est.q1.2
pro2.final.hessian.q1.2 <- pro2.final.mod.q1.2$hessian
pro2.final.vcov.q1.2 <- solve(pro2.final.hessian.q1.2)
pro2.final.vcov.q1.2
diag(pro2.final.vcov.q1.2)
pro2.final.se.q1.2 <- sqrt(diag(pro2.final.vcov.q1.2))
pro2.final.se.q1.2
pro2.final.maxLL.q1.2 <- pro2.final.mod.q1.2$minimum
pro2.final.maxLL.q1.2 # MAximum Log-likelihood value
round(pro2.final.maxLL.q1.2,4) # MAx LL
pro2.final.gradient.q1.2 <- pro2.final.mod.q1.2$gradient
pro2.final.gradient.q1.2 # Gradients

## For alpha # Q1- case-2

alpha.est.q1.2 <- pro2.final.est.q1.2[1]
alpha.LowerCI.q1.2 <- alpha.est.q1.2 - 1.96*pro2.final.se.q1.2[1]
alpha.UpperCI.q1.2 <- alpha.est.q1.2 + 1.96*pro2.final.se.q1.2[1]
c(alpha.est.q1.2,alpha.LowerCI.q1.2,alpha.UpperCI.q1.2,alpha.UpperCI.q1.2-alpha.LowerCI.q1.2)
#round(alpha.LowerCI.q1.2,4)


## For beta # Q1- case-2

beta.est.q1.2 <- pro2.final.est.q1.2[2]
beta.LowerCI.q1.2 <- beta.est.q1.2 - 1.96*pro2.final.se.q1.2[2]
beta.UpperCI.q1.2 <- beta.est.q1.2 + 1.96*pro2.final.se.q1.2[2]
c(beta.est.q1.2,beta.LowerCI.q1.2,beta.UpperCI.q1.2,beta.UpperCI.q1.2-beta.LowerCI.q1.2)


## For gamma # Q1- case-2

gamma.est.q1.2 <- pro2.final.est.q1.2[3]
gamma.LowerCI.q1.2 <- gamma.est.q1.2 - 1.96*pro2.final.se.q1.2[3]
gamma.UpperCI.q1.2 <- gamma.est.q1.2 + 1.96*pro2.final.se.q1.2[3]
c(gamma.est.q1.2,gamma.LowerCI.q1.2,gamma.UpperCI.q1.2,gamma.UpperCI.q1.2-gamma.LowerCI.q1.2)


## For theta # # Q1- case-2

theta.est.q1.2 <- pro2.final.est.q1.2[4]
theta.LowerCI.q1.2 <- theta.est.q1.2 - 1.96*pro2.final.se.q1.2[4]
theta.UpperCI.q1.2 <- theta.est.q1.2 + 1.96*pro2.final.se.q1.2[4]
c(theta.est.q1.2,theta.LowerCI.q1.2,theta.UpperCI.q1.2,theta.UpperCI.q1.2-theta.LowerCI.q1.2)


## For variance ## # Q1- case-2

var.est.q1.2 <- pro2.final.est.q1.2[5]
var.LowerCI.q1.2 <- var.est.q1.2 - 1.96*pro2.final.se.q1.2[5]
var.UpperCI.q1.2 <- var.est.q1.2 + 1.96*pro2.final.se.q1.2[5]
c(var.est.q1.2,var.LowerCI.q1.2,var.UpperCI.q1.2,var.UpperCI.q1.2-var.LowerCI.q1.2)




CI.Int.q1.2 <- c()
CI.Int.q1.2 <- cbind(Alpha=as.numeric(c(alpha.est.q1.2,alpha.LowerCI.q1.2,alpha.UpperCI.q1.2,alpha.UpperCI.q1.2-alpha.LowerCI.q1.2)),Beta=as.numeric(c(beta.est.q1.2,beta.LowerCI.q1.2,beta.UpperCI.q1.2,beta.UpperCI.q1.2-beta.LowerCI.q1.2)),Gamma=as.numeric(c(gamma.est.q1.2,gamma.LowerCI.q1.2,gamma.UpperCI.q1.2,gamma.UpperCI.q1.2-gamma.LowerCI.q1.2)),Theta=as.numeric(c(theta.est.q1.2,theta.LowerCI.q1.2,theta.UpperCI.q1.2,theta.UpperCI.q1.2-theta.LowerCI.q1.2)),Var=as.numeric(c(var.est.q1.2,var.LowerCI.q1.2,var.UpperCI.q1.2,var.UpperCI.q1.2-var.LowerCI.q1.2)))
row.names(CI.Int.q1.2) <- c("Point Estimate","Lower_95%_CI","Upper_95% CI","Width of CI")
CI.Int.q1.2
round(CI.Int.q1.2,4)


### Choice of Final Estimates

pro2.final.gradient.q1.2<=pro2.final.gradient.q1.1

# Case-2 is the Final MLE Estimates


# 1.d) Plot of data overlaid with Expectation function, and 95% confidence band.

# Case-1 : Question-1

# Fitted Expectation Function

alpha.est.q1.1*XYdata.2[,1]^(beta.est.q1.1*XYdata.2[,1]^(-gamma.est.q1.1))
par.q1.1 <- c(alpha.est.q1.1,beta.est.q1.1,gamma.est.q1.1,theta.est.q1.1,var.est.q1.1)
par.q1.1
mu.i.est.q1.1 <- mu.i.func(X = XYdata.2[,1],par = par.q1.1)
mu.i.est.q1.1

# Standard error of Expectation function

der.mat.q1.1 <- der.func(X = XYdata.2[,1],par = par.q1.1[1:3])
der.mat.q1.1
vcov.mat.q1.1 <- pro2.final.vcov.q1.1[1:3,1:3]
vcov.mat.q1.1
der.mat.q1.1%*%vcov.mat.q1.1%*%t(der.mat.q1.1)

diag(der.mat.q1.1%*%vcov.mat.q1.1%*%t(der.mat.q1.1))
mu.i.var.q1.1 <- diag(der.mat.q1.1%*%vcov.mat.q1.1%*%t(der.mat.q1.1))
mu.i.var.q1.1

mu.i.se.q1.1 <- sqrt(mu.i.var.q1.1)
mu.i.se.q1.1

mu.i.LowerCI.q1.1 <- mu.i.est.q1.1 - 1.96*mu.i.se.q1.1
mu.i.UpperCI.q1.1 <- mu.i.est.q1.1 + 1.96*mu.i.se.q1.1

cbind(mu.i.LowerCI.q1.1,mu.i.est.q1.1,mu.i.UpperCI.q1.1)

windows(width = 10,height = 8)
plot(x = s601.proj2.data$x,y = s601.proj2.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.est.q1.1[order(s601.proj2.data$x)],col="red",lwd=2)
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y =  mu.i.LowerCI.q1.1[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.UpperCI.q1.1[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
dev.off()


# Case-2 : Question-1

# Fitted Expectation Function

alpha.est.q1.2*XYdata.2[,1]^(beta.est.q1.2*XYdata.2[,1]^(-gamma.est.q1.2))
par.q1.2 <- c(alpha.est.q1.2,beta.est.q1.2,gamma.est.q1.2,theta.est.q1.2,var.est.q1.2)
par.q1.2
mu.i.est.q1.2 <- mu.i.func(X = XYdata.2[,1],par = par.q1.2)
mu.i.est.q1.2

# Standard error of Expectation function

der.mat.q1.2 <- der.func(X = XYdata.2[,1],par = par.q1.2[1:3])
der.mat.q1.2
vcov.mat.q1.2 <- pro2.final.vcov.q1.2[1:3,1:3]
vcov.mat.q1.2
der.mat.q1.2%*%vcov.mat.q1.2%*%t(der.mat.q1.2)

diag(der.mat.q1.2%*%vcov.mat.q1.2%*%t(der.mat.q1.2))
mu.i.var.q1.2 <- diag(der.mat.q1.2%*%vcov.mat.q1.2%*%t(der.mat.q1.2))
mu.i.var.q1.2

mu.i.se.q1.2 <- sqrt(mu.i.var.q1.2)
mu.i.se.q1.2

mu.i.LowerCI.q1.2 <- mu.i.est.q1.2 - 1.96*mu.i.se.q1.2
mu.i.UpperCI.q1.2 <- mu.i.est.q1.2 + 1.96*mu.i.se.q1.2

cbind(mu.i.LowerCI.q1.2,mu.i.est.q1.2,mu.i.UpperCI.q1.2)

windows(width = 10,height = 8)
plot(x = s601.proj2.data$x,y = s601.proj2.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.est.q1.2[order(s601.proj2.data$x)],col="red",lwd=2)
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y =  mu.i.LowerCI.q1.2[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.UpperCI.q1.2[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
dev.off()



### 1.e) Large-scale model component: Summary of Parameter Results, Fitted Curve, 95% Confidence Band

### Results from Part-1 Q.1. (Rounded off to 4 places of decimal)

### GLS Estimates and CI

#Parameters      alpha  beta   gamma
#Point Estimate 0.5063 1.2220 0.1614
#Lower 95% CI   0.3389 1.0198 0.1519
#Upper 95% CI   0.6736 1.4241 0.1708
# Width of CI   0.3347 0.4043 0.0189  

# plot(x = s601.proj2.data$x,y = s601.proj2.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")
#    lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.est[order(s601.proj2.data$x)],col="red",lwd=2)
#    lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y =  mu.i.LowerCI[order(s601.proj2.data$x)],col="green",lty=2,lwd=1.5)
#    lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.UpperCI[order(s601.proj2.data$x)],col="green",lty=2,lwd=1.5)
 
  


### Results from Part-2 Q.1. (Rounded off to 4 places of decimal)

# Case-1

round(CI.Int.q1.1,4)[,1:3]

head(cbind(mu.i.LowerCI.q1.1,mu.i.est.q1.1,mu.i.UpperCI.q1.1))

# Case-2

round(CI.Int.q1.2,4)[,1:3]

head(cbind(mu.i.LowerCI.q1.2,mu.i.est.q1.2,mu.i.UpperCI.q1.2))




### 1.f) Small-scale model component: Summary of Parameter Results

### Results from Part-1 Q.1. (Rounded off to 4 places of decimal)

### theta, variance - Optimal estimating equations and MOM

#var-MOM = 1.0149
#var-optimal=0.6295

#theta-optimal=0.5422



### Results from Part-2 Q.1. (Rounded off to 4 places of decimal)

# Case-1

round(CI.Int.q1.1,4)[,4:5]


# Case-2

round(CI.Int.q1.2,4)[,4:5]



### Q.2.


## Log-composite Likelihood

"Pro2.Lhood.CLE"<- function(X,Y,pars){
  alpha <- pars[1]
  beta <- pars[2]
  gamma <- pars[3]
  theta <- pars[4]
  sigma.sq <- pars[5]
  n <- length(Y)
  mu.i <- alpha*X^(beta*(X^((-1)*gamma)))
  
  # cat("Mu.i",mu.i,fill=T)
  # cat("Y[-1]",Y[-1],fill=T)
  # cat("Mean of Y_i 2 to n",mu.i[-1]+theta*(Y[-n]-mu.i[-n]),fill=T)
  # cat("Y[1]",Y[1],fill=T)
  # cat("mu.i.1",mu.i[1],fill=T)
  # cat("sigma.sq",sigma.sq,"sd",sqrt(sigma.sq),fill=T)
  # cat("dnorm.1",dnorm(x = Y[1],mean = mu.i[1],sd = sqrt(sigma.sq)),fill=T)
  # cat("dnorm.2",dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq)),fill=T)
  # cat("prod.dnorm.2",prod(dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq))),fill=T)
  
  Lhood <- prod(dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq)))
  logLhood <- log(x = Lhood,base = exp(1))
  
  logLhood.1 <- ((-1)*((n-1)/2))*log(2*pi*sigma.sq,base = exp(1))-sum((Y[-1]-(mu.i[-1]+theta*(Y[-n]-mu.i[-n])))^2/(2*sigma.sq))
  
  logLhood.2 <- sum(log(dnorm(x = Y[-1],mean = mu.i[-1]+theta*(Y[-n]-mu.i[-n]),sd = sqrt(sigma.sq))))
  
  sum.2.n <- 0
  for (i in 2:n){
    sum.2.n <- sum.2.n + (Y[i]-(mu.i[i]+theta*(Y[i-1]-mu.i[i-1])))^2
  }
  
  logLhood.3 <- (-(n-1)/2)*log(2*pi*sigma.sq)- sum.2.n/(2*sigma.sq)
  
  return (list(Lhood=Lhood,logLhood=logLhood,logLhood.1=logLhood.1,logLhood.2=logLhood.2,logLhood.3=logLhood.3))
}

"Pro2.negLhood.CLE"<-function(X,Y,pars){
  
  #neg.Lhood <- (-1)*Pro2.Lhood.MLE(X,Y,pars)$logLhood
  neg.Lhood <- (-1)*Pro2.Lhood.CLE(X,Y,pars)$logLhood.1
  #neg.Lhood <- (-1)*Pro2.Lhood.MLE(X,Y,pars)$logLhood.2
  #neg.Lhood <- (-1)*Pro2.Lhood.MLE(X,Y,pars)$logLhood.3
  
  return (neg.Lhood)
}





## Fit Model using log-Composite Likelihood

# Initial Values - Q.2. Case-1

pro2.final.init.q2.1 <- c(0.2222222,	2,	0.2222222,	0.99,	0.61) # Initial values
pro2.final.init.q2.1
round(pro2.final.init.q2.1,4)

# Initial Values - Q.2. Case-2

pro2.final.init.q2.2 <- c(0.2222222,	2.0000000,	0.2222222,	0.1,	0.5) # Initial values
pro2.final.init.q2.2

## Q.2. Case-1: Composite Likelihood Estimation

Pro2.Lhood.CLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.final.init.q2.1)

optim(par = pro2.final.init.q2.1,fn = Pro2.negLhood.CLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")

optim(par = pro2.final.init.q2.1,fn = Pro2.negLhood.CLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

nlm(f = Pro2.negLhood.CLE,p = pro2.final.init.q2.1,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])



## Q.2. Case-2: Composite Likelihood Estimation

Pro2.Lhood.CLE(X = XYdata.2[,1],Y = XYdata.2[,2],pars = pro2.final.init.q2.2)

optim(par = pro2.final.init.q2.2,fn = Pro2.negLhood.CLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")

optim(par = pro2.final.init.q2.2,fn = Pro2.negLhood.CLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

nlm(f = Pro2.negLhood.CLE,p = pro2.final.init.q2.2,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])


# Final Choice of Maximum Composite Likelihood Estimates.

# We will go by Case-1 estimates, since both optim() and nlm() converged correctly with that initial value combination.

# But, among optim(Nelder-mead), optim(BFGS) and nlm(), we may have to see which one gives narrower CI for the parameters, and also the plot of expectation function.

# For now, we will start with nlm(), since we chose this method for MLE estimates, as well.

pro2.final.mod.q2.1 <- nlm(f = Pro2.negLhood.CLE,p = pro2.final.init.q2.1,hessian = T,X=XYdata.2[,1],Y=XYdata.2[,2])
pro2.final.mod.q2.1
pro2.final.est.q2.1 <- pro2.final.mod.q2.1$estimate
pro2.final.est.q2.1
#pro2.final.hessian.q2.1 <- pro2.final.mod.q2.1$hessian
#pro2.final.vcov.q2.1 <- solve(pro2.final.hessian.q2.1)
#pro2.final.vcov.q2.1
#diag(pro2.final.vcov.q2.1)
#pro2.final.se.q2.1 <- sqrt(diag(pro2.final.vcov.q2.1))
#pro2.final.se.q2.1
pro2.final.maxLL.q2.1 <- pro2.final.mod.q2.1$minimum
pro2.final.maxLL.q2.1 # MAximum Log-likelihood value
round(pro2.final.maxLL.q2.1,4)
pro2.final.gradient.q2.1 <- pro2.final.mod.q2.1$gradient
pro2.final.gradient.q2.1 # Gradients

## Inference of Parameters using Composite Likelihood

## Following Inference with respect to Final Model Selected (pro2.final.mod.q2.1) for composite Likelihood

"CLE.theta._i" <- function(X,Y,par,maxf){
  est <- nlm(f = maxf,p = par,hessian = T,X=X,Y=Y)$estimate
  return (est)
}

CLE.cov <- matrix(data = 0,nrow = length(pro2.final.est.q2.1),ncol = length(pro2.final.est.q2.1))
CLE.cov
theta.full.CLE <- pro2.final.est.q2.1
theta.full.CLE
CLE.cov.11 <- c()
for (i in 1:nrow(XYdata.2)){
  theta._i.CLE <- CLE.theta._i(X = XYdata.2[-i,1],Y = XYdata.2[-i,2],par = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE)
  cat("Estimates_i",theta._i.CLE,fill=T,"\n")
  diff.est.1 <-  theta._i.CLE - theta.full.CLE
  diff.est.1.trans <- t(diff.est.1)
  cat("Difference",diff.est.1,fill=T,"\n")
  cat(diff.est.1.trans,"\n")
  matrix.cov.i <- as.matrix(diff.est.1%*%diff.est.1.trans)
  cat("Matrix.prod",fill=T,"\n")
  print(matrix.cov.i)
  CLE.cov <- CLE.cov + matrix.cov.i
  cat("CLE.cov.sum","\n")
  print(CLE.cov)
  CLE.cov.11[i] <- CLE.cov[1,1]
}

CLE.cov.11
CLE.theta._i(X = XYdata.2[-3,1],Y = XYdata.2[-3,2],par = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE)
#CLE.theta._i(X = XYdata.2[-4,1],Y = XYdata.2[-4,2],par = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE)
CLE.theta._i(X = XYdata.2[-6,1],Y = XYdata.2[-6,2],par = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE)
# optim(X = XYdata.2[-6,1],Y = XYdata.2[-6,2],par = pro2.final.init.q2.1,fn =  Pro2.negLhood.CLE,method = "BFGS",hessian = T)
# optim(X = XYdata.2[-3,1],Y = XYdata.2[-3,2],par = pro2.final.init.q2.1,fn =  Pro2.negLhood.CLE,method = "BFGS",hessian = T)
# optim(X = XYdata.2[-4,1],Y = XYdata.2[-4,2],par = pro2.final.init.q2.1,fn =  Pro2.negLhood.CLE,method = "BFGS",hessian = T)

CLE.cov
((nrow(XYdata.2)-1)/nrow(XYdata.2))*CLE.cov
diag(((nrow(XYdata.2)-1)/nrow(XYdata.2))*CLE.cov)
sqrt(diag(((nrow(XYdata.2)-1)/nrow(XYdata.2))*CLE.cov))

CLE.cov.final.q2.1 <- ((nrow(XYdata.2)-1)/nrow(XYdata.2))*CLE.cov
CLE.cov.final.q2.1

#solve(CLE.cov.final.q2.1)

pro2.final.se.q2.1 <- sqrt(diag(CLE.cov.final.q2.1))
pro2.final.se.q2.1




## For alpha # Q2- case-1

alpha.est.q2.1 <- pro2.final.est.q2.1[1]
alpha.LowerCI.q2.1 <- alpha.est.q2.1 - 1.96*pro2.final.se.q2.1[1]
alpha.UpperCI.q2.1 <- alpha.est.q2.1 + 1.96*pro2.final.se.q2.1[1]
c(alpha.est.q2.1,alpha.LowerCI.q2.1,alpha.UpperCI.q2.1,alpha.UpperCI.q2.1-alpha.LowerCI.q2.1)

## For beta # Q2- case-1

beta.est.q2.1 <- pro2.final.est.q2.1[2]
beta.LowerCI.q2.1 <- beta.est.q2.1 - 1.96*pro2.final.se.q2.1[2]
beta.UpperCI.q2.1 <- beta.est.q2.1 + 1.96*pro2.final.se.q2.1[2]
c(beta.est.q2.1,beta.LowerCI.q2.1,beta.UpperCI.q2.1,beta.UpperCI.q2.1-beta.LowerCI.q2.1)


## For gamma # Q2- case-1

gamma.est.q2.1 <- pro2.final.est.q2.1[3]
gamma.LowerCI.q2.1 <- gamma.est.q2.1 - 1.96*pro2.final.se.q2.1[3]
gamma.UpperCI.q2.1 <- gamma.est.q2.1 + 1.96*pro2.final.se.q2.1[3]
c(gamma.est.q2.1,gamma.LowerCI.q2.1,gamma.UpperCI.q2.1,gamma.UpperCI.q2.1-gamma.LowerCI.q2.1)


## For theta # # Q2- case-1

theta.est.q2.1 <- pro2.final.est.q2.1[4]
theta.LowerCI.q2.1 <- theta.est.q2.1 - 1.96*pro2.final.se.q2.1[4]
theta.UpperCI.q2.1 <- theta.est.q2.1 + 1.96*pro2.final.se.q2.1[4]
c(theta.est.q2.1,theta.LowerCI.q2.1,theta.UpperCI.q2.1,theta.UpperCI.q2.1-theta.LowerCI.q2.1)


## For variance ## # Q2- case-1

var.est.q2.1 <- pro2.final.est.q2.1[5]
var.LowerCI.q2.1 <- var.est.q2.1 - 1.96*pro2.final.se.q2.1[5]
var.UpperCI.q2.1 <- var.est.q2.1 + 1.96*pro2.final.se.q2.1[5]
c(var.est.q2.1,var.LowerCI.q2.1,var.UpperCI.q2.1,var.UpperCI.q2.1-var.LowerCI.q2.1)




CI.Int.q2.1 <- c()
CI.Int.q2.1 <- cbind(Alpha=as.numeric(c(alpha.est.q2.1,alpha.LowerCI.q2.1,alpha.UpperCI.q2.1,alpha.UpperCI.q2.1-alpha.LowerCI.q2.1)),Beta=as.numeric(c(beta.est.q2.1,beta.LowerCI.q2.1,beta.UpperCI.q2.1,beta.UpperCI.q2.1-beta.LowerCI.q2.1)),Gamma=as.numeric(c(gamma.est.q2.1,gamma.LowerCI.q2.1,gamma.UpperCI.q2.1,gamma.UpperCI.q2.1-gamma.LowerCI.q2.1)),Theta=as.numeric(c(theta.est.q2.1,theta.LowerCI.q2.1,theta.UpperCI.q2.1,theta.UpperCI.q2.1-theta.LowerCI.q2.1)),Var=as.numeric(c(var.est.q2.1,var.LowerCI.q2.1,var.UpperCI.q2.1,var.UpperCI.q2.1-var.LowerCI.q2.1)))
row.names(CI.Int.q2.1) <- c("Point Estimate","Lower_95%_CI","Upper_95% CI","Width of CI")
CI.Int.q2.1
round(CI.Int.q2.1,4)


### Q2: Case-1: Inference of Composite Likelihood: Expectation Function with respect to the 

# Fitted Expectation Function

par.q2.1 <- pro2.final.est.q2.1
par.q2.1
mu.i.est.q2.1 <- mu.i.func(X = XYdata.2[,1],par = par.q2.1)
mu.i.est.q2.1

# Standard error of Expectation function

der.mat.q2.1 <- der.func(X = XYdata.2[,1],par = par.q2.1[1:3])
der.mat.q2.1
vcov.mat.q2.1 <- CLE.cov.final.q2.1[1:3,1:3]
vcov.mat.q2.1
der.mat.q2.1%*%vcov.mat.q2.1%*%t(der.mat.q2.1)

diag(der.mat.q2.1%*%vcov.mat.q2.1%*%t(der.mat.q2.1))
mu.i.var.q2.1 <- diag(der.mat.q2.1%*%vcov.mat.q2.1%*%t(der.mat.q2.1))
mu.i.var.q2.1

mu.i.se.q2.1 <- sqrt(mu.i.var.q2.1)
mu.i.se.q2.1

mu.i.LowerCI.q2.1 <- mu.i.est.q2.1 - 1.96*mu.i.se.q2.1
mu.i.UpperCI.q2.1 <- mu.i.est.q2.1 + 1.96*mu.i.se.q2.1

cbind(mu.i.LowerCI.q2.1,mu.i.est.q2.1,mu.i.UpperCI.q2.1)

windows(width = 10,height = 8)
plot(x = s601.proj2.data$x,y = s601.proj2.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(mu.i.LowerCI.q2.1),max(mu.i.UpperCI.q2.1)))
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.est.q2.1[order(s601.proj2.data$x)],col="red",lwd=2)
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y =  mu.i.LowerCI.q2.1[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.i.UpperCI.q2.1[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
dev.off()



###########################

# Checking Jackknife Matrix using Estimates from optim BFGS, Nelder-Mead, based on Initial Values of Q.2 Case-1


pro2.final.mod.q2.1.1 <- optim(par = pro2.final.init.q2.1,fn = Pro2.negLhood.CLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "Nelder-Mead")

pro2.final.mod.q2.1.2 <- optim(par = pro2.final.init.q2.1,fn = Pro2.negLhood.CLE,X=XYdata.2[,1],Y=XYdata.2[,2],hessian = T,method = "BFGS")

"CLE.theta._i.1" <- function(X,Y,par,maxf,method){
  est1 <- NULL
  est2 <- NULL
  if (method=="BFGS"){
    est2 <- optim(fn = maxf,par = par,hessian = T,X=X,Y=Y,method="BFGS")$par
  } else {
  est1 <- optim(fn = maxf,par = par,hessian = T,X=X,Y=Y,method="Nelder-Mead")$par}
  return (list(est.NM=est1,est.BFGS=est2))
}

"Jackknife.cov" <- function(X,Y,init,maxf,opt.func,opt.method=NULL){
  
CLE.cov.mat <- matrix(data = 0,nrow = length(init),ncol = length(init))
full.CLE <- c()
chk <- 0
if (opt.func=="optim" && is.null(opt.method)==F && (opt.method=="Nelder-Mead" || opt.method=="BFGS")){
  if (opt.method=="Nelder-Mead"){
    full.CLE <- CLE.theta._i.1(X = X,Y = Y,par = init,maxf = maxf,method = "Nelder-Mead")$est.NM
    chk <- 1
  } else {
    full.CLE <- CLE.theta._i.1(X = X,Y = Y,par = init,maxf = maxf,method = "BFGS")$est.BFGS
    chk <- 2
  }
} else if (opt.func=="nlm"){
  full.CLE <- CLE.theta._i(X = XYdata.2[,1],Y = XYdata.2[,2],par = init,maxf = maxf)
  chk <- 3
} else {
  stop()
}
n <- length(X)
theta._i.CLE.est <- c()

for (i in 1:n){
  if (chk==1){
    theta._i.CLE.est <- CLE.theta._i.1(X = X[-i],Y = Y[-i],par = init,maxf = maxf,method = "Nelder-Mead")$est.NM
  } else if (chk==2){
    theta._i.CLE.est <- CLE.theta._i.1(X = X[-i],Y = Y[-i],par = init,maxf = maxf,method = "BFGS")$est.BFGS
  } else if (chk==3){
    theta._i.CLE.est <- CLE.theta._i(X = X[-i],Y = Y[-i],par = init,maxf = maxf)
  } else {
    break
  }
  #cat("Estimates_i",theta._i.CLE,fill=T,"\n")
  diff.1 <-  theta._i.CLE.est - full.CLE
  diff.1.trans <- t(diff.1)
  #cat("Difference",diff.est.1,fill=T,"\n")
  #cat(diff.est.1.trans,"\n")
  matrix.i <- as.matrix(diff.1%*%diff.1.trans)
  #cat("Matrix.prod",fill=T,"\n")
  #print(matrix.cov.i)
  CLE.cov.mat <- CLE.cov.mat + matrix.i
  #cat("CLE.cov.sum","\n")
  #print(CLE.cov)
  #CLE.cov.11[i] <- CLE.cov[1,1]
}

CLE.cov.final<- ((n-1)/n)*CLE.cov.mat
CLE.cov.final.q2.1
return (list(covMat=CLE.cov.final,cov.diag=diag(CLE.cov.final),se.diag=sqrt(diag(CLE.cov.final))))
}

Jackknife.cov(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "nlm")

Jackknife.cov(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "BFGS")

Jackknife.cov(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "Nelder-Mead")

"CLE.Infer" <- function(X,Y,init,maxf,opt.func,opt.method=NULL){
  par.est <- c()
  se.est <- c()
  if (opt.func=="optim" && is.null(opt.method)==F && (opt.method=="Nelder-Mead" || opt.method=="BFGS")){
    if (opt.method=="Nelder-Mead"){
      par.est <- CLE.theta._i.1(X = X,Y = Y,par = init,maxf = maxf,method = "Nelder-Mead")$est.NM
      se.est <- Jackknife.cov(X = X,Y = Y,init = init,maxf = maxf,opt.func = "optim",opt.method = "Nelder-Mead")$se.diag
      #chk <- 1
    } else {
      par.est <- CLE.theta._i.1(X = X,Y = Y,par = init,maxf = maxf,method = "BFGS")$est.BFGS
      se.est <- Jackknife.cov(X = X,Y = Y,init = init,maxf = maxf,opt.func = "optim",opt.method = "BFGS")$se.diag
      #chk <- 2
    }
  } else if (opt.func=="nlm"){
    par.est <- CLE.theta._i(X = XYdata.2[,1],Y = XYdata.2[,2],par = init,maxf = maxf)
    se.est <- Jackknife.cov(X = X,Y = Y,init = init,maxf = maxf,opt.func = "nlm")$se.diag
    #chk <- 3
  } else {
    stop()
  }
  LowerCI <- c()
  UpperCI <- c()
  WidthCI <- c()
  for (i in 1:length(par.est)){
    LowerCI[i] <- par.est[i] - 1.96*se.est[i]
    UpperCI[i] <- par.est[i] + 1.96*se.est[i]
    WidthCI[i] <- UpperCI[i] - LowerCI[i]
  }
  
  CI.Int.q2 <- c()
  CI.Int.q2 <- cbind(alpha=c(par.est[1],LowerCI[1],UpperCI[1],WidthCI[1]),beta=c(par.est[2],LowerCI[2],UpperCI[2],WidthCI[2]),gamma=c(par.est[3],LowerCI[3],UpperCI[3],WidthCI[3]),theta=c(par.est[4],LowerCI[4],UpperCI[4],WidthCI[4]),var=c(par.est[5],LowerCI[5],UpperCI[5],WidthCI[5]))
  
  return (list(CI.Int=CI.Int.q2,CI.Int.round=round(CI.Int.q2,4)))
} 

CLE.Infer(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "nlm")

CLE.Infer(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "BFGS")

CLE.Infer(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "Nelder-Mead") # Variance Lower CI is negative


"Mu.Exp.func" <- function(X,Y,init,maxf,opt.func,opt.method=NULL){
  
  par.est <- c()
  vcov.mat <- c()
  der.mat <- c()
  mu.est <- c()
  
  if (opt.func=="optim" && is.null(opt.method)==F && (opt.method=="Nelder-Mead" || opt.method=="BFGS")){
    if (opt.method=="Nelder-Mead"){
      par.est <- CLE.theta._i.1(X = X,Y = Y,par = init,maxf = maxf,method = "Nelder-Mead")$est.NM
      vcov.mat <- Jackknife.cov(X = X,Y = Y,init = init,maxf = maxf,opt.func = "optim",opt.method = "Nelder-Mead")$covMat
      
      #chk <- 1
    } else {
      par.est <- CLE.theta._i.1(X = X,Y = Y,par = init,maxf = maxf,method = "BFGS")$est.BFGS
      vcov.mat <- Jackknife.cov(X = X,Y = Y,init = init,maxf = maxf,opt.func = "optim",opt.method = "BFGS")$covMat
      #chk <- 2
    }
  } else if (opt.func=="nlm"){
    par.est <- CLE.theta._i(X = XYdata.2[,1],Y = XYdata.2[,2],par = init,maxf = maxf)
    vcov.mat <- Jackknife.cov(X = X,Y = Y,init = init,maxf = maxf,opt.func = "nlm")$covMat
    #chk <- 3
  } else {
    stop()
  }
  der.mat <- der.func(X = X,par = par.est[1:3])
  mu.est <-  mu.i.func(X = X,par = par.est)
  mu.var <-  diag(der.mat%*%vcov.mat[1:3,1:3]%*%t(der.mat))
  mu.se <- sqrt(mu.var)
  
  mu.LowerCI <- mu.est - 1.96*mu.se
  mu.UpperCI <- mu.est + 1.96*mu.se
  
  mu.res <- cbind(LowerCI.mu=mu.LowerCI,MU.est=mu.est,UpperCI.mu=mu.UpperCI)
  
  plot(x = s601.proj2.data$x,y = s601.proj2.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(mu.LowerCI),max(mu.UpperCI)))
  lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.est[order(s601.proj2.data$x)],col="red",lwd=2)
  lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y =  mu.LowerCI[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
  lines(x = s601.proj2.data$x[order(s601.proj2.data$x)],y = mu.UpperCI[order(s601.proj2.data$x)],col="green",lty=1,lwd=1.5)
  
  return (list(muResult=mu.res))
}

Mu.Exp.func(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "nlm")

Mu.Exp.func(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "BFGS")

Mu.Exp.func(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.1,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "Nelder-Mead")



#############################


## Testing Jackknife Matrix using Q2: Case-2 Initial values

# for (a in 1:25){print(CLE.theta._i(X = XYdata.2[-a,1],Y = XYdata.2[-a,2],par = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE))}

Mu.Exp.func(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "nlm")

Mu.Exp.func(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "BFGS")

Mu.Exp.func(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "Nelder-Mead")

CLE.Infer(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "nlm")

CLE.Infer(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "BFGS") # Variance Negative LowerCI

CLE.Infer(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "Nelder-Mead")

Jackknife.cov(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "nlm")

Jackknife.cov(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "BFGS") # Absurd

Jackknife.cov(X = XYdata.2[,1],Y = XYdata.2[,2],init = pro2.final.init.q2.2,maxf = Pro2.negLhood.CLE,opt.func = "optim",opt.method = "Nelder-Mead")



##########

## To Do: Alternative Method for Computation of Confidence Intervals of Parameters in the case of Composite Likelihood, using Parametric Bootstrap.


###########




### Q.3.


# Non-Linear Regression Model from Part-1

pars.Part_1.q3 <- initvals[[1]]
pars.Part_1.q3
mod.GLS.Part_1.q3 <- Q4.GLS.1(par = pars.Part_1.q3,XYdata = XYdata.2) # GLS Model
mod.GLS.Part_1.q3
mod.GLS.Part_1.q3$pars # estimates

mu.i.GLS.q3 <- mu.i.func(X = XYdata.2[,1],par =mod.GLS.Part_1.q3$pars) # Fitted Y.hat or mu
mu.i.GLS.q3

mod.GLS.Part_1.resid <- XYdata.2[,2] - mu.i.GLS.q3
mod.GLS.Part_1.resid # Raw Residuals or W_i

## Loading the given simulated Dataset of K-S Statistic and Partial Lag-1 Autocorrelation coefficient.

proj2.data.q3 <- read.table(file = "reference_distributions.txt",header = T)
proj2.data.q3
nrow(proj2.data.q3)
str(proj2.data.q3)



# 3.a) diagnostic plot using generalized residuals from your fitted nonlinear regression model from Part-1

# Calculation of Generalized Residuals

X.q3 <- XYdata.2[,2]  # Original Response Data (say, X)

# Generalized Residuals - Version-1
gen.res <- pnorm(q = XYdata.2[,2],mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)) # CDF of X (Fx(x))
gen.res  #CDF of Fitted Y.hat
plot(gen.res)
abline(h=mean(gen.res))
plot(ecdf(gen.res),verticals = T)
Y.q3 <- gen.res # Y = Fx(x)
F.y.q3 <- cumsum(table(Y.q3))/sum(table(Y.q3)) # CDF of Fy(y)
F.y.q3
Y.q3.ecdf <- ecdf(Y.q3)
Y.q3.ecdf(Y.q3)
hist(X.q3)
hist(sort(X.q3))
plot(F.y.q3)
#plot(y = F.y.q3,x = Y.q3)
plot(x = sort(Y.q3),y = F.y.q3)
plot(x = Y.q3,y = Y.q3.ecdf(Y.q3))

plot(x = sort(Y.q3),y = F.y.q3)
lines(x = sort(Y.q3),y = Y.q3.ecdf(Y.q3))

par(mfrow=c(1,3))
hist(sort(X.q3))
plot(x = sort(Y.q3),y = F.y.q3)
plot(x = Y.q3,y = Y.q3.ecdf(Y.q3))
par(mfrow=c(1,1))


### Generalized Residuals - 2nd version

gen.res.1 <- pnorm(q = mod.GLS.Part_1.resid,mean = 0,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)) # CDF of RAW Residuals
gen.res.1
identical(x = gen.res,y = gen.res.1)

plot(pnorm(q = XYdata.2[,2],mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)))
mean(pnorm(q = XYdata.2[,2],mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)))
#curve(expr = pnorm(q = XYdata.2[,2],mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)),from = -Inf,to = Inf,n = nrow(XYdata.2))
curve(expr = pnorm,from = min(XYdata.2[,2]),to = max(XYdata.2[,2]),n = nrow(XYdata.2))
curve(expr = pnorm,from = -max(XYdata.2[,2]),to = max(XYdata.2[,2]),n = nrow(XYdata.2))

q3.unif <- runif(n = 25,min = 0.0001,max = 0.9999)
q3.unif.cdf <- punif(q = q3.unif,min = 0,max = 1)
plot(q3.unif.cdf)
#pnorm(q = runif(n = 25,min = 0.0001,max = 0.9999),mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2))
plot.ecdf(q3.unif.cdf)




# 3.b) Identify the model I produced simulated data sets from.

### Joint Distribution of Fitted Model F(y|theta.hat)






# 3.c.) how the Kolmogorov test statistic values were computed


# Procedure for K-S Test Statistic # Comparing Empirical CDF and Theoretical CDF of Response Y_i

# Empirical CDF : Version-1 using ecdf
emp.cdf.q3 <- ecdf(x = XYdata.2[,2])
emp.cdf.q3
summary(emp.cdf.q3)
str(emp.cdf.q3)
emp.cdf.q3(XYdata.2[,2]) # emprirical cdf
plot(emp.cdf.q3(XYdata.2[,2]),verticals=T)
plot.ecdf(x = XYdata.2[,2])
plot.ecdf(x = XYdata.2[,2],verticals = T)
plot(emp.cdf.q3,verticals=T)
summary.stepfun(emp.cdf.q3)
knots(Fn = emp.cdf.q3)

# emprical distribution: Version-2 Manually

sum.ecdf.q3 <- c()
Q3.Y.unique <- unique(x = XYdata.2[,2])
#Q3.Y.unique <- sort(x = XYdata.2[,2])
Q3.Y.unique
for (i5 in 1:length(Q3.Y.unique)){
  sum.q3 <- 0
  for (j5 in 1:nrow(XYdata.2)){
    sum.q3 <- ifelse(XYdata.2[j5,2]<=Q3.Y.unique[i5],sum.q3 +1,sum.q3+0) 
    
  }
  sum.ecdf.q3[i5] <- sum.q3/nrow(XYdata.2)
}
sum.ecdf.q3
plot(y = sum.ecdf.q3[order(Q3.Y.unique)],x = Q3.Y.unique,type='l')
plot(y = sum.ecdf.q3,x = Q3.Y.unique,type='l')
identical(x = sum.ecdf.q3,y = emp.cdf.q3(XYdata.2[,2]))


# emprical distribution: Version-2 Manually (Alternative)

sum.ecdf.q3.1 <- c()
#Q3.Y.unique <- unique(x = XYdata.2[,2])
Q3.Y.unique.1 <- sort(x = XYdata.2[,2])
Q3.Y.unique.1
q3.len <- length(Q3.Y.unique.1)

for (i5 in 1:q3.len){
  sum.ecdf.q3.1[i5] <- sum(Q3.Y.unique.1<=XYdata.2[i5,2])/q3.len
  }
sum.ecdf.q3.1


plot(y = sum.ecdf.q3.1,x = (XYdata.2[,2]),type='l')
#plot(y = sum.ecdf.q3,x = Q3.Y.unique,type='l')
identical(x = sum.ecdf.q3.1,y = emp.cdf.q3(XYdata.2[,2]))


# Theoretical Distribution Function.

Q3.Y.unique <- unique(x = XYdata.2[,2])
Q3.Y.unique
theor.dist.q3 <- pnorm(q = Q3.Y.unique,mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2))
theor.dist.q3
plot(theor.dist.q3)
plot(theor.dist.q3[order(Q3.Y.unique)])

# Computation of K-S Statistic Manually : Version-1

d.plus <- max(sum.ecdf.q3-theor.dist.q3)
d.plus
d.minus <- max(theor.dist.q3-sum.ecdf.q3)
d.minus

k.s.val.q3 <- max(d.plus,d.minus)
k.s.val.q3

ks.test(x = XYdata.2[,2],y = "pnorm")
(ks.test(x = XYdata.2[,2],y = "pnorm")$statistic)
ks.test(x = XYdata.2[,2],y = "pnorm")$method
ks.test(x = XYdata.2[,2],y = proj2.data.q3$ks)

prod(theor.dist.q3) # Joint CDF of Y (Response)
prod(gen.res) # Joint CDF of Y (Response)

prod(dnorm(x = XYdata.2[,2],mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)))
q3.jointden.f <- function(X,mean,var){
  n <- length(X)
  const.term <- (2*pi*var)^(-n/2)
  variable.term <- exp((-1/(2*var))*sum((X-mean)^2))
  joint.density <- const.term*variable.term
  return (joint.density)
}

q3.jointden.f(X = XYdata.2[,2],mean = mu.i.GLS.q3,var = mod.GLS.Part_1.q3$sigma.2)

q3.integrand <- function(x){
  ans <- exp((-1/2)*(x^2))
  return (ans)
}

ind.cdf.q3 <- c()
for (i in 1:nrow(XYdata.2)){
  ind.cdf.q3[i] <- (1/sqrt(2*pi))*integrate(f = q3.integrand,lower = -Inf,upper = ((XYdata.2[,2]-mu.i.GLS.q3)/sqrt(mod.GLS.Part_1.q3$sigma.2))[i])$value
}
ind.cdf.q3 # Manual Computation of CDF of Y_i

cbind(manual.norm.cdf=ind.cdf.q3,pnorm.cdf=gen.res,pnorm.cdf1=theor.dist.q3) # comparing the Manual and pnorm CDFs

prod(ind.cdf.q3) # Manual computation of Joint CDF
prod(theor.dist.q3) # Joint CDF of Y (Response) using pnorm
prod(gen.res) # Joint CDF of Y (Response) using pnorm

max(max(sum.ecdf.q3-prod(theor.dist.q3)),max(prod(theor.dist.q3)-sum.ecdf.q3))


### Comparing Empirical Dist and Theoretical CDF

plot(x = XYdata.2[,2],y = theor.dist.q3[order(XYdata.2[,2])])
plot(theor.dist.q3)
plot(sort(theor.dist.q3))
lines((sort(sum.ecdf.q3.1)))
cbind(XYdata.2[,2],theor.dist.q3)


#### Computation of K-S Statistic Manually : Version-2

max(abs(ind.cdf.q3-sum.ecdf.q3.1))
max(abs(prod(theor.dist.q3)-sum.ecdf.q3.1))
max(abs(prod(ind.cdf.q3)-sum.ecdf.q3.1))


#### Computation of K-S Statistic: Version-3

ks.test(x = XYdata.2[,2],y = mu.i.GLS.q3)



# 3.d) Complete the model assessment and explain what it means (that is, interpret the results)

# Simulation based assessment based on K-S Statistic computed between Empirical CDF and Theoretical CDF of Response Y_i.

n.q3 <- nrow(proj2.data.q3)
n.q3

p.ks <- sum(k.s.val.q3>=proj2.data.q3$ks)/n.q3
p.ks

round(sum((-0.04)>=proj2.data.q3$ks)/n.q3,10)
sum((1)>=proj2.data.q3$ks)/n.q3


summary(proj2.data.q3$ks)

sum((ks.test(x = XYdata.2[,2],y = "pnorm")$statistic)>=proj2.data.q3$ks)

sum((ks.test(x = XYdata.2[,2],y = "pnorm")$statistic)<=proj2.data.q3$ks)


# Procedure for K-S Test Statistic # Comparing Empirical CDF and Theoretical CDF of Raw Residuals w_i or Generalized Residuals Fi(Y_i).

##### Updated K-S Statistic version

ks.res.val.1 <- ks.test(x = mod.GLS.Part_1.resid,y = "pnorm") # Comparing Raw Residuals with pnorm
ks.res.val.2 <- ks.test(x = gen.res,y = "punif") # Comparing Generalized Residuals with Uniform (0,1)
ks.res.val.1
as.numeric(ks.res.val.1$statistic)
ks.res.val.2
as.numeric(ks.res.val.2$statistic)

sum(as.numeric(ks.res.val.1$statistic)>=proj2.data.q3$ks)/n.q3 # As in Course Notes, for p-value = 0.1085714
sum(as.numeric(ks.res.val.2$statistic)>=proj2.data.q3$ks)/n.q3 # As in Course Notes, for p-value = 0.09714286


sum(as.numeric(ks.res.val.1$statistic)<=proj2.data.q3$ks)/n.q3 # As in Course Notes, for p-value = 0.8914286
sum(as.numeric(ks.res.val.2$statistic)<=proj2.data.q3$ks)/n.q3 # As in Course Notes, for p-value = 0.9028571


ks.res.val.3 <- ks.test(x = gen.res.1,y = "punif")
ks.res.val.3
sum(as.numeric(ks.res.val.2$statistic)>=proj2.data.q3$ks)/n.q3 # As in Course Notes, for p-value = 0.09714286







### Simulation based assessment based on AR(1) Lag-1 Partial Auto-correlation values.

summary(proj2.data.q3$rho)

# Calculation of Lag-1 Partial Autocorrelation of Raw Residuals W_i

acf(x = mod.GLS.Part_1.resid,lag.max = 1,type = "partial",plot = F)
pacf(x = mod.GLS.Part_1.resid,lag.max = 1,plot = F)
#unlist(pacf(x = mod.GLS.Part_1.resid,lag.max = 1,plot = F))
pacf.GLS.res.q3 <- (pacf(x = mod.GLS.Part_1.resid,lag.max = 1,plot = F))$acf
as.numeric(pacf.GLS.res.q3)

pacf(x = mod.GLS.Part_1.resid,plot = F)
pacf(x = mod.GLS.Part_1.resid,plot = T)

pacf(x = pnorm(q = XYdata.2[,2],mean = mu.i.GLS.q3,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)),lag.max = 1,plot=F)


# Calculation of Lag-1 Partial Autocorrelation of Generalized Residuals r_i


pacf.gen.res.q3 <- (pacf(x = gen.res,lag.max = 1,plot = F))$acf
as.numeric(pacf.gen.res.q3)




## Calculation of p-value of Discrepancy Measure of PACF(1)


### Using Raw Residuals (W_i)
p.acf <- sum(as.numeric(pacf.GLS.res.q3)>=proj2.data.q3$rho)/n.q3
p.acf # P-value = 1

p.acf.1 <- sum(as.numeric(pacf.GLS.res.q3)<=proj2.data.q3$rho)/n.q3
p.acf.1 # P-value = 0

min(p.acf,p.acf.1)  # p-value = 0


### Using Generalized Residuals (r_i)
p.acf.2 <- sum(as.numeric(pacf.gen.res.q3)>=proj2.data.q3$rho)/n.q3
p.acf.2 # P-value = 1

p.acf.3 <- sum(as.numeric(pacf.gen.res.q3)<=proj2.data.q3$rho)/n.q3
p.acf.3 # P-value = 0

min(p.acf.2,p.acf.3) # p-value = 0

#######################################

rnorm(n = 25,mean = 0,sd = sqrt(mod.GLS.Part_1.q3$sigma.2))
pacf(x = rnorm(n = 25,mean = 0,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)),lag.max = 1,plot = F)

sum(as.numeric(pacf(x = rnorm(n = 25,mean = 0,sd = sqrt(mod.GLS.Part_1.q3$sigma.2)),lag.max = 1,plot = F)$acf)>=proj2.data.q3$rho)/n.q3


# Model-assessment based on Lag-1 PACF of AR(1) Model of Raw Residuals # So, PACF on epsilon_i of Original Model

arima(x = mod.GLS.Part_1.resid,order = c(1,0,0))
str(arima(x = mod.GLS.Part_1.resid,order = c(1,0,0)))
arima(x = mod.GLS.Part_1.resid,order = c(1,0,0))$residuals
mod.GLS.Part_1.resid
acf(x = arima(x = mod.GLS.Part_1.resid,order = c(1,0,0))$residuals,lag.max = 1,plot = F,type = "partial")
pacf(x = arima(x = mod.GLS.Part_1.resid,order = c(1,0,0))$residuals,lag.max = 1,plot = F)
str(pacf(x = arima(x = mod.GLS.Part_1.resid,order = c(1,0,0))$residuals,lag.max = 1,plot = F))
sum(as.numeric(pacf(x = arima(x = mod.GLS.Part_1.resid,order = c(1,0,0))$residuals,lag.max = 1,plot = F)$acf)>=proj2.data.q3$rho)/n.q3 #p-value = 0.8971429
sum(as.numeric(pacf(x = arima(x = mod.GLS.Part_1.resid,order = c(1,0,0))$residuals,lag.max = 1,plot = F)$acf)<=proj2.data.q3$rho)/n.q3 #p-value = 0.1028571


save.image("601_Project2.RData")





