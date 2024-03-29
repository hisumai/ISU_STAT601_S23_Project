
### STAT-601 Spring-2023-Project-Part-1 by SAMIPAN MAJUMDER
### Date: 11-April-2023.

setwd("D:/Dell_Laptop/Desktop/ISU US 19aug_or_ 25Aug2021and09Oct2021 onward/STAT 601 Advanced Statistical Methods S2022/Spring2023/Project_Part-1/")
getwd()


### Q.1.

## Load Dataset

s601.proj1.data <- read.table(file = "projdat.txt",header = T,fill = T)
str(s601.proj1.data)
print(s601.proj1.data)
nrow(s601.proj1.data)
plot(x = s601.proj1.data$x,y = s601.proj1.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")

# 1.a)


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

## Loading XYdata

XYdata.1 <- cbind(s601.proj1.data$x,s601.proj1.data$y)
XYdata.1
nrow(XYdata.1)

# Initial Values of Parameters - (alpha, beta, gamma)

# Initial values : Case-1

Q4.pars.1 <- as.numeric(c(alpha=1,beta=2,gamma=3))
Q4.pars.1

#Q4.GLS.1(par = Q4.pars.1,XYdata = XYdata.1)



X.1 <- XYdata.1[,1]
Y.1 <- XYdata.1[,2]
alpha.init.1 <- Q4.pars.1[1]
beta.init.1 <- Q4.pars.1[2]
gamma.init.1 <- Q4.pars.1[3]
mu.i.init.1 <- alpha.init.1*X.1^(beta.init.1*X.1^(-gamma.init.1))
dmu_i.dalpha.init.1 <- X.1^(beta.init.1*X.1^(-gamma.init.1))
dmu_i.dbeta.init.1 <-  mu.i.init.1*(X.1^(-gamma.init.1))*log(X.1,base = exp(1))
dmu_i.dgamma.init.1 <- mu.i.init.1*log(mu.i.init.1/alpha.init.1,base=exp(1))*log(X.1^(-1),base=exp(1))
V.j.init.1 <- matrix(data = c(dmu_i.dalpha.init.1,dmu_i.dbeta.init.1,dmu_i.dgamma.init.1),nrow =length(X.1),ncol = 3)
V.j.init.1

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

wts.func(X = XYdata.1[,1],par = Q4.pars.1)
der.func(X = XYdata.1[,1],par = Q4.pars.1)
mu.i.func(X = XYdata.1[,1],par = Q4.pars.1)

#nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=Q4.pars.1, fctn=ps[1]*xmat^(ps[2]*xmat^(-ps[3])), ders=V.j.init.1, wts=diag(length(xmat)))

#gauss.newton(xmat=XYdata.1[,1], y=XYdata.1[,2], ps=Q4.pars.1, fctn=ps[1]*xmat^(ps[2]*xmat^(-ps[3])), ders=V.j.init.1, wts=diag(length(xmat)))

#nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=Q4.pars.1, fctn=mu.i.func, ders=der.func, wts=wts.func)


# Initial values : Case-2

#Q4.GLS.1(par = c(1.5,3.5,0),XYdata = XYdata.1)

#nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=c(1.5,3.5,0), fctn=mu.i.func, ders=der.func, wts=wts.func)


# Initial values : Case-3

gamma.init.2 <- 0

mod.init.1 <- lm(log(Y.1,base=exp(1)) ~ log(X.1,base=exp(1)) )
summary(mod.init.1)
mod.init.1.coeff <- as.numeric(mod.init.1$coefficients)
mod.init.1.coeff

beta.init.2 <- mod.init.1.coeff[2]
beta.init.2

alpha.init.2 <- exp(mod.init.1.coeff[1])
alpha.init.2

Q4.pars.3 <- as.numeric(c(alpha=alpha.init.2,beta=beta.init.2,gamma=gamma.init.2))
Q4.pars.3

#Q4.GLS.1(par = Q4.pars.3,XYdata = XYdata.1)

#nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=Q4.pars.3, fctn=mu.i.func, ders=der.func, wts=wts.func)



# Initial values : Case-4

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
         mod.GLS <- Q4.GLS.1(par = pars.4,XYdata = XYdata.1)
         mod.nonlin <- nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=pars.4, fctn=mu.i.func, ders=der.func, wts=wts.func)
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


Q4.GLS.1(par = initvals[[2]],XYdata = XYdata.1)
nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=initvals[[2]], fctn=mu.i.func, ders=der.func, wts=wts.func)

#res.data <- data.frame()
alphaGLS <- c()
betaGLS <- c()
gammaGLS <- c()
varGLS <- c()
alphaNL <- c()
betaNL <- c()
gammaNL <- c()
varNL <- c()

for (i in 1:length(initvals)){
  cat("Initial Parameter values:",initvals[[i]],fill=T)
  mod.GLS.1 <- Q4.GLS.1(par = initvals[[i]],XYdata = XYdata.1)
  alphaGLS[i] <- mod.GLS.1$pars[1]
  betaGLS[i] <- mod.GLS.1$pars[2]
  gammaGLS[i] <- mod.GLS.1$pars[3]
  varGLS[i] <- mod.GLS.1$sigma.2
  mod.nonlin.1 <- nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=initvals[[i]], fctn=mu.i.func, ders=der.func, wts=wts.func)
  alphaNL <- mod.nonlin.1$bs[1]
  betaNL <- mod.nonlin.1$bs[2]
  gammaNL <- mod.nonlin.1$bs[3]
  varNL <- mod.nonlin.1$sshat[1]
}
init.res.data <- cbind(alphaGLS,betaGLS,gammaGLS,varGLS,alphaNL,betaNL,gammaNL,varNL)
init.res.data

### Result so far - version-1: 
# Parameters of Systematic component are alpha, beta, gamma. We choose a sequence of values for each parameter from 0 to 2, with length of such a sequence of values being 10. Now, we brute-force through all such combinations of initial values, and apply the GLS function and Non-Linear function, and found that 66 of such combinations are valid.

# Among the 66 combinations of initial values of parameters, all of them converged to the same final parameter estimates, and all such combinations of initial values, arrived at same value of error variance (sigma.square), using both the GLS and Non-Linear function. What may differ, between the two method functions, is the computation of Variance-Covariance Matrix, for each combination of initial parameter values.

### Parameter estimates via GLS method

# alpha= 0.506264, beta=  1.221994, gamma = 0.1613655

# Moment-based Error sigma.square variance = 1.014938













# 1.b)

mod.GLS.2 <- Q4.GLS.1(par = initvals[[2]],XYdata = XYdata.1)
mod.GLS.2$pars
mod.GLS.2$sigma.2
mod.GLS.2$sigma.Matrix
mod.GLS.2$vcov.Mat

mod.nonlin.2 <- nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=initvals[[2]], fctn=mu.i.func, ders=der.func, wts=wts.func)
mod.nonlin.2$bs
mod.nonlin.2$sshat
mod.nonlin.2$covb

mod.GLS.3 <- Q4.GLS.1(par = initvals[[5]],XYdata = XYdata.1)
mod.GLS.3$pars
mod.GLS.3$sigma.2
mod.GLS.3$sigma.Matrix
mod.GLS.3$vcov.Mat

mod.nonlin.3 <- nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=initvals[[5]], fctn=mu.i.func, ders=der.func, wts=wts.func)
mod.nonlin.3$bs
mod.nonlin.3$sshat
mod.nonlin.3$covb

Q4.GLS.1(par = initvals[[1]],XYdata = XYdata.1)$vcov.Mat

nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=initvals[[1]], fctn=mu.i.func, ders=der.func, wts=wts.func)$covb

### Approximate 95% CI for parameters in Systematic component

mod.GLS.4 <- Q4.GLS.1(par = initvals[[1]],XYdata = XYdata.1)
mod.GLS.4$pars
mod.GLS.4$sigma.2 # MOM-variance = 1.014938
mod.GLS.4$vcov.Mat
round(mod.GLS.4$sigma.2,4) # 1.0149
sqrt(mod.GLS.4$sigma.2)
round(sqrt(mod.GLS.4$sigma.2),4)

mod.nonlin.4 <- nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=initvals[[1]], fctn=mu.i.func, ders=der.func, wts=wts.func)
mod.nonlin.4$bs
mod.nonlin.4$sshat
mod.nonlin.4$covb

## For alpha

alpha.est <- mod.GLS.4$pars[1]
  alpha.LowerCI <- alpha.est - 1.96*sqrt(mod.GLS.4$vcov.Mat[1,1]/nrow(XYdata.1))
  alpha.UpperCI <- alpha.est + 1.96*sqrt(mod.GLS.4$vcov.Mat[1,1]/nrow(XYdata.1))


  ## For beta
  
  beta.est <- mod.GLS.4$pars[2]
  beta.LowerCI <- beta.est - 1.96*sqrt(mod.GLS.4$vcov.Mat[2,2]/nrow(XYdata.1))
  beta.UpperCI <- beta.est + 1.96*sqrt(mod.GLS.4$vcov.Mat[2,2]/nrow(XYdata.1))
  
  
  ## For gamma
  
  gamma.est <- mod.GLS.4$pars[3]
  gamma.LowerCI <- gamma.est - 1.96*sqrt(mod.GLS.4$vcov.Mat[3,3]/nrow(XYdata.1))
  gamma.UpperCI <- gamma.est + 1.96*sqrt(mod.GLS.4$vcov.Mat[3,3]/nrow(XYdata.1))


  CI.Int.1 <- cbind(Alpha=as.numeric(c(alpha.est,alpha.LowerCI,alpha.UpperCI)),Beta=as.numeric(c(beta.est,beta.LowerCI,beta.UpperCI)),Gamma=c(gamma.est,gamma.LowerCI,gamma.UpperCI))
  row.names(CI.Int.1) <- c("Point Estimate","Lower_95%_CI","Upper_95% CI")
  CI.Int.1
  round(CI.Int.1,4)
  
  

# 1.c)

  # Fitted Expectation Function
  
  mu.i.est <- alpha.est*XYdata.1[,1]^(beta.est*XYdata.1[,1]^(-gamma.est))
  mu.i.est
  
  # Standard error of Expectation function
  
  mod.GLS.4 <- Q4.GLS.1(par = initvals[[1]],XYdata = XYdata.1)
  mod.GLS.4$pars
  mod.GLS.4$sigma.2
  mod.GLS.4$vcov.Mat
  mod.GLS.4$Derivatives.Mat
  
  nonlin(xmat=XYdata.1[,1], ys=XYdata.1[,2], ps=initvals[[1]], fctn=mu.i.func, ders=der.func, wts=wts.func)
  
  der.mat <- mod.GLS.4$Derivatives.Mat
  vcov.mat <- mod.GLS.4$vcov.Mat
  der.mat%*%vcov.mat%*%t(der.mat)
  
  mu.i.var <- diag(der.mat%*%vcov.mat%*%t(der.mat))
  mu.i.var
  
  mu.i.se <- sqrt(mu.i.var)
  mu.i.se
  
  mu.i.LowerCI <- mu.i.est - 1.96*mu.i.se
  mu.i.UpperCI <- mu.i.est + 1.96*mu.i.se
  
  cbind(mu.i.LowerCI,mu.i.est,mu.i.UpperCI)
  
  windows(width = 10,height = 8)
  plot(x = s601.proj1.data$x,y = s601.proj1.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")
  lines(x = s601.proj1.data$x[order(s601.proj1.data$x)],y = mu.i.est[order(s601.proj1.data$x)],col="red",lwd=2)
  lines(x = s601.proj1.data$x[order(s601.proj1.data$x)],y =  mu.i.LowerCI[order(s601.proj1.data$x)],col="green",lty=2,lwd=1.5)
  lines(x = s601.proj1.data$x[order(s601.proj1.data$x)],y = mu.i.UpperCI[order(s601.proj1.data$x)],col="green",lty=2,lwd=1.5)
  dev.off()
  
  #### Second version of standard error of Mean Expectation function
  
  mu.i.se.1 <- sqrt(mu.i.var/nrow(XYdata.1))
  mu.i.se.1
  
  mu.i.LowerCI.1 <- mu.i.est - 1.96*mu.i.se.1
  mu.i.UpperCI.1 <- mu.i.est + 1.96*mu.i.se.1
  
  cbind(mu.i.LowerCI.1,mu.i.est,mu.i.UpperCI.1)
  
  windows(width = 10,height = 8)
  plot(x = s601.proj1.data$x,y = s601.proj1.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")
  lines(x = s601.proj1.data$x[order(s601.proj1.data$x)],y = mu.i.est[order(s601.proj1.data$x)],col="red",lwd=2)
  lines(x = s601.proj1.data$x[order(s601.proj1.data$x)],y =  mu.i.LowerCI.1[order(s601.proj1.data$x)],col="green",lty=2,lwd=1.5)
  lines(x = s601.proj1.data$x[order(s601.proj1.data$x)],y = mu.i.UpperCI.1[order(s601.proj1.data$x)],col="green",lty=2,lwd=1.5)
  dev.off()
  
  
  
# 1.d)

# Raw Residuals
  
  cbind(mu.i.est,mod.nonlin.4$yhat)
  res.raw <- s601.proj1.data$y - mu.i.est
  res.raw
  
  
# Standardized Residuals
  
  V.mat <- der.mat
  W.mat <- diag(nrow(XYdata.1))
  H.mat <- (W.mat^(0.5))%*%V.mat%*%solve(t(V.mat)%*%W.mat%*%V.mat)%*%t(V.mat)%*%(W.mat^(0.5))
  H.mat
  
  sigma.sq <- mod.GLS.4$sigma.2
  g2.func <- 1
  h <- diag(H.mat)
  h
  
  res.std <- res.raw/sqrt(sigma.sq*g2.func*(1-h))
  res.std
  mean(res.std) # 0.009415311
  round(mean(res.std),1) # 0
  
  cbind(res.std,mod.nonlin.4$stdres)
  
  windows(width = 10,height = 8)
  plot(x = s601.proj1.data$x,y = res.std,type='b',xlab = "Covariate(X):Phosphorus content of Soil",ylab="Standardized Residuals",col="red")
  abline(h = round(mean(res.std),1))
  dev.off()
  
  
#### Q.2.)


# 2.d)

  w.i <- res.raw
  w.i
  n.1 <- length(w.i)
  n.1
  ## Naive function
  
  theta.naive <- sum(w.i)/sum(w.i[1:(n.1-1)])
  theta.naive
  round(theta.naive,4)
  
  ## Optimal function
  
  theta.optimal <- sum(w.i[2:n.1]*w.i[1:(n.1-1)])/sum(w.i[1:(n.1-1)]^2)
  theta.optimal
  round(theta.optimal,4)
  
  # sigma^2
  
  ## Version-1
  
  sigma.sq.1.1 <- ((1-theta.naive^2)*sum(w.i^2))/n.1 # Naive
  sigma.sq.1.1
  
  sigma.sq.1.2 <- ((1-theta.optimal^2)*sum(w.i^2))/n.1 # Optimal
  sigma.sq.1.2
  
  ## Version-2
  
  sigma.sq.2.1 <- (sum((w.i[2:n.1]-theta.naive*w.i[1:(n.1-1)])^2))/n.1 # Naive
  sigma.sq.2.1
  
  sigma.sq.2.2 <- (sum((w.i[2:n.1]-theta.optimal*w.i[1:(n.1-1)])^2))/n.1 # Optimal
  sigma.sq.2.2
  
  ## Version-3
  
  sigma.sq.3.1 <- (((1-theta.naive^2)*sum(w.i[1:(n.1-1)]^2))/(n.1-1))*(((sum((w.i[2:n.1]-theta.naive*w.i[1:(n.1-1)])^2))/n.1)/((sum((w.i[2:n.1]-theta.naive*w.i[1:(n.1-1)])^2))/n.1)) # Naive
  sigma.sq.3.1
  
  sigma.sq.3.2 <- (((1-theta.optimal^2)*sum(w.i[1:(n.1-1)]^2))/(n.1-1))*(((sum((w.i[2:n.1]-theta.optimal*w.i[1:(n.1-1)])^2))/n.1)/((sum((w.i[2:n.1]-theta.optimal*w.i[1:(n.1-1)])^2))/n.1)) # Optimal
  sigma.sq.3.2
  
  Var.compare <- cbind(Naive=as.numeric(c(sigma.sq.1.1,sigma.sq.2.1,sigma.sq.3.1)),Optimal=as.numeric(c(sigma.sq.1.2,sigma.sq.2.2,sigma.sq.3.2)))
  row.names(Var.compare) <- c("Var.est.v1","Var.est.v2","Var.est.v3")
  Var.compare
  round(Var.compare,4)
  
  var.final <- c(GLS.var=sigma.sq,Naive.var=sigma.sq.2.1,Optimal.var=sigma.sq.2.2)
  var.final
  round(var.final,4)
  
  

# 2.e)

#M <- 2500
#alpha <- 0.05


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

boot.sample(theta = theta.naive,var = sigma.sq.2.1,n = nrow(XYdata.1),type = "naive")

boot.sample(theta = theta.naive,var = sigma.sq.2.1,n = nrow(XYdata.1),type = "optimal") #Using naive theta and variance

boot.sample(theta = theta.naive,var = sigma.sq.2.1,n = nrow(XYdata.1),type = "")

boot.sample(theta = theta.naive,var = sigma.sq.2.1,n = nrow(XYdata.1),type = "skhdfks")

boot.sample(theta = theta.optimal,var = sigma.sq.2.2,n = nrow(XYdata.1),type = "optimal")

boot.sample(theta = theta.naive,var = sigma.sq.2.1,n = nrow(XYdata.1),type = "naive",M=3000,alpha = 0.10)



#   for (i in 1:M){
#     w.i.boot <- c()
#     e.i.boot <- c()
#     w.0 <- rnorm(n = 1,mean = 0,sd = sqrt(sigma.sq.2.1/(1-theta.naive^2)))
#     e.i.boot[1] <- rnorm(n = 1,mean = 0,sd = sqrt(sigma.sq.2.1))
#     w.i.boot[1] <- theta.naive*w.0 + e.i.boot[1]
#     
#     for (j in 2:n.1){
#       e.i.boot[j] <- rnorm(n = 1,mean = 0,sd = sqrt(sigma.sq.2.1))
#       w.i.boot[j] <- theta.naive*w.i.boot[j-1] + e.i.boot[j]
#     }
# 
#   theta.star[i] <- sum(w.i.boot)/sum(w.i.boot[1:(n.1-1)])
#   sigma.sq.star[i] <- (sum((w.i.boot[2:n.1]-theta.star[i]*w.i.boot[1:(n.1-1)])^2))/n.1
# }
# theta.star
# sigma.sq.star
# 
# bias.theta <- sum(theta.star-theta.naive)/M
# bias.theta
# 
# theta.star[order(theta.star)][floor((M+1)*(alpha/2))]
# quantile(x = theta.star,probs = c(0.025,0.975))
# 
# bias.sigma.sq <- sum(sigma.sq.star-sigma.sq.2.1)/M
# bias.sigma.sq



save.image("601_Project1.RData")

