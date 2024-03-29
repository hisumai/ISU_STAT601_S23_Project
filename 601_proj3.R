### STAT-601 Spring-2023-Project-Part-3 by SAMIPAN MAJUMDER
### Date: 1-May-2023.

### When sourcing, it may result in an error, if some objects don't exist in the current environment.

# TO DO: HPD intervals; another combination of proposed variances.


#tryCatch(expr = {withCallingHandlers(expr = {

setwd("D:/Dell_Laptop/Desktop/ISU US 19aug_or_ 25Aug2021and09Oct2021 onward/STAT 601 Advanced Statistical Methods S2022/Spring2023/Project-Part-3/")
getwd()




## Load Dataset

s601.proj3.data <- read.table(file = "projdat.txt",header = T,fill = T)
str(s601.proj3.data)
print(s601.proj3.data)
nrow(s601.proj3.data)
plot(x = s601.proj3.data$x,y = s601.proj3.data$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption")

## Loading XYdata - Part-3 Project

XYdata.3 <- s601.proj3.data
str(XYdata.3)
print(XYdata.3)



####################################################
## Functions/Methods from Part-1 and Part-2 Project

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
          mod.GLS <- Q4.GLS.1(par = pars.4,XYdata = XYdata.3)
          mod.nonlin <- nonlin(xmat=XYdata.3[,1], ys=XYdata.3[,2], ps=pars.4, fctn=mu.i.func, ders=der.func, wts=wts.func)
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


"CLE.theta._i" <- function(X,Y,par,maxf){
  est <- nlm(f = maxf,p = par,hessian = T,X=X,Y=Y)$estimate
  return (est)
}

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


q3.jointden.f <- function(X,mean,var){
  n <- length(X)
  const.term <- (2*pi*var)^(-n/2)
  variable.term <- exp((-1/(2*var))*sum((X-mean)^2))
  joint.density <- const.term*variable.term
  return (joint.density)
}

q3.integrand <- function(x){
  ans <- exp((-1/2)*(x^2))
  return (ans)
}


#######################################################


#-------------------------------------------------------

## Metropolis-Hastings within Gibbs Sampling Markov chain Monte Carlo (MCMC) Algorithm Functions, shared by Prof. Kaiser of (STAT-601).

freundpost<-function(pars,dat,priorpars,proposalpars,B,M){
  #overall Gibbs algorithm for extended freundlich model with AR1 errors
  #three Gibbs steps, one is Metropolis for alpha, beta
  #one is Metropolis for  gamma
  #other is Metropolis for theta, sigma
  #pars is vector of starting values alpha, beta, gamma, theta, sigma
  #dat has names $x and $y
  #priorpars is vector M,V,A,B,L,P,S
  #proposalpars is vector valp,vbet,vgam,vthet,vsig
  #B is burn-in
  #M is number of MC iterations to collect after burn-in
  ind1<-NULL; ind2<-NULL; ind3<-NULL
  alps<-NULL; bets<-NULL; gams<-NULL; thets<-NULL; sigs<-NULL
  cpars<-pars
  cnt<-0
  repeat{
    cnt<-cnt+1
    #metropolis for alpha, beta
    abstar<-abproposals(cpars,proposalpars)
    abres<-abmetrop(cpars,abstar,priorpars,proposalpars,dat)
    if(cnt>B) ind1<-c(ind1,abres[1])
    newpars1<-abres[-1]
    #metropolis for gamma
    gamstar<-gproposal(newpars1,proposalpars)
    gres<-gmetrop(newpars1,gamstar,priorpars,proposalpars,dat)
    if(cnt>B) ind2<-c(ind2,gres[1])
    newpars2<-gres[-1]
    #metropolis for theta, sigma
    tsstar<-tsproposals(newpars2,proposalpars)
    tsres<-tsmetrop(newpars2,tsstar,priorpars,proposalpars,dat)
    if(cnt>B) ind3<-c(ind3,tsres[1])
    newpars<-tsres[-1]
    if(cnt>B){
      alps<-c(alps,newpars[1])
      bets<-c(bets,newpars[2])
      gams<-c(gams,newpars[3])
      thets<-c(thets,newpars[4])
      sigs<-c(sigs,newpars[5])}
    cpars<-newpars
    if(cnt==B+M) break
  }
  res<-data.frame(ind1=ind1,ind2=ind2,ind3=ind3,alp=alps,bet=bets,gam=gams,thet=thets,sig=sigs)
  return(res)
}
#-------------------------------------------------------------------------
prioralp<-function(alp,M,V){
  #compute prior pdf for alpha
  #truncated normal (at zero)
  #
  num<-dnorm(alp,M,sqrt(V))
  den<-pnorm(alp,M,sqrt(V))
  palp<-num/den
  return(palp)
}
#------------------------------------------------
priorbet<-function(bet,A,B){
  #compute prior pdf for beta
  #gamma
  pbet<-dgamma(bet,A,B)
  return(pbet)
}
#------------------------------------------------
priorgam<-function(gam,L){
  #compute prior for gamma
  #exponential
  pgam<-dgamma(gam,1,L)
  return(pgam)
}
#---------------------------------------------
priorthet<-function(thet,P){
  #compute prior for theta
  #Laplace mean 0 and variance 2 phi^2
  #truncated at -1 and 1
  pthet<-(1/(2*P))*exp(-abs(thet)/P)
  pthet2<-function(thet,phi){(1/(2*phi))*exp(-abs(thet)/phi)}
  den<-integrate(pthet2,lower=-1,upper=1,phi=P)$value
  pthet3<-pthet/den
  return(pthet3)
}
#----------------------------------------------
priorsig<-function(sig,S){
  #compute prior for sigma
  #uniform (0,S)
  psig<-(1/S)*((sig>0) & (sig<S))
  return(psig)
}
#---------------------------------------------
datmod<-function(dat,pars){
  #compute joint density for observed responses
  #pars is vector alpha, beta, gamma, theta, sigma (not sigma2)
  #dat has names $x and $y
  #
  ys<-dat$y
  xs<-dat$x
  n<-length(ys)
  thet<-pars[4]; sig<-pars[5]
  mus<-freundlich(xs,pars)
  zs<-ys-mus
  zsm1<-c(0,zs[1:(n-1)])
  k<-(n/2)*log(2*pi*sig^2)
  eterm<-(0.5/sig^2)*sum((zs-thet*zsm1)^2)
  joint<-exp(-1*k-eterm)
  return(joint)
}
#----------------------------------------------------
freundlich<-function(xs,pars){
  #compute extended Freundlich response curve
  #pars is alpha, beta, gamma
  #xs are covariates (e.g., phosphorus concentration in soil)
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]
  epart<-bet*xs^(-gam)
  fs<-alp*xs^epart
  #res<-data.frame(x=xs,f=fs)
  res<-fs
  return(res)
}
#--------------------------------------------------------
abproposals<-function(pars,proposalpars){
  #proposal values for alpha, beta
  #proposalpars is vector valp, vbet, vgam, vthet, vsig
  #all from random walks
  #all restricted to be positive
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  valp<-proposalpars[1]; vbet<-proposalpars[2]
  zalp<-rnorm(1,0,sqrt(valp))
  zbet<-rnorm(1,0,sqrt(vbet))
  alpstar<-alp+zalp
  betstar<-bet+zbet
  if(alpstar<0) alpstar<-alp
  if(betstar<0) betstar<-bet
  res<-c(alpstar,betstar,gam,thet,sig)
  return(res)
}
#-------------------------------------------------------
gproposal<-function(pars,proposalpars){
  #proposal value for gamma
  #proposalpars is vector valp, vbet, vgam, vthet, vsig
  #from random walk
  #restricted to be positive
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  vgam<-proposalpars[3]
  zgam<-rnorm(1,0,sqrt(vgam))
  gamstar<-gam+zgam
  if(gamstar<0) gamstar<-gam
  res<-c(alp,bet,gamstar,thet,sig)
  return(res)
}
#-------------------------------------------------------
tsproposals<-function(pars,proposalpars){
  #proposal values for theta and sigma
  #both random walks
  #theta restricted to be on the interval (-1,1)
  #and sigma restricted to be positive 
  #pars is vector alpha, beta, gamma, theta, sigma
  #proposalpars is vector valp, vbet, vgam, vthet, vsig
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  vthet<-proposalpars[4]; vsig<-proposalpars[5]
  zthet<-rnorm(1,0,sqrt(vthet))
  zsig<-rnorm(1,0,sqrt(vsig))
  thetstar<-thet+zthet
  sigstar<-sig+zsig
  if((thetstar< -1) | (thetstar>1)) thetstar<-thet
  if(sigstar<0) sigstar<-sig
  res<-c(alp,bet,gam,thetstar,sigstar)
  return(res)
}
#------------------------------------------------------------
abmetrop<-function(pars,stars,priorpars,vars,dat){
  #conduct metropolis  for alpha, beta
  #pars is vector alpha, beta, gamma, theta, sigma -- current values
  #stars are jump proposals (alphastar, betastar, gamma, theta, sigma)
  #as produced by function abproposals
  #priorpars is vector M,V (for alpha), A,B (for beta), L (for gamma), PHI (for theta), S (for sigma)
  #vars are variances for proposal distributions
  #all computations are in log scale
  #then compbine into Metropolis acceptance probability for alpha, beta
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  alpstar<-stars[1]; betstar<-stars[2]
  valp<-vars[1]; vbet<-vars[2]
  M<-priorpars[1]; V<-priorpars[2]; A<-priorpars[3]; B<-priorpars[4]
  pialp<-log(prioralp(alp,M,V))
  pialpstar<-log(prioralp(alpstar,M,V))
  pibet<-log(priorbet(bet,A,B))
  pibetstar<-log(priorbet(betstar,A,B))
  fypar<-log(datmod(dat,pars))
  fystar<-log(datmod(dat,stars))
  metprob<-exp(pialpstar+pibetstar+fystar-pialp-pibet-fypar)
  accept<-min(c(1,metprob))
  ustar<-runif(1,0,1)
  if(ustar<=accept){ newpars<-stars; ind1<-1}
  if(ustar>accept){ newpars<-pars; ind1<-0}
  res<-c(ind1,newpars)
  return(res)
}
#-------------------------------------------------------
gmetrop<-function(pars,stars,priorpars,vars,dat){
  #conduct metropolis  for gamma
  #pars is vector alpha, beta, gamma, theta, sigma -- current values
  #stars are jump proposals (alpha, beta, gammastar, theta, sigma)
  #as produced by function gproposal
  #priorpars is vector M,V (for alpha), A,B (for beta), L (for gamma), PHI (for theta), S (for sigma)
  #vars are variances for proposal distributions
  #all computations are in log scale
  #then compbine into Metropolis acceptance probability for gamma
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  gamstar<-stars[3]
  vgam<-vars[3]
  L<-priorpars[5]
  pigam<-log(priorgam(gam,L))
  pigamstar<-log(priorgam(gamstar,L))
  fypar<-log(datmod(dat,pars))
  fystar<-log(datmod(dat,stars))
  qgamstar<-(-0.5/vgam)*(gamstar-gam)^2
  qgam<-(-0.5/vgam)*(gam-gamstar)^2
  metprob<-exp(pigamstar+fystar-pigam-fypar)
  accept<-min(c(1,metprob))
  ustar<-runif(1,0,1)
  if(ustar<=accept){ newpars<-stars; ind1<-1}
  if(ustar>accept){ newpars<-pars; ind1<-0}
  res<-c(ind1,newpars)
  return(res)
}
#--------------------------------------------------------------
tsmetrop<-function(pars,stars,priorpars,vars,dat){
  #compute proposal densities for theta and sigma
  #pars is vector alpha, beta, gamma, theta, sigma -- current values
  #stars are jump proposals (alpha, beta, gamma, thetastar, sigmastar)
  #with the star values produced by tsproposals
  #priorpars is vector M V (for alpha), A B (for beta), L (for gamma),  PHI (for theta), S (for sigma)
  #vars are variances for proposal distributions
  #all computations are in log scale
  #then compbine into Metropolis acceptance probability for theta and sigma
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  thetstar<-stars[1]; sigstar<-stars[2]
  vthet<-vars[4]; vsig<-vars[5]
  PHI<-priorpars[6]; S<-priorpars[7]
  pithet<-log(priorthet(thet,PHI))
  pisig<-log(priorsig(sig,S))
  pithetstar<-log(priorthet(thetstar,PHI))
  pisigstar<-log(priorsig(sigstar,S))
  fypar<-log(datmod(dat,pars))
  fystar<-log(datmod(dat,stars))
  metprob<-exp(pithetstar+pisigstar+fystar-pithet-pisig-fypar)
  accept<-min(c(1,metprob))
  ustar<-runif(1,0,1)
  if(ustar<=accept){ newpars<-stars; ind2<-1}
  if(ustar>accept){ newpars<-pars; ind2<-0}
  res<-c(ind2,newpars)
  return(res)
}
#--------------------------------------------------------------


#-------------------------------------------------------

## Possible Set of Initial Values of Parameters.


# alpha	beta	gamma	theta	var
# 0.2222222	2.0000000	0.2222222	0	0.5
# 0.2222222	2.0000000	0.2222222	0.1	0.5 # Selected
# 0.22	2.00	0.22	0.1	0.5
# 0.5063	1.2220	0.1614	0.0589	1.0149 # Selected
# 0.5063	1.2220	0.1614	0.5422	0.629
# 0.5063	1.2220	0.1614	0.0589	0.8378
# 0.2222222	2	0.2222222	0.99	0.61     # Selected




set.seed(123)

## Question-1

### Prior Distribution Parameters for MCMC Gibbs Method.

prior.pars <- c(0.5,10,0.01,0.01,15,0.35,20)  # Order M,V,A,B,L,P,S  # Mean(M),Variance(V) of Normal Prior, Alpha(A),Beta(B) of Gamma Prior, Lambda(L) of Exponential Prior, "P" from Laplacian Prior, "S" from Uniform(0,S) Prior.
prior.pars

### Default Variance Values for the Proposal Distribution.

proposal.vars <- c(1,1,1,1,1) # var(alpha), var(beta), var(gamma), var(theta), var(sigma).
proposal.vars

#### Q.1. Case-1 Initial value  when M=10,000

q1.pars.init.1 <- c(0.2222222,2.0000000,0.2222222,0.1,sqrt(0.5))
q1.pars.init.1

q1.res1 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 0,M = 10000)
q1.res1
dim(q1.res1)

#### Q.1. Case-2 Initial value when M=10,000

q1.pars.init.2 <- c(0.2222222,	2,	0.2222222,	0.99,	sqrt(0.61))
q1.pars.init.2

q1.res2 <- freundpost(pars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 0,M = 10000)
q1.res2

#### Q.1. Case-3 Initial value when M=10,000

q1.pars.init.3 <- c(0.5063,	1.2220,	0.1614,	0.0589,	sqrt(1.0149))
q1.pars.init.3

q1.res3 <- freundpost(pars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 0,M = 10000)
q1.res3


### Trace Plot # Case when M= 10,000

plot(x = 1:10000,y = q1.res1$alp,type='l',col="red")
lines(x = 1:10000,y = q1.res2$alp,col="blue")
lines(x = 1:10000,y = q1.res3$alp,col="green")

plot(x = 1:10000,y = q1.res1$ind1)

plot(x = 1:100,y = q1.res1$alp[1:100],type='l',col="red")
lines(x = 1:100,y = q1.res2$alp[1:100],col="blue")
lines(x = 1:100,y = q1.res3$alp[1:100],col="green")

plot(x = 1:1000,y = q1.res1$alp[1:1000],type='l',col="red")
lines(x = 1:1000,y = q1.res2$alp[1:1000],col="blue")
lines(x = 1:1000,y = q1.res3$alp[1:1000],col="green")



#### Q.1. Case-1 Initial value, when M=50,000

q1.pars.init.1 <- c(0.2222222,2.0000000,0.2222222,0.1,sqrt(0.5))
q1.pars.init.1

q1.res1.1 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 0,M = 50000)
q1.res1.1
dim(q1.res1.1)

#### Q.1. Case-2 Initial value when M=50,000

q1.pars.init.2 <- c(0.2222222,	2,	0.2222222,	0.99,	sqrt(0.61))
q1.pars.init.2

q1.res2.1 <- freundpost(pars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 0,M = 50000)
q1.res2.1

#### Q.1. Case-3 Initial value when M=50,000

q1.pars.init.3 <- c(0.5063,	1.2220,	0.1614,	0.0589,	sqrt(1.0149))
q1.pars.init.3

q1.res3.1 <- freundpost(pars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 0,M = 50000)
q1.res3.1


### Trace Plot # Case when M= 50,000

plot(x = 1:50000,y = q1.res1.1$alp,type='l',col="red")
lines(x = 1:50000,y = q1.res2.1$alp,col="blue")
lines(x = 1:50000,y = q1.res3.1$alp,col="green")


plot(x = 1:1000,y = q1.res1.1$alp[1:1000],type='l',col="red")
lines(x = 1:1000,y = q1.res2.1$alp[1:1000],col="blue")
lines(x = 1:1000,y = q1.res3.1$alp[1:1000],col="green")


plot(x = 1:5000,y = q1.res1.1$alp[1:5000],type='l',col="red")
lines(x = 1:5000,y = q1.res2.1$alp[1:5000],col="blue")
lines(x = 1:5000,y = q1.res3.1$alp[1:5000],col="green")


plot(x = 1:10000,y = q1.res1.1$alp[1:10000],type='l',col="red")
lines(x = 1:10000,y = q1.res2.1$alp[1:10000],col="blue")
lines(x = 1:10000,y = q1.res3.1$alp[1:10000],col="green")

plot(x = 1:20000,y = q1.res1.1$alp[1:20000],type='l',col="red")
lines(x = 1:20000,y = q1.res2.1$alp[1:20000],col="blue")
lines(x = 1:20000,y = q1.res3.1$alp[1:20000],col="green")


plot(x = 1:50000,y = q1.res1.1$alp,type='l',col="red")
lines(x = 1:50000,y = q1.res2.1$alp,col="blue")
lines(x = 1:50000,y = q1.res3.1$alp,col="green")
abline(h=median(q1.res1.1$alp),col="red")
abline(h=median(q1.res2.1$alp),col="blue")
abline(h=median(q1.res3.1$alp),col="green")

acf(x = q1.res1.1$alp)
pacf(x = q1.res1.1$alp)
acf(x = q1.res2.1$alp)
pacf(x = q1.res2.1$alp)
acf(x = q1.res3.1$alp)
pacf(x = q1.res3.1$alp)


## Trace plots for M=50000 

iter.vec <- c(100,500,1000,5000,10000,50000)
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  plot(x = 1:iter.vec[i],y = q1.res1.1$alp[1:iter.vec[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q1.res2.1$alp[1:iter.vec[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q1.res3.1$alp[1:iter.vec[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()



plot(x = 1:50000,y = q1.res1.1$bet,type='l',col="red",main = "For Beta")
lines(x = 1:50000,y = q1.res2.1$bet,col="blue")
lines(x = 1:50000,y = q1.res3.1$bet,col="green")


plot(x = 1:50000,y = q1.res1.1$gam,type='l',col="red",main = "For Gamma")
lines(x = 1:50000,y = q1.res2.1$gam,col="blue")
lines(x = 1:50000,y = q1.res3.1$gam,col="green")

plot(x = 1:50000,y = q1.res1.1$thet,type='l',col="red",main = "For Theta")
lines(x = 1:50000,y = q1.res2.1$thet,col="blue")
lines(x = 1:50000,y = q1.res3.1$thet,col="green")

plot(x = 1:50000,y = q1.res1.1$sig,type='l',col="red",main = "For Sigma")
lines(x = 1:50000,y = q1.res2.1$sig,col="blue")
lines(x = 1:50000,y = q1.res3.1$sig,col="green")


hist(x = q1.res1.1$alp)
hist(x=q1.res2.1$alp,add=T) # Overlaid to previous hist
hist(x = q1.res3.1$alp,add=T) # Overlaid to previous hist

## Plotting the Histograms of the 3 Chains M=50000, B=0
windows(width = 10,height=10)
par(mfrow=c(3,1))
hist(x = q1.res1.1$alp,main="Distribution of 1st Chain")
hist(x=q1.res2.1$alp,main="Distribution of 2nd Chain")
hist(x = q1.res3.1$alp,main="Distribution of 3rd Chain")
par(mfrow=c(1,1))
dev.off()


## Plotting the Histograms of the 3 Chains M=10000, B=0
windows(width = 10,height=10)
par(mfrow=c(3,1))
hist(x = q1.res1$alp,main="Distribution of 1st Chain")
hist(x=q1.res2$alp,main="Distribution of 2nd Chain")
hist(x = q1.res3$alp,main="Distribution of 3rd Chain")
par(mfrow=c(1,1))
dev.off()






#### Gelman-Rubin Statistic for Convergence

  post.alpha.val1 <- data.frame(cbind(Chain1=q1.res1$alp,Chain2=q1.res2$alp,Chain3=q1.res3$alp))  # M=10,000
  post.alpha.val1
  nrow(post.alpha.val1)

"gelman.rubin.func" <- function(nchains,nsamples_M,val.df){
  if (!is.data.frame(val.df)){
    cat("Val.df should be a Dataframe",fill=T,"\n")
    stop()
  } else if (nsamples_M!=nrow(val.df) || nchains!=ncol(val.df)){
    cat("nsamples_M should match with Number of Rows in val.df.","And, nchains should = Number of Columns in Val.df dataframe.",fill=T,"\n")
    stop()
  }    
  
  var.chains <- apply(X = val.df,MARGIN = 2,FUN=var)
  W <- sum(var.chains)/nchains
  Mean.Chain <- apply(X = val.df, MARGIN = 2,FUN = mean)
 # Gr.Mean <- mean(val.df[,])
  Gr.Mean <- mean(Mean.Chain)/length(Mean.Chain)
  B <- (nsamples_M/(nchains-1))*sum((Mean.Chain-Gr.Mean)^2)
  var.par <- ((nsamples_M-1)/nsamples_M)*W + (1/nsamples_M)*B
  R <- sqrt(var.par/W)
  
  
  return (list(gelman.rubin.stat=R))
}

gelman.rubin.func(nchains = 3,nsamples_M = 10000,val.df = post.alpha.val1) #Gelman-Rubin=2.235088 #M=10,000, #B=0

#gelman.rubin.func(nchains = 5,nsamples_M = 1000,val.df = initvals)

gelman.rubin.func(nchains = 2,nsamples_M = 10000,val.df = post.alpha.val1[,c(1,3)]) # 1.979898

gelman.rubin.func(nchains = 2,nsamples_M = 10000,val.df = post.alpha.val1[,c(2,3)]) # 2.028189

gelman.rubin.func(nchains = 2,nsamples_M = 10000,val.df = post.alpha.val1[,c(1,2)]) #1.998421



post.alpha.val2 <- data.frame(cbind(Chain1=q1.res1.1$alp,Chain2=q1.res2.1$alp,Chain3=q1.res3.1$alp))  # M=50,000
post.alpha.val2
nrow(post.alpha.val2)

gelman.rubin.func(nchains = 3,nsamples_M = 50000,val.df = post.alpha.val2) #Gelman-Rubin=2.145072 #M=50,000, #B=0

gelman.rubin.func(nchains = 2,nsamples_M = 50000,val.df = post.alpha.val2[,c(1,3)]) # 1.925675

gelman.rubin.func(nchains = 2,nsamples_M = 50000,val.df = post.alpha.val2[,c(2,3)]) # 1.907829

gelman.rubin.func(nchains = 2,nsamples_M = 50000,val.df = post.alpha.val2[,c(1,2)]) #1.944248


"conv.diag.func" <- function(initpars, prior.pars, data,B=0,nchains=3,proposalpars,M.iter){
  #M.iter <- c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000),seq(from=20000,to=50000,by=10000))
  
  #df.chain <- data.frame(matrix(nrow=0,ncol=nchains))
  
  if (typeof(initpars)!="list" || class(initpars)!="list"){
    initpars <- list(initpars)
  }
  if (length(initpars)!=nchains){
    cat("Number of Initial Parameter vectors must match with number of chains; OR, the Input Should be a LIST of multiple vectors",fill=T,"\n")
    stop()
  }
  
  
  df.res <- data.frame()
  
  #list.chain <- list()
  
  for (i in 1:length(M.iter)){
    df.chain <- data.frame()
    df.chain <- data.frame(matrix(nrow=M.iter[i],ncol=nchains))
    for (j in 1:nchains){
      chain <-  freundpost(pars = initpars[[j]],dat = data,priorpars = prior.pars,proposalpars = proposalpars,B = 0,M = M.iter[i])$alp
      
      #print(chain)
      #print(col.name)
      #df.chain <- ifelse(test = ncol(df.chain)>0,yes = cbind(df.chain,chain),no = chain)
     # df.chain <- cbind(df.chain,chain)
      #list.chain[[j]] <- chain
      df.chain[,j] <- chain
    }
    col.name <- paste("Chain",1:nchains,sep="")
    colnames(df.chain)<-col.name
    #df.chain <- as.data.frame(cbind(list.chain))
    #print(df.chain)
    conv.val <- gelman.rubin.func(nchains = nchains,nsamples_M = M.iter[i],val.df = df.chain)$gelman.rubin.stat
    df.res <- rbind(df.res,c(M.iter[i],conv.val))
    cat("Iter=",df.res[nrow(df.res),1],"Gel-Rub_Stat",df.res[nrow(df.res),2],fill=T,"\n")
   # df.chain <- df.chain[0,]
  }
  colnames(df.res) <- c("Num_of_Iterations","Gelman-Rubin statistic")
  
  return (df.res)
}

# Using Case-3 Initial Values
#conv.res1 <- conv.diag.func(initpars = q1.pars.init.3,data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars,B = 0,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000)))
if (exists(x = "conv.res1")){
conv.res1
min(conv.res1$`Gelman-Rubin statistic`)
which.min(conv.res1$`Gelman-Rubin statistic`)
conv.res1[which.min(conv.res1$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
#         90               1.767779
#Num_of_Iterations Gelman-Rubin statistic
#        500               1.712358
}

# Using Case-3 Initial Values
#conv.res2 <- conv.diag.func(initpars = q1.pars.init.3,data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars,B = 0,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000),seq(from=20000,to=50000,by=10000)))
if (exists(x = "conv.res2")){
conv.res2
conv.res2[which.min(conv.res2$`Gelman-Rubin statistic`),]
#Num_of_Iterations Gelman-Rubin statistic
#         90               1.767779
}

### Gelman-Rubin based on combining three vectors of initial parameters ### M till 10000, B=0
conv.res3 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars,B = 0,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000)))
conv.res3
conv.res3[which.min(conv.res3$`Gelman-Rubin statistic`),]
#Num_of_Iterations Gelman-Rubin statistic
#6                60               1.643763

### Gelman-Rubin based on combining three vectors of initial parameters ### M till 50000, B=0
#conv.res4 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars,B = 0,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000),seq(from=20000,to=50000,by=10000)))
if (exists(x = "conv.res4")){
conv.res4
conv.res4[which.min(conv.res4$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 3                30               1.526561
}

### New Version of Joint Posterior Model, with B = 90,M=50,000


#### Q.1. Case-1 Initial value when M=50,000, B= 90
q1.pars.init.1
q1.res1.2 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 90,M = 50000)
q1.res1.2
dim(q1.res1.2)

#### Q.1. Case-2 Initial value when M=50,000, B=90

q1.pars.init.2
q1.res2.2 <- freundpost(pars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 90,M = 50000)
q1.res2.2

#### Q.1. Case-3 Initial value when M=50,000, B=90

q1.pars.init.3
q1.res3.2 <- freundpost(pars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 90,M = 50000)
q1.res3.2



## Trace plots for M=50000 , B=90

iter.vec <- c(100,500,1000,5000,10000,50000)
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  plot(x = 1:iter.vec[i],y = q1.res1.2$alp[1:iter.vec[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q1.res2.2$alp[1:iter.vec[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q1.res3.2$alp[1:iter.vec[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()


### Performing separate execution of Joint Posterior Model for M= 100,500,1000,5000,10000,50000 and B= 90, for three chains.

#q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 90,M = 50000

"joint.post.func" <- function(niter,initpars,dat,priorpars,proposalpars,B){
  res.df <- data.frame()
  
  for (i in 1:length(niter)){
    res.df  <- freundpost(pars = initpars,dat = dat,priorpars = priorpars,proposalpars = proposalpars,B = B,M = niter[i])
    if (i==1){
      res.df1 <- matrix(nrow=niter[i],ncol=ncol(res.df))
      res.df1 <- res.df
    } else {
      res.df1 <- rbind(res.df1,res.df)
    }
  }
  res.df1 <- as.data.frame(res.df1)
  
  return (res.df1)
}

# Q.1. Combined Joint Posterior for Different Number of samples b= 90, Case-1 Initial pars.
#q1.res1.3 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 90)
if (exists(x = "q1.res1.3")){
q1.res1.3
dim(q1.res1.3)
}
# Q.1. Combined Joint Posterior for Different Number of samples b= 90, Case-2 Initial pars.
#q1.res2.3 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 90)
if (exists(x = "q1.res2.3")){
q1.res2.3
dim(q1.res2.3)
}
# Q.1. Combined Joint Posterior for Different Number of samples b= 90, Case-3 Initial pars.
#q1.res3.3 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 90)
if (exists(x = "q1.res3.3")){
q1.res3.3
dim(q1.res3.3)
}

## Trace plots for M=c(100,500,1000,5000,10000,50000) , B=90

if (exists(x = "q1.res1.3") && exists(x = "q1.res2.3") && exists(x = "q1.res3.3")){
iter.vec <- c(100,500,1000,5000,10000,50000)
iter.vec1 <- cumsum(iter.vec)
iter.vec1
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec1)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  #print((((i-1)==0)?(1):(iter.vec[i-1]+1)))
  print((ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)))
 # print((ifelse(test = (i-1)==0,yes = iter.vec[i],no = iter.vec[i-1]+1)))
  plot(x = 1:iter.vec[i],y = q1.res1.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q1.res2.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q1.res3.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()
}


### New Version of Joint Posterior Model, with B = 500,M=50,000


#### Q.1. Case-1 Initial value when M=50,000, B= 500
q1.pars.init.1
q1.res1.4 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 500,M = 50000)
q1.res1.4
#dim(q1.res1.2)

#### Q.1. Case-2 Initial value when M=50,000, B=500

q1.pars.init.2
q1.res2.4 <- freundpost(pars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 500,M = 50000)
q1.res2.4

#### Q.1. Case-3 Initial value when M=50,000, B=90

q1.pars.init.3
q1.res3.4 <- freundpost(pars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 500,M = 50000)
q1.res3.4


## Trace plots for M=50000 , B=500

iter.vec <- c(100,500,1000,5000,10000,50000)
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  plot(x = 1:iter.vec[i],y = q1.res1.4$alp[1:iter.vec[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q1.res2.4$alp[1:iter.vec[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q1.res3.4$alp[1:iter.vec[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()


# Q.1. Combined Joint Posterior for Different Number of samples b= 500, Case-1 Initial pars.
q1.res1.5 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 500)
q1.res1.5


# Q.1. Combined Joint Posterior for Different Number of samples b= 500, Case-2 Initial pars.
q1.res2.5 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 500)
q1.res2.5


# Q.1. Combined Joint Posterior for Different Number of samples b= 500, Case-3 Initial pars.
q1.res3.5 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 500)
q1.res3.5

## Trace plots for M=c(100,500,1000,5000,10000,50000) , B=500

iter.vec <- c(100,500,1000,5000,10000,50000)
iter.vec1 <- cumsum(iter.vec)
iter.vec1
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec1)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  #print((((i-1)==0)?(1):(iter.vec[i-1]+1)))
  print((ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)))
  # print((ifelse(test = (i-1)==0,yes = iter.vec[i],no = iter.vec[i-1]+1)))
  plot(x = 1:iter.vec[i],y = q1.res1.5$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q1.res2.5$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q1.res3.5$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()


# Gelman-Rubin stats of multiple iterations with B=500

conv.res5 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars,B = 500,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000)))
conv.res5
conv.res5[which.min(conv.res5$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 1                10               1.293614



# Gelman-Rubin stats of multiple iterations with B=5000 (Based on previous Trace Plots)

#conv.res6 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars,B = 5000,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000),seq(from=20000,to=50000,by=10000)))
if (exists(x = "conv.res6")){
conv.res6
conv.res6[which.min(conv.res6$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 8                80                1.47476
}




# Gelman-Rubin stats of multiple iterations with B=10000 (Based on previous Trace Plots)

#conv.res7 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars,B = 10000,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000),seq(from=20000,to=50000,by=10000)))
if (exists(x = "conv.res7")){
conv.res7
conv.res7[which.min(conv.res7$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 7                70               1.552158
}

####################################################
# Final Question-1 Summary from conv.diag.func onwards######

# Goal :For the subsequent analysis, we used the Gelman-Rubin Statistics, Trace Plots to determine the appropriate burn-in period (B)

# For now, we have assumed the proposed variances as (1,1,1,1,1), and the Burn-in Period (B) = 0.

#Using Case-3 Initial Values (Number of Samples till 10,000)
# Num_of_Iterations Gelman-Rubin statistic
#         90               1.767779
#Num_of_Iterations Gelman-Rubin statistic
#        500               1.712358

#Using Case-3 Initial Values (Number of Samples till 50,000)
#Num_of_Iterations Gelman-Rubin statistic
#         90               1.767779

# Gel-rub Statistic for three chains combined M=10,000,B=0
#Num_of_Iterations Gelman-Rubin statistic
#6                60               1.643763

# Gel-rub Statistic for three chains combined M=50,000,B=0
# Num_of_Iterations Gelman-Rubin statistic
# 3                30               1.526561


# Separate Joint Posterior Models for 3 Markov chains, with B = 90,M=50,000 (Trace Plots below)

## Trace plots for M=c(100,500,1000,5000,10000,50000) , B=90

if (exists(x = "q1.res1.3") && exists(x = "q1.res2.3") && exists(x = "q1.res3.3")){
iter.vec <- c(100,500,1000,5000,10000,50000)
iter.vec1 <- cumsum(iter.vec)
iter.vec1
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec1)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  #print((((i-1)==0)?(1):(iter.vec[i-1]+1)))
  print((ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)))
  # print((ifelse(test = (i-1)==0,yes = iter.vec[i],no = iter.vec[i-1]+1)))
  plot(x = 1:iter.vec[i],y = q1.res1.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q1.res2.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q1.res3.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()
}

# Separate Joint Posterior Models for 3 Markov chains, with B = 500,M=50,000 (Trace Plots below)
## Trace plots for M=c(100,500,1000,5000,10000,50000),B=500

iter.vec <- c(100,500,1000,5000,10000,50000)
iter.vec1 <- cumsum(iter.vec)
iter.vec1
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec1)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  #print((((i-1)==0)?(1):(iter.vec[i-1]+1)))
  print((ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)))
  # print((ifelse(test = (i-1)==0,yes = iter.vec[i],no = iter.vec[i-1]+1)))
  plot(x = 1:iter.vec[i],y = q1.res1.5$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q1.res2.5$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q1.res3.5$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()

# Gel-rub Statistic for three chains combined M=10,000,B=500
# Num_of_Iterations Gelman-Rubin statistic
# 1                10               1.293614

#Gel-rub Statistic for three chains combined M=50,000,B=5000
# Num_of_Iterations Gelman-Rubin statistic
# 8                80                1.47476

#Gel-rub Statistic for threechains combined M=50,000,B=10000
# Num_of_Iterations Gelman-Rubin statistic
# 7                70               1.552158


# Based on the above summaries, we initially thought of using a Burn-in Period of (B) = 90 0r 500, based on the minimum Gelman-rubin Statistic values among different number of iterations (or, samples) from the joint posterior model, but, the trace plots always hinted that proper mixing of chains have happened from 5000 to 10,000 samples onwards. And, between 5000 and 10,000, the minimum Gelman-Rubin statistic values is obtained in the case of 5000 samples. So, we have decided to stick to 5000 as the Burn-in(B) Period.

# But the Gelman-Rubin Statistic values that have been obtained so far, none of them were below 1.1 or 1.2, which is the ideal cutoff for that statistic. Part of the reason can be the choice of variance parameters for the proposed distribution. And another reason, may be the choice of Initial values of parameters, as well. But, these were the same initial values that were used for Likelihood estimation, as well.


######################################################










#####Gelman-Rubin Statistic Diagnosic #######

# If multiple Markov chains are appropriately mixing in the stationary distribution, then, Gel-man-Rubin Statistic (Rk) <1.2 or Rk<1.1.

################################################



###---------------------------

## Subsequent Initial Values of Parameters TO Be Tested

#0.2222222	2.0000000	0.2222222	0	0.5
#0.22	2.00	0.22	0.1	0.5
#0.5063	1.2220	0.1614	0.5422	0.629
#0.5063	1.2220	0.1614	0.0589	0.8378

## Not so Good Initial values
# 0	0.8	0	0.99	0.4008
# 0.4	0	0.4	0.594	1.6002
# 0.8	1.2	1.6	0.99	0.8006

## Full MLE Results from Part-2
## c(0.7467, 1.0209, 0.1568, 0.5597, 0.6250)

## CLE Results from Part-2
## c(1.1099, 0.7871, 0.1430, 0.5772, 0.6437)


df.init.pars <- data.frame(rbind(c(0.7467, 1.0209, 0.1568, 0.5597, 0.6250),c(1.1099, 0.7871, 0.1430, 0.5772, 0.6437),c(0.2222222,	2.0000000,	0.2222222,	0,	0.5),c(0.22,	2.00,	0.22,	0.1,	0.5),c(0.5063,	1.2220,	0.1614,	0.5422,	0.629),c(0.5063,	1.2220,	0.1614,	0.0589,	0.8378),c(0,	0.8,	0,	0.99,	0.4008),c(0.4,	0,	0.4,	0.594,	1.6002),c(0.8,	1.2,	1.6,	0.99,	0.8006)))
names(df.init.pars) <- c("alpha.init","beta.init","gamma.init","theta.init","sigma.init")
df.init.pars$sigma.init <- sqrt(df.init.pars$sigma.init)
df.init.pars








###-----------------------


#### To Do ##########

# Goal: Monitoring Convergence
# 
# Refer to https://econ.pages.code.wm.edu/414/notes/docs/convergence_diagnostics.html
# 
# Goal: Inference (HPD Intervals)
# 
# Refer: https://econ.pages.code.wm.edu/414/notes/docs/inference.html


##########################








## Question - 2

# Initial values for Markov chain-1 M=50,000
# 0.2222222	2.0000000	0.2222222	0.1	0.5 # q1.res1.1

# Initial values for Markov chain-2 M=50,000
# 0.2222222	2	0.2222222	0.99	0.61     # q1.res2.1 

# Initial values for Markov chain-3 M=50,000
# 0.5063	1.2220	0.1614	0.0589	1.0149 # q1.res3.1




#### Q2: Case-1 : proposal.vars=(1,1,1,1,1)

## ind1 : (alpha,beta), ind2 : (gamma), ind3: (theta, sigma)

# Markov chain-1  Case-1 Prop.var
sum(q1.res1.1$ind1==1)/50000  # 0.10606
sum(q1.res1.1$ind2==1)/50000  # 0.47432
sum(q1.res1.1$ind3==1)/50000  # 0.05316

round(100*sum(q1.res1.1$ind1==1)/50000,2)
round(100*sum(q1.res1.1$ind2==1)/50000,2)
round(100*sum(q1.res1.1$ind3==1)/50000,2)  



# Markov chain-2
sum(q1.res2.1$ind1==1)/50000  # 0.10654
sum(q1.res2.1$ind2==1)/50000  # 0.46518
sum(q1.res2.1$ind3==1)/50000  # 0.0644

# Markov chain-3
sum(q1.res3.1$ind1==1)/50000  # 0.1073
sum(q1.res3.1$ind2==1)/50000  # 0.46076
sum(q1.res3.1$ind3==1)/50000  # 0.07296



#### Q2: Case-2 : proposal.vars

proposal.vars.1 <- c(1.5,1.5,0.5,2,2) # var(alpha), var(beta), var(gamma), var(theta), var(sigma).
proposal.vars.1

#### Case-1 Initial value (q1.pars.init.1)  when M=50,000

q2.res1 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.1,B = 0,M = 50000)
q2.res1
dim(q2.res1)


# Markov chain-1 : Case-2 Prop.var
sum(q2.res1$ind1==1)/50000  # 0.10606
sum(q2.res1$ind2==1)/50000  # 0.47432
sum(q2.res1$ind3==1)/50000  # 0.05316

round(100*sum(q2.res1$ind1==1)/50000,2)
round(100*sum(q2.res1$ind2==1)/50000,2)
round(100*sum(q2.res1$ind3==1)/50000,2)

#optim()


# Markov chain-1 : Case-3 Prop.var

propvar3 <- seq(from=0.01,to=1,by=0.01)

# for (i in 1:length(propvar3)){
#   for (j in 1:length(propvar3)){
#     for (k in 1:length(propvar3)){
#       for (l in 1:length(propvar3)){
#         for (m in 1:length(propvar3)){
#           prop.var <- c(propvar3[i],propvar3[j],propvar3[k],propvar3[l],propvar3[m])
#           joint.post <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = prop.var,B = 0,M = 1000)
#           prop.ind1 <- sum(joint.post[,1])/1000
#           prop.ind2 <- sum(joint.post[,2])/1000
#           prop.ind3 <- sum(joint.post[,3])/1000
#           if (prop.ind1>=0.20 && prop.ind1<=0.40 && prop.ind2>=0.20 && prop.ind2<=0.40 && prop.ind3>=0.20 && prop.ind3<=0.40){
#             cat("Posterior Variances",prop.var,fill=T,"\n")
#           }
#         }
#       }
#     }
#   }
# }

"prop.var.func" <- function(from1=0.01,to1=2,by1=0.01, initpars, prior.pars, data,B=0,M=100,niter=1000,dur.secs=120){
  propvar <- seq(from=from1,to=to1,by=by1)
  
 # df1 <- data.frame(matrix(nrow = 0,ncol=length(initpars),dimnames = list(NULL,c(paste("Var.Par",1:length(initpars),sep = "")))))
  #df.res <- data.frame(matrix(nrow = 0,ncol=length(initpars),dimnames = list(NULL,c(paste("Var.Par",1:length(initpars),sep = "")))))
  
  df1 <- data.frame(matrix(nrow=0,ncol=length(initpars)))
 # colnames(df1) <- c(paste("Var.Par",1:length(initpars),sep = ""))
  df.res <- data.frame(matrix(nrow=0,ncol=length(initpars)))
 # colnames(df.res) <- c(paste("Var.Par",1:length(initpars),sep = ""))
  
  cnt <- 0
  start <- 0
  end <- 0
  start <- Sys.time()
  repeat{
  while (TRUE){
  par.var <- sample(x = propvar,size = length(initpars),replace = T)
  if (nrow(df1)>0){
  if (!identical(x = par.var,y = df1[1:nrow(df1),])){
    df1 <- rbind(df1, par.var)
    break
  }
  } else {
    df1 <- rbind(df1, par.var)
    break
  }
  }
  jo.post <- freundpost(pars = initpars,dat = data,priorpars = prior.pars,proposalpars = par.var,B = B,M = M)
  prop.i1 <- sum(jo.post[,1])/M
  prop.i2 <- sum(jo.post[,2])/M
  prop.i3 <- sum(jo.post[,3])/M
  if (prop.i1>=0.20 && prop.i1<=0.40 && prop.i2>=0.20 && prop.i2<=0.40 && prop.i3>=0.20 && prop.i3<=0.40){
    #par.var0 <- setNames(par.var,paste("Var.Par",1:length(initpars),sep = ""))
    df.res <- rbind(df.res,par.var)
    cat("Posterior Variances",par.var,fill=T,"\n")
  }
  cnt <- cnt+1
  
  if (cnt==niter){
    cat("Number of default iterations ",niter,"exceeded",fill=T,"\n")
    break}
  end <- Sys.time()
  cat("Total Time elpased till Iteration:",cnt,"=",difftime(time1 = end,time2 = start,units = "secs")[[1]],"seconds",fill=T,"\n")
  #if ((end-start)>=dur.secs){
    if ((difftime(time1 = end,time2 = start,units = "secs")[[1]])>=dur.secs){  
    cat("Default Time Limit (in seconds)",dur.secs,"exceeded",fill=T,"\n")
    break
  }
  }
  
  if (ncol(df1)>0){
    colnames(df1) <- c(paste("Var.Par",1:length(initpars),sep = ""))
  }
  if (ncol(df.res)>0){
    colnames(df.res) <- c(paste("Var.Par",1:length(initpars),sep = ""))
  }
  
  
  return (list(proposal.var.init=df1,Proposal.Variances.final=df.res,size.init.propvar=nrow(df1),size.final.prop.var=nrow(df.res)))
}

#prop.var.func(initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3)

#prop.var.func(initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,M = 50000) # Result=None 

#prop.var.func(initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,M = 10000) # Result=None


q2.var.res.1 <- prop.var.func(initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3)
q2.var.res.1

q2.prop.var <- q2.var.res.1$Proposal.Variances.final
q2.prop.var

# Latex Export of Proposed Variances for 100 samples.

print(xtable::xtable(x = head(q2.prop.var), type = "latex", tabular.environment="longtable"), file = "Q2_PropVar_v1.tex")

######################

# for (i in 1:nrow(q2.prop.var)){
#   cat("Posterior Variances",unlist(q2.prop.var[i,]),fill=T,"\n")
#   jnt.post <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = as.numeric(q2.prop.var[i,]),B = 0,M = 50000)
#   prop.in1 <- sum(jnt.post[,1])/50000
#   prop.in2 <- sum(jnt.post[,2])/50000
#   prop.in3 <- sum(jnt.post[,3])/50000
#   if (prop.in1>=0.20 && prop.in1<=0.40 && prop.in2>=0.20 && prop.in2<=0.40 && prop.in3>=0.20 && prop.in3<=0.40){
#     #df.res <- rbind(df.res,par.var)
#     cat("Selected Post_var",q2.prop.var[i,],fill=T,"\n")
#   }  
#   
# }   ##### No Result (None of the proposed variances qualified when M=50000)


########### Q.2. Case-4 Proposal.vars ##############

# Based on the Results of Q.1., especially since the computation of Gelman-Rubin Statistic, we found that Burn-in Period is 500, although the trace plots showed a somewhat stationary picture around 5000 iterations. We compared the Gelman-Rubin statistic values from 100 to 50000 iterations, and most of the time the values were around 2.0, with some values between 1.5 and 2.0, since every computation is based on random generation of samples. In all these computations in Q.1., we have assumed the posterior variance parameter to be all 1s.

# In addition in Q.1., based on trace plots, we also decided to use the Burn-in Period (B) as 5000 or 10,000. But, even using that Burn-in Period, and keeping the variances of Proposed distribution to "1"s, the Gelman-Rubin statistic never got lowered below 1.1 or 1.2, with multiple different iterations (or, number of joint posterior samples "M") ranging from 100 to 10,000 or 50,000.


"accept.prob.func" <- function(df.lst){
  
  if (class(df.lst)!="list" || typeof(df.lst)!="list"){
    cat("The Input has to be a List of Dataframe(s)",fill=T,"\n")
    stop()
  }
  else {
    for ( i in 1:length(df.lst)){
  if (!is.data.frame(df.lst[[i]])){
    cat("Individual Entries of the List has to be a Dataframe",fill=T,"\n")
    stop()
  }
  }
  }
  for ( i in 1:length(df.lst)){
    df0 <- data.frame()
    r <- nrow(df.lst[[i]])
    c <- ncol(df.lst[[i]])
    cat("r=",r,"c=",c,fill=T,"\n")
    corr.col <- c()
    #print(df.lst[[i]])
    df0 <- df.lst[[i]]
    for (j in 1:c){
   #   if (all((df.lst[[i]])[,j]==0 & (df.lst[[i]])[,j]==1)){
 #   if (any(df0[,j]==0 & df0[,j]==1) && all(!df0[,j]<0 & (!df0[,j]>0 && !df0[,j]<1) & !df0[,j]>1)){
      if (any(df0[,j]==0 | df0[,j]==1)){  
        corr.col <- c(corr.col,j)    
      }
    }
   # print(corr.col)
    df.lst[[i]] <- df.lst[[i]][,corr.col]
    cat("Reduced","r=",nrow(df.lst[[i]]),"c=",ncol(df.lst[[i]]),fill=T,"\n")
  }
  
  
  
  acc.res <- list()
  
  #print(df.lst)
  
  for (i in 1:length(df.lst)){
    tmp.res <- c()
    df0 <- data.frame()
    df0 <- df.lst[[i]]
    for (j in 1:ncol(df.lst[[i]])){
      #acc.prob <- sum(df.lst[[i]][,j]==1)/length(df.lst[[i]][,j])
      acc.prob <- sum(df0[,j]==1)/length(df0[,j])
      lab.txt <- paste("ind.prob",j,sep="")
      tmp.res <- c(tmp.res,acc.prob)
      names(tmp.res)[j] <- lab.txt
    }
  #  tmp.res <- c(tmp.res,total.prob=sum(df.lst[[i]][,1:ncol(df.lst[[i]])]==1)/length(df.lst[[i]][,1:ncol(df.lst[[i]])]))
    tmp.res <- c(tmp.res,total.prob=sum(df0[,1:ncol(df0)]==1)/(nrow(df0)*ncol(df0)))
    #print(tmp.res)
    
    acc.res[[i]] <- tmp.res
    #print(acc.res[[i]])
  }
  
  if (length(acc.res)>0){
    names(acc.res) <- c(paste("Df",1:length(acc.res),sep = ""))
  }
  #print(acc.res)
  return (acc.res)
}


#### Acceptance Probabilities with Original Prop var (1,1,1,1,1)

accept.prob.func(df.lst = list(q1.res1,q1.res2,q1.res3))

accept.prob.func(df.lst = list(q1.res1.1,q1.res2.1,q1.res3.1))

accept.prob.func(df.lst = list(q1.res1.2,q1.res2.2,q1.res3.2))

if (exists(x = "q1.res1.3") && exists(x = "q1.res2.3") && exists(x = "q1.res3.3")){
accept.prob.func(df.lst = list(q1.res1.3,q1.res2.3,q1.res3.3))
}

accept.prob.func(df.lst = list(q1.res1.4,q1.res2.4,q1.res3.4))

accept.prob.func(df.lst = list(q1.res1.5,q1.res2.5,q1.res3.5))

#accept.prob.func(df.lst = list(q1.res1.5,q1res2))
# accept.prob.func(df.lst = list(q1.res1.5,q1.res2))
# accept.prob.func(df.lst = list(q1.res1.5,q1.res2.1,q1.res3.5))


# So far, in Q.2., what we have seen is that, keeping the variances of proposed distribution as (1,1,1,1,1), for all the joint posterior model results with the different combinations of initial parameters, and with different number of iterations, we have noticed that the acceptance probabilities for index-1 is less than 0.20, for index-2 its greater than 0.40, and for index-3, its less than 0.20 also.


## So now, we will try to change the values of the variance parameters of the proposed distribution, and check if the acceptance probabilities of individual indices are within the accepted limits (0.20 to 0.40) or not.

# We will consider the First Markov Chain from Case-1 Initial Parameters, with M= 50,000, and with two options of B = 5000 or 10,000 (based on trace plots), since the Gelman-Rubin statistic did not go below 1.1 or 1.2.


# Initial Posterior-Variance proposal.vars=(1,1,1,1,1)

# Initial values for Markov chain-1 M=50,000
# 0.2222222	2.0000000	0.2222222	0.1	sqrt(0.5)

###### Q.2. Case-4.1 M=50,000, B= 5000

q1.pars.init.1 # Case-1 Initial values of parameters

proposal.vars # Initial Variances of Proposed Distribution

q2.res2.0 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars,B = 5000,M = 50000)
q2.res2.0
accept.prob.func(df.lst = list(q2.res2.0))
#ind.prob1  ind.prob2  ind.prob3 total.prob 
#0.1053800  0.4773400  0.0484600  0.2103933 
accept.prob.func(df.lst = list(q2.res2.0))$Df1
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.1053800  0.4773400  0.0484600  0.2103933 

windows(width = 10,height = 10)
par(mfrow=c(3,2))
plot(x = 1:50000,y = q2.res2.0$alp,type='l',col="red",main = "For alpha")
# lines(x = 1:50000,y = q1.res2.1$bet,col="blue")
# lines(x = 1:50000,y = q1.res3.1$bet,col="green")

plot(x = 1:50000,y = q2.res2.0$bet,type='l',col="red",main = "For Beta")

plot(x = 1:50000,y = q2.res2.0$gam,type='l',col="red",main = "For Gamma")

plot(x = 1:50000,y = q2.res2.0$thet,type='l',col="red",main = "For Theta")

plot(x = 1:50000,y = q2.res2.0$sig,type='l',col="red",main = "For Sigma")
par(mfrow=c(1,1))

dev.off()
diag(var(q2.res2.0[,4:8]))
var(q2.res2.0$alp)


### New Variances of Proposed distribution: Version-1
proposal.vars.2 <- c(0.40,0.5,0.03,0.2,0.15) # var(alpha), var(beta), var(gamma), var(theta), var(sigma).
proposal.vars.2

q2.res2.1 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.2,B = 5000,M = 50000)
q2.res2.1
accept.prob.func(df.lst = list(q2.res2.1)) # Ind.prob2 right
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.1193600  0.3499000  0.0859000  0.1850533
diag(var(q2.res2.1[,4:8]))


### New Variances of Proposed distribution: Version-2
proposal.vars.3 <- c(0.04,0.5,0.03,0.02,0.015) # var(alpha), var(beta), var(gamma), var(theta), var(sigma).
proposal.vars.3

q2.res2.2 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.3,B = 5000,M = 50000)
q2.res2.2
accept.prob.func(df.lst = list(q2.res2.2)) # Ind.prob2 right
#ind.prob1  ind.prob2  ind.prob3 total.prob 
#0.18136    0.34518    0.19124    0.23926 
diag(var(q2.res2.2[,4:8]))



### New Variances of Proposed distribution: Version-3
proposal.vars.4 <- c(0.04,0.05,0.03,0.002,0.015) # var(alpha), var(beta), var(gamma), var(theta), var(sigma).
proposal.vars.4

q2.res2.3 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.4,B = 5000,M = 50000)
q2.res2.3
accept.prob.func(df.lst = list(q2.res2.3)) # Ind.prob2 right
#ind.prob1  ind.prob2  ind.prob3 total.prob 
#0.20730    0.34640    0.23728    0.26366 
diag(var(q2.res2.3[,4:8]))

windows(width = 10,height = 10)
par(mfrow=c(3,2))
plot(x = 1:50000,y = q2.res2.3$alp,type='l',col="red",main = "For alpha")
# lines(x = 1:50000,y = q1.res2.1$bet,col="blue")
# lines(x = 1:50000,y = q1.res3.1$bet,col="green")
plot(x = 1:50000,y = q2.res2.3$bet,type='l',col="red",main = "For Beta")
plot(x = 1:50000,y = q2.res2.3$gam,type='l',col="red",main = "For Gamma")
plot(x = 1:50000,y = q2.res2.3$thet,type='l',col="red",main = "For Theta")
plot(x = 1:50000,y = q2.res2.3$sig,type='l',col="red",main = "For Sigma")
par(mfrow=c(1,1))
dev.off()

# Q.2. Gelman-Rubin statistic of Different Iterations, with B=5000, and Case-4, Version-3 of Proposed Variances
#q2.conv.res.1 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars.4,B = 5000,M.iter = c(seq(from=10,to=100,by=10),seq(from=200,to=500,by=100),seq(from=1000,to=10000,by=1000),seq(from=20000,to=50000,by=10000)))
if (exists(x = "q2.conv.res.1")){
q2.conv.res.1
q2.conv.res.1[which.min(q2.conv.res.1$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 15              1000               1.530241
}

# Q.2. Gelman-Rubin statistic of Different Iterations, with B=5000, and Case-4, Version-3 of Proposed Variances
#q2.conv.res.2 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars.4,B = 5000,M.iter = c(seq(from=5000,to=10000,by=1000),seq(from=20000,to=50000,by=10000)))
if (exists(x = "q2.conv.res.2")){
q2.conv.res.2
q2.conv.res.2[which.min(q2.conv.res.2$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 5              9000                1.95191
}

# Q.2. Gelman-Rubin statistic of Different Iterations, with B=10000, and Case-4, Version-3 of Proposed Variances
#q2.conv.res.3 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars.4,B = 10000,M.iter = c(seq(from=10000,to=50000,by=10000)))
if (exists(x = "q2.conv.res.3")){
q2.conv.res.3
q2.conv.res.3[which.min(q2.conv.res.3$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 2             20000                2.07693
}


# Q.2. Combined Joint Posterior for Different Number of samples b= 5000,(Q1.Case-3 Initial pars), using Case-4 version-3 Proposal vars
#q2.joint.post.res.1 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.4,B = 5000)
if (exists(x = "q2.joint.post.res.1")){
q2.joint.post.res.1
accept.prob.func(df.lst = list(q2.joint.post.res.1))
# $Df1  # For Total 66600 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.1995045  0.3531982  0.2407958  0.2644995 
accept.prob.func(df.lst = list(q2.joint.post.res.1[16601:66600,]))
# $Df1 # Last 50,000 samples
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.19462    0.34924    0.23728    0.26038 
}

# Q.2. Combined Joint Posterior for Different Number of samples b= 5000,(Q1.Case-1 Initial pars), using Case-4 version-3 Proposal vars
#q2.joint.post.res.2 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.4,B = 5000)
if (exists(x = "q2.joint.post.res.2")){
q2.joint.post.res.2
accept.prob.func(df.lst = list(q2.joint.post.res.2))
# $Df1  # For Total 66600 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.2041441  0.3581381  0.2224324  0.2615716 
accept.prob.func(df.lst = list(q2.joint.post.res.2[16601:66600,]))
# $Df1 # Last 50,000 samples
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.2075600  0.3517400  0.2293800  0.2628933 
}


# Q.2. Combined Joint Posterior for Different Number of samples b= 5000,(Q1.Case-2 Initial pars), using Case-4 version-3 Proposal vars
#q2.joint.post.res.3 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.4,B = 5000)
if (exists(x = "q2.joint.post.res.3")){
q2.joint.post.res.3
accept.prob.func(df.lst = list(q2.joint.post.res.3))
# $Df1  # For Total 66600 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.2015315  0.3596096  0.2232282  0.2614565 
accept.prob.func(df.lst = list(q2.joint.post.res.3[16601:66600,]))
# $Df1 # Last 50,000 samples
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.1977600  0.3546400  0.2323800  0.2615933 
}




## Q.2. Trace plots for M=c(100,500,1000,5000,10000,50000) , B=5000  using Case-4 version-3 Proposal vars

if (exists("q2.joint.post.res.1") && exists("q2.joint.post.res.2") && exists("q2.joint.post.res.3")){
iter.vec <- c(100,500,1000,5000,10000,50000)
iter.vec1 <- cumsum(iter.vec)
iter.vec1
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec1)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  #print((((i-1)==0)?(1):(iter.vec[i-1]+1)))
  print((ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)))
  # print((ifelse(test = (i-1)==0,yes = iter.vec[i],no = iter.vec[i-1]+1)))
  plot(x = 1:iter.vec[i],y = q2.joint.post.res.1$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q2.joint.post.res.2$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q2.joint.post.res.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()
}




### New Variances of Proposed distribution: Version-4
proposal.vars.5 <- c(0.004,0.005,0.03,0.002,0.015) # var(alpha), var(beta), var(gamma), var(theta), var(sigma).
proposal.vars.5

q2.res2.4 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.5,B = 5000,M = 50000)
q2.res2.4
accept.prob.func(df.lst = list(q2.res2.4)) # Ind.prob2 right
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.40932    0.36982    0.19250    0.32388 
diag(var(q2.res2.4[,4:8]))


# Q.2. Gelman-Rubin statistic of Different Iterations, with B=5000, and Case-4, Version-4 of Proposed Variances
#q2.conv.res.4 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars.5,B = 5000,M.iter = c(seq(from=5000,to=10000,by=1000),seq(from=20000,to=50000,by=10000)))
if (exists("q2.conv.res.4")){
q2.conv.res.4
q2.conv.res.4[which.min(q2.conv.res.4$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 9             40000               1.898834
}

# Q.2. Gelman-Rubin statistic of Different Iterations, with B=10000, and Case-4, Version-4 of Proposed Variances
#q2.conv.res.5 <- conv.diag.func(initpars = list(q1.pars.init.1,q1.pars.init.2,q1.pars.init.3),data = XYdata.3,prior.pars = prior.pars,proposalpars = proposal.vars.5,B = 10000,M.iter = c(seq(from=10000,to=50000,by=10000)))
if (exists("q2.conv.res.5")){
q2.conv.res.5
q2.conv.res.5[which.min(q2.conv.res.5$`Gelman-Rubin statistic`),]
# Num_of_Iterations Gelman-Rubin statistic
# 1             10000               1.937835
}


# Q.2. Combined Joint Posterior for Different Number of samples b= 5000,(Q1.Case-2 Initial pars), using Case-4 version-4 Proposal vars
#q2.joint.post.res.4 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.2,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.5,B = 5000)
if (exists("q2.joint.post.res.4")){
q2.joint.post.res.4
accept.prob.func(df.lst = list(q2.joint.post.res.4))
# $Df1  # For Total 66600 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.4193694  0.3483483  0.2095045  0.3257407 
accept.prob.func(df.lst = list(q2.joint.post.res.4[16601:66600,]))
# #Last 50,000 rows 
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.40072    0.34974    0.19292    0.31446 
}


# Q.2. Combined Joint Posterior for Different Number of samples b= 5000,(Q1.Case-3 Initial pars), using Case-4 version-4 Proposal vars
#q2.joint.post.res.5 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.3,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.5,B = 5000)
if (exists("q2.joint.post.res.5")){
q2.joint.post.res.5
accept.prob.func(df.lst = list(q2.joint.post.res.5))
# $Df1  # For Total 66600 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.4216517  0.3678979  0.2186486  0.3360661 
accept.prob.func(df.lst = list(q2.joint.post.res.5[16601:66600,]))
# Last 50,000 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.4239800  0.3590000  0.2052400  0.3294067 
}


"prop.var.func.v1" <- function(num=5,from=0,to=1,initpars, prior.pars, data,B=0,M=100,niter=1000,dur.secs=120){
  
  
 df1 <- data.frame(matrix(nrow=0,ncol=length(initpars)))
  df.res <- data.frame(matrix(nrow=0,ncol=length(initpars)))
  
  propvar <- c()
  cnt <- 0
  start <- 0
  end <- 0
  start <- Sys.time()
  repeat{
    while (TRUE){
    #propvar <- round(runif(n = num,min = from,max = to),3)
      propvar <- runif(n = num,min = from,max = to)
      if (nrow(df1)>0){
        if (!identical(x = propvar,y = df1[1:nrow(df1),])){
          df1 <- rbind(df1, propvar)
          break
        }
      } else {
        df1 <- rbind(df1, propvar)
        break
      }
    }
    cat("Proposal Vars=",propvar,fill=T,"\n")
    jo.post <- freundpost(pars = initpars,dat = data,priorpars = prior.pars,proposalpars = propvar,B = B,M = M)
    prop.i1 <- sum(jo.post[,1]==1)/M
    prop.i2 <- sum(jo.post[,2]==1)/M
    prop.i3 <- sum(jo.post[,3]==1)/M
    cat("prop1=",prop.i1,"prop2",prop.i2,"prop3",prop.i3,fill=T,"\n")
    cat("prop1=",sum(jo.post[,1]==1)/(length(jo.post[,1])),"prop2",sum(jo.post[,2]==1)/(length(jo.post[,2])),"prop3",sum(jo.post[,3]==1)/(length(jo.post[,3])),fill=T,"\n")
    if (prop.i1>=0.20 && prop.i1<=0.40 && prop.i2>=0.20 && prop.i2<=0.40 && prop.i3>=0.20 && prop.i3<=0.40){
      
      df.res <- rbind(df.res,propvar)
      cat("Posterior Variances",propvar,fill=T,"\n")
    } 
    else if (prop.i1>=0.20 && prop.i1<=0.40 && prop.i2>=0.20 && prop.i2<=0.40){
      df.res <- rbind(df.res,c(propvar[1:3],rep(NA,2)))
    }
    else if (prop.i3>=0.20 && prop.i3<=0.40 && prop.i2>=0.20 && prop.i2<=0.40){
      df.res <- rbind(df.res,c(rep(NA,2),propvar[3:5]))
    }
    else if (prop.i3>=0.20 && prop.i3<=0.40 && prop.i1>=0.20 && prop.i1<=0.40){
    df.res <- rbind(df.res,c(propvar[1:2],NA,propvar[4:5]))
    }
    else if (prop.i1>=0.20 && prop.i1<=0.40){
      df.res <- rbind(df.res,c(propvar[1:2],rep(NA,3)))
    }
    else if (prop.i2>=0.20 && prop.i2<=0.40){
      df.res <- rbind(df.res,c(rep(NA,2),propvar[3],rep(NA,2)))
    }
    else if (prop.i3>=0.20 && prop.i3<=0.40){
      df.res <- rbind(df.res,c(rep(NA,3),propvar[4:5]))
    }
    
    cnt <- cnt+1
    
    if (cnt==niter){
      cat("Number of default iterations ",niter,"exceeded",fill=T,"\n")
      break}
    end <- Sys.time()
    cat("Total Time elpased till Iteration:",cnt,"=",difftime(time1 = end,time2 = start,units = "secs")[[1]],"seconds",fill=T,"\n")
    # if ((difftime(time1 = end,time2 = start,units = "secs")[[1]])>=dur.secs){  
    #   cat("Default Time Limit (in seconds)",dur.secs,"exceeded",fill=T,"\n")
    #   break
    # }
  }
  
  if (ncol(df1)>0){
    colnames(df1) <- c(paste("Var.Par",1:length(initpars),sep = ""))
  }
  if (ncol(df.res)>0){
    colnames(df.res) <- c(paste("Var.Par",1:length(initpars),sep = ""))
  }
  
  
  return (list(proposal.var.init=df1,Proposal.Variances.final=df.res,size.init.propvar=nrow(df1),size.final.prop.var=nrow(df.res)))
}

### New Variances of Proposed distribution: Version-5

#prop.var.func.v1(initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,B = 5000,M = 50000) # No Result

#prop.var.func.v1(to = 0.5,initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,B = 5000,M = 50000) # No Result

#prop.var.func.v1(to = 0.1,initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,B = 5000,M = 50000) # No Result

#prop.var.func.v1(from = 0.001,to = 0.1,initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,B = 5000,M = 50000,niter = 10) # No Result

#prop.var.func.v1(to = 0.1,initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,B = 5000,M = 50000,niter = 10) 

## Output from one of the runs of above function

# $Proposal.Variances.final
# Var.Par1   Var.Par2    Var.Par3  Var.Par4    Var.Par5
# 1          NA         NA 0.021713165        NA          NA
# 2 0.021438119 0.01542436 0.035516988        NA          NA
# 3 0.088259557 0.01231290          NA        NA          NA
# 4 0.008802654 0.03825724 0.052718347        NA          NA
# 5          NA         NA          NA 0.0226799 0.008153796
# 6          NA         NA 0.019400881        NA          NA
# 7          NA         NA 0.009531331        NA          NA


#q2.res2.5 <- prop.var.func.v1(to = 0.1,initpars = q1.pars.init.1,prior.pars = prior.pars,data = XYdata.3,B = 5000,M = 50000,niter = 10)
#q2.res2.5
# $Proposal.Variances.final
# Var.Par1    Var.Par2    Var.Par3 Var.Par4 Var.Par5
# 1         NA          NA 0.054021535       NA       NA
# 2 0.06496289 0.009217158          NA       NA       NA
# 3         NA          NA 0.006377275       NA       NA
# 4         NA          NA 0.031039415       NA       NA
# 5         NA          NA 0.001309950       NA       NA
# 6 0.04836955 0.074471408 0.049856769       NA       NA
# 7         NA          NA 0.078873819       NA       NA
# 8 0.09319848 0.025585282 0.058015121       NA       NA


### New Variances of Proposed distribution: Version-5
proposal.vars.6 <- c(0.008802654, 0.03825724, 0.052718347, 0.0226799, 0.008153796) # var(alpha), var(beta), var(gamma), var(theta), var(sigma). # These values are obtained from enumeration function above (Combination of Results of the Previous two attempts shown above).
proposal.vars.6

q2.joint.post.res.6 <- joint.post.func(niter = c(100,500,1000,5000,10000,50000),initpars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.6,B = 5000)
accept.prob.func(df.lst = list(q2.joint.post.res.6))
# $Df1  # For Total 66600 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.4216517  0.3678979  0.2186486  0.3360661 
accept.prob.func(df.lst = list(q2.joint.post.res.6[16601:66600,]))
# Last 50,000 Rows
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.2069800  0.3925800  0.1641000  0.2545533 


###################################################

##Final Question-2 Summary of Case-4 analysis onwards

# Goal: To tweak the proposed variance values in order to get the acceptance probabilities from the joint posterior model in and around 0.2 to 0.4

# Results: It will comprise of acceptance probabilities (primarily), and sometimes with trace plots, and Gelman- statistics, for M= 50,000 samples, and B=5000 or 10,000 for different combinations of proposed variances.

#Proposed variances V-4.2 (0.04,0.5,0.03,0.02,0.015)
# Case-1 Initial Values B=5000, M=50,000
#ind.prob1  ind.prob2  ind.prob3 total.prob 
#0.18136    0.34518    0.19124    0.23926 

# Proposed Variances V-4.4 (0.004,0.005,0.03,0.002,0.015)
### Case-1 Initial Values B=5000 M=50,000
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.40932    0.36982    0.19250    0.32388
### Case-2 Initial values B=5000 M=50,000 (Last Rows)
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.40072    0.34974    0.19292    0.31446 
### Case-3 Initial values B=5000 M=50,000 (Last Rows)
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.4239800  0.3590000  0.2052400  0.3294067 

### B=5000 (All three initial values) M=50,000
# Num_of_Iterations Gelman-Rubin statistic
# 9             40000               1.898834
### B=10000 (All three initial values) M=50,000
# Num_of_Iterations Gelman-Rubin statistic
# 1             10000               1.937835



#Proposed Variances V-4.3 (0.04,0.05,0.03,0.002,0.015)
### Case-1 Initial Values B=5000, M=50,000
#ind.prob1  ind.prob2  ind.prob3 total.prob 
#0.20730    0.34640    0.23728    0.26366 
### Case-1 Initial Values B=5000, M=50,000 (Last Rows)
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.2075600  0.3517400  0.2293800  0.2628933 
### Case-2 Initial Values B=5000, M=50,000 (Last Rows)
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.1977600  0.3546400  0.2323800  0.2615933 
### Case-3 Initial Values B=5000, M=50,000 (Last Rows)
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.19462    0.34924    0.23728    0.26038   


### B=5000 (All three initial values) M=50,000
# Num_of_Iterations Gelman-Rubin statistic
# 15              1000               1.530241
### B=5000 (All three initial values) M=50,000 (from 5000)
# Num_of_Iterations Gelman-Rubin statistic
# 5              9000                1.95191
### B=10000 (All three initial values) M=50,000 (from 10000)
# Num_of_Iterations Gelman-Rubin statistic
# 2             20000                2.07693


## Q.2. Trace plots for M=c(100,500,1000,5000,10000,50000) , B=5000  using Case-4 version-3 Proposal vars

if (exists("q2.joint.post.res.1") && exists("q2.joint.post.res.2") && exists("q2.joint.post.res.3")){
iter.vec <- c(100,500,1000,5000,10000,50000)
iter.vec1 <- cumsum(iter.vec)
iter.vec1
windows(width = 10,height=10)
par(mfrow=c(3,2))
for (i in 1:length(iter.vec1)){
  lbl.txt <- paste("Trace_plot_of_alpha for",iter.vec[i],"iterations",sep=" ")
  #print((((i-1)==0)?(1):(iter.vec[i-1]+1)))
  print((ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)))
  # print((ifelse(test = (i-1)==0,yes = iter.vec[i],no = iter.vec[i-1]+1)))
  plot(x = 1:iter.vec[i],y = q2.joint.post.res.1$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],type='l',col="red",xlab="No. of Iterations",ylab="Posterior_alpha",main=lbl.txt)
  lines(x = 1:iter.vec[i],y = q2.joint.post.res.2$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="blue")
  lines(x = 1:iter.vec[i],y = q2.joint.post.res.3$alp[(ifelse(test = (i-1)==0,yes = 1,no = iter.vec1[i-1]+1)):iter.vec1[i]],col="green")
}
par(mfrow=c(1,1))
dev.off()
}

### Trace Plot of Individual Parameters of Case-4 Version-3 of Proposed Variances, using Case-1 Initial values.
windows(width = 10,height = 10)
par(mfrow=c(3,2))
plot(x = 1:50000,y = q2.res2.3$alp,type='l',col="red",main = "For alpha")
# lines(x = 1:50000,y = q1.res2.1$bet,col="blue")
# lines(x = 1:50000,y = q1.res3.1$bet,col="green")
plot(x = 1:50000,y = q2.res2.3$bet,type='l',col="red",main = "For Beta")
plot(x = 1:50000,y = q2.res2.3$gam,type='l',col="red",main = "For Gamma")
plot(x = 1:50000,y = q2.res2.3$thet,type='l',col="red",main = "For Theta")
plot(x = 1:50000,y = q2.res2.3$sig,type='l',col="red",main = "For Sigma")
par(mfrow=c(1,1))
dev.off()


# To summarize the results in Q.2., we have decided to choose the Proposed Variance V-4.3 (0.04,0.05,0.03,0.002,0.015), with a Burn-in Period (B) = 5000, for M=50,000 samples. The Burn-in period has been chosen based on trace plots, and also the Gelman-Rubin statistic value was lower than that in the case of B=10000, for the same proposed variances. Also, the acceptance probabilities using this version of proposed variances were within the range of 0.20 to 0.40, for all the three Markov chains.

## So, the subsequent results in Q.3. and Q.4. in the subsequent analysis section, have been done using the Joint Posterior Model for Q1 Case-1 Initial values, Q2 V-4.3 Proposed Variances, B=5000. M=50,000.


#####################################################






## Question -3

# 5-number summary # M=50,000

# Markov Chain-1
summary(q1.res1.1$alp)
summary(q1.res1.1$bet)
summary(q1.res1.1$gam)
summary(q1.res1.1$thet)
summary(q1.res1.1$sig)
hist(q1.res1.1$alp)
hist(q1.res1.1$bet)
hist(q1.res1.1$gam)
hist(q1.res1.1$thet)
hist(q1.res1.1$sig)
q3.res1 <- rbind(Alpha=summary(q1.res1.1$alp),Beta=summary(q1.res1.1$bet),Gamma=summary(q1.res1.1$gam),Theta=summary(q1.res1.1$thet),Sigma=summary(q1.res1.1$sig))
q3.res1
density(q1.res1.1$alp)$y
str(density(q1.res1.1$alp))
plot(density(q1.res1.1$alp))
#plot(x = q1.res1.1$alp,y = density(q1.res1.1$alp)$y)
plot(density(q1.res1.1$alp)$y)
hist(x = q1.res1.1$alp,main = "distribution of alpha")

hist(x = q1.res1.1$alp,probability = T,main = "distribution of alpha")
lines(density(q1.res1.1$alp),col="blue")
lines(density(q1.res1.1$alp,adjust = 2),col="red")
legend("topright",legend = c("Original Density","Adjusted Density"),col = c("blue","red"),lty = c(1,1))


# Markov Chain-2
summary(q1.res2.1$alp)
summary(q1.res2.1$bet)
summary(q1.res2.1$gam)
summary(q1.res2.1$thet)
summary(q1.res2.1$sig)

# Markov Chain-3
summary(q1.res3.1$alp)
summary(q1.res3.1$bet)
summary(q1.res3.1$gam)
summary(q1.res3.1$thet)
summary(q1.res3.1$sig)


## 95% Credible Intervals # M=50,000

# Chain-1
quantile(x = q1.res1.1$alp,probs = c(0.025,0.975))
quantile(x = q1.res1.1$bet,probs = c(0.025,0.975))
quantile(x = q1.res1.1$gam,probs = c(0.025,0.975))
quantile(x = q1.res1.1$thet,probs = c(0.025,0.975))
quantile(x = q1.res1.1$sig,probs = c(0.025,0.975))
hist(q1.res1.1$alp,main="95% credible Intervals for alpha")
abline(v = quantile(x = q1.res1.1$alp,probs = c(0.025,0.975)),col=c("red","blue"))
legend("topright",legend = c("Lower_95%_Cred_CI","Upper_95%_Cred_CI"),col = c("red","blue"),lty = c(1,1),cex=0.5)

# Chain-2
quantile(x = q1.res2.1$alp,probs = c(0.025,0.975))
quantile(x = q1.res2.1$bet,probs = c(0.025,0.975))
quantile(x = q1.res2.1$gam,probs = c(0.025,0.975))
quantile(x = q1.res2.1$thet,probs = c(0.025,0.975))
quantile(x = q1.res2.1$sig,probs = c(0.025,0.975))

# Chain-3
quantile(x = q1.res3.1$alp,probs = c(0.025,0.975))
quantile(x = q1.res3.1$bet,probs = c(0.025,0.975))
quantile(x = q1.res3.1$gam,probs = c(0.025,0.975))
quantile(x = q1.res3.1$thet,probs = c(0.025,0.975))
quantile(x = q1.res3.1$sig,probs = c(0.025,0.975))




## Correlations in the Joint Posterior, when M=50,000

# Chain-1

cor(q1.res1.1[,-c(1:3)]) # Correlation
cov(q1.res1.1[,-c(1:3)]) # Covariance
diag(cov(q1.res1.1[,-c(1:3)])) # Variance
solve(cov(q1.res1.1[,-c(1:3)])) # Inverse of Covariance Matrix
diag(solve(cov(q1.res1.1[,-c(1:3)])))

# Chain-2

cor(q1.res2.1[,-c(1:3)]) # Correlation
cov(q1.res2.1[,-c(1:3)]) # Covariance
diag(cov(q1.res2.1[,-c(1:3)])) # Variance

# Chain-3

cor(q1.res3.1[,-c(1:3)]) # Correlation
cov(q1.res3.1[,-c(1:3)]) # Covariance
diag(cov(q1.res3.1[,-c(1:3)])) # Variance

"freq.CI" <- function(pars,se){
  
  if (!is.data.frame(pars)){
    cat("Parameter Input Should be a Dataframe from the Joint Posterior","\n")
    stop()
  } 
  else if (length(se)!=ncol(pars)){
    cat("Parameter Dataframe Columns Should be Equal to length of Standard Error vector","\n")
    stop()
  }
  median.flag <- FALSE
  count<- 0
  mean.res <- c()
  med.res <- c()
  
  repeat{
  mean.est <- c()
  mean.lowerCI <- c()
  mean.upperCI <- c()
  mean.widthCI <- c()
  for ( i in 1:ncol(pars)){
    mean.est[i] <- ifelse(test = median.flag==FALSE,yes = mean(pars[,i]),no = median(pars[,i]))
    mean.lowerCI[i] <- mean.est[i] - 1.96*se[i]
    mean.upperCI[i] <- mean.est[i] + 1.96*se[i]
    mean.widthCI[i] <- mean.upperCI[i] - mean.lowerCI[i]
  }
  if (median.flag==F){
    mean.res <- cbind(LowerCI=mean.lowerCI,Estimate=mean.est,UpperCI=mean.upperCI,WidthCI=mean.widthCI)
   # mean.res <- cbind(c(names(pars)),mean.res)
   # print(row.names(mean.res))
  #  print(names(pars))
    row.names(mean.res) <- c(names(pars))
    median.flag=TRUE
    count <- count+1
   # print(mean.res)
  }
  else {
    med.res <- cbind(LowerCI=mean.lowerCI,Estimate=mean.est,UpperCI=mean.upperCI,WidthCI=mean.widthCI)
   # med.res <- cbind(c(names(pars)),med.res)
    row.names(med.res) <- c(names(pars))
    count <- count+1
   # print(med.res)
  }
  if (count==2){break}
  }
  
  return (list(MeanCI=round(mean.res,4),MedianCI=round(med.res,4)))
}

freq.CI(pars = q1.res1.1[,-c(1:3)],se = c(0.6809, 0.3023, 0.0529, 0.1597, 0.1313))


"post.infer.func" <- function(X,Y,final.pars){
  r <- 0
  c <- 0
  
  if (is.data.frame(final.pars)){
    r <- dim(final.pars)[1]
    c <- dim(final.pars)[2]
  } else {
    cat("Parameter Input Should be a Dataframe from the Joint Posterior","\n")
    stop()
  }
  
  corr.col <- c()
  
  
  for (i in 1:c){
    if (any(final.pars[,i]!=0 & final.pars[,i]!=1)){
      corr.col <- c(corr.col,i)    
    }
  }
  final.pars <- final.pars[,corr.col]
  
  cat("First Few Rows of Dataframe of Parameters",fill=T,"\n")
  print(head(final.pars))
  
  #-------------------------------------------------
  
  
  summary.res1 <- rbind(Alpha=summary(final.pars$alp),Beta=summary(final.pars$bet),Gamma=summary(final.pars$gam),Theta=summary(final.pars$thet),Sigma=summary(final.pars$sig))
  
  credCI <- rbind(Alpha=quantile(x = final.pars$alp,probs = c(0.025,0.975)),Beta=quantile(x = final.pars$bet,probs = c(0.025,0.975)),Gamma=quantile(x = final.pars$gam,probs = c(0.025,0.975)),Theta=quantile(x = final.pars$thet,probs = c(0.025,0.975)),Sigma=quantile(x = final.pars$sig,probs = c(0.025,0.975)))
  credCI <- cbind(credCI,width.int=c(credCI[1,2]-credCI[1,1],credCI[2,2]-credCI[2,1],credCI[3,2]-credCI[3,1],credCI[4,2]-credCI[4,1],credCI[5,2]-credCI[5,1]))
  
  corr <- cor(final.pars)
  
  covar <- cov(final.pars)
  
  var.pars <- diag(covar)
  se.pars <- sqrt(var.pars)
  
  var.se.pars <- rbind(Variance=var.pars,Std.Error=se.pars)
  
  MeanCI.res <- freq.CI(pars = final.pars,se = se.pars)$MeanCI
  MedianCI.res <- freq.CI(pars = final.pars,se = se.pars)$MedianCI
  
  #---------------Plots-------------------
  windows(width=10,height=10)
  par(mfrow=c(5,2))
  for (i in 1:ncol(final.pars)){
  #par(mfrow=c(1,2))
    label.txt1 <- paste("distribution of",names(final.pars)[i],sep = " ")
  hist(x = final.pars[,i],probability = T,main = label.txt1)
  lines(density(final.pars[,i]),col="blue")
  lines(density(final.pars[,i],adjust = 2),col="red")
  legend("topright",legend = c("Original Density","Adjusted Density"),col = c("blue","red"),lty = c(1,1),cex=0.5)
  
  label.txt2 <- paste("95% credible Intervals for",names(final.pars)[i],sep = " ")
  hist(final.pars[,i],main=label.txt2)
 # abline(v = quantile(x = q1.res1.1$alp,probs = c(0.025,0.975)),col=c("red","blue"))
  abline(v=as.numeric(c(credCI[i,1],credCI[i,2])),col=c("red","blue"))
  legend("topright",legend = c("Lower_95%_Cred_CI","Upper_95%_Cred_CI"),col = c("red","blue"),lty = c(1,1),cex=0.5)
  
  #par(mfrow=c(1,1))
  
  }
  par(mfrow=c(1,1))
  #----------------------------------------
  
  
  return (list(num.summary = round(summary.res1,4),cred.CI=round(credCI,4),joint.post.corr=round(corr,4),joint.post.cov=round(covar,4),var_se_pars=round(var.se.pars,4),freq.MeanCI=MeanCI.res,freq.MedianCI=MedianCI.res))
}


post.infer.func(X = XYdata.3$x,Y = XYdata.3$y,final.pars = q1.res1.1)

post.infer.func(X = XYdata.3$x,Y = XYdata.3$y,final.pars = q1.res2.1)

post.infer.func(X = XYdata.3$x,Y = XYdata.3$y,final.pars = q1.res3.1)


### Export Inference Results into Latex

q3infer.res1 <- post.infer.func(X = XYdata.3$x,Y = XYdata.3$y,final.pars = q1.res1.1)
q3infer.res1

print(xtable::xtable(x = q3infer.res1$num.summary, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_v1.tex")

print(xtable::xtable(x = q3infer.res1$cred.CI, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_v2.tex")

print(xtable::xtable(x = q3infer.res1$joint.post.corr, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_v3.tex")

print(xtable::xtable(x = q3infer.res1$var_se_pars, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_v4.tex")


#xtable::print.xtableList(xtable::xtableList(x = q3infer.res1, type = "latex", tabular.environment="longtable",colnames.format="multiple"), file = "Q3_Infer_v2.tex")



#### Subsequent Analysis of Question-3 based on Final Summaries of question-1 and Question-2

# #Proposed Variances V-4.3 (0.04,0.05,0.03,0.002,0.015) proposal.vars.4

# Q1 Case-1 Initial Value Parameters (0.2222222 2.0000000 0.2222222 0.1000000 0.7071068)  q1.pars.init.1 

##B=5000,M=50,000
prior.pars # Prior Parameters
proposal.vars.4 # Proposed Variances
q1.pars.init.1 # Initial Value of Parameters

q3.res2 <- freundpost(pars = q1.pars.init.1,dat = XYdata.3,priorpars = prior.pars,proposalpars = proposal.vars.4,B = 5000,M = 50000)
q3.res2 # # Here, sigma has original values.

accept.prob.func(df.lst = list(q3.res2))
### Attempt-1
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.1935600  0.3543000  0.2218000  0.2565533 
### Attempt-2
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.19044    0.36366    0.20076    0.25162 
### Attempt-3
# ind.prob1  ind.prob2  ind.prob3 total.prob 
# 0.19848    0.34700    0.22864    0.25804 

q3.res2.1 <- q3.res2
q3.res2.1
str(q3.res2.1)
q3.res2.1$sig <- q3.res2.1$sig^2
q3.res2.1 # Here, sigma has squared Values

q3infer.res2 <- post.infer.func(X = XYdata.3$x,Y = XYdata.3$y,final.pars = q3.res2.1)
q3infer.res2  # Here, sigma has squared Values
q3infer.res2$cred.CI
q3infer.res2$freq.MeanCI


### HPD Intervals 


### Version using R-Packages
#library(gmodels)
#library(Rmisc)
#sapply(q3.res2.1[, 4:8], Rmisc::CI)
#Rmisc::
#gmodels::
#?sapply

#library(BayesTwin)


##### To DO: HPD Intervals

## Refer the following:

# (Bayesian HPD Interval) (Important)
# google - calculate HPD intervals in R
# 
# (Bayesian HPD Interval) (Important)
# https://stats.stackexchange.com/questions/381520/how-can-i-estimate-the-highest-posterior-density-interval-from-a-set-of-x-y-valu
# 
# (Bayesian HPD Interval) (Important)
# https://www.reddit.com/r/AskStatistics/comments/xtbjxp/how_are_hpd_intervals_actually_computed/
#   
#   (Bayesian HPD Interval) (Important)
# https://search.r-project.org/CRAN/refmans/BayesTwin/html/HPD.html

## Refer to STAT 520 Lecture notes, as well.


### Export Inference Results into Latex --New Version Q3

print(xtable::xtable(x = q3infer.res2$num.summary, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_new_v1.tex")

print(xtable::xtable(x = q3infer.res2$cred.CI, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_new_v2.tex")

print(xtable::xtable(x = q3infer.res2$joint.post.corr, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_new_v3.tex")

print(xtable::xtable(x = q3infer.res2$var_se_pars, type = "latex", tabular.environment="longtable",digits=4), file = "Q3_Infer_new_v4.tex")




## Question-4


### Scatter-plot with Median Response function and Mean Response function

# With Respect to Markov-Chain-1,M=50,000 (q1.res1.1)

summary(q1.res1.1$alp)
median(q1.res1.1$alp)
mean(q1.res1.1$alp)

# Median Response function

post.par.est.median <- c(median(q1.res1.1$alp),median(q1.res1.1$bet),median(q1.res1.1$gam),median(q1.res1.1$thet),median(q1.res1.1$sig))
post.par.est.median

freundlich(xs = XYdata.3$x,pars = post.par.est.median)
mu.i.func(X = XYdata.3$x,par = post.par.est.median)
identical(x = freundlich(xs = XYdata.3$x,pars = post.par.est.median),y = mu.i.func(X = XYdata.3$x,par = post.par.est.median))

med.resp.est <- mu.i.func(X = XYdata.3$x,par = post.par.est.median)
med.resp.est

# Confidence Bands for Median Response function

q4.med.der <- der.func(X = XYdata.3$x,par = post.par.est.median[1:3])
q4.med.cov <- cov(q1.res1.1[,-c(1:3)])[1:3,1:3]
q4.med.vcov <- q4.med.der%*%q4.med.cov%*%t(q4.med.der)
q4.med.var <- diag(q4.med.vcov)
q4.med.se <- sqrt(q4.med.var)
q4.med.LowerCI <- med.resp.est-1.96*q4.med.se
q4.med.upperCI <- med.resp.est+1.96*q4.med.se
cbind(q4.med.LowerCI,med.resp.est,q4.med.upperCI)

plot(x = XYdata.3$x,y = XYdata.3$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(q4.med.LowerCI),max(q4.med.upperCI)),main="Median Response function")
lines(x = XYdata.3$x[order(XYdata.3$x)],y = med.resp.est[order(XYdata.3$x)],col="red",lwd=2)
lines(x = XYdata.3$x[order(XYdata.3$x)],y =  q4.med.LowerCI[order(XYdata.3$x)],col="blue",lty=1,lwd=1.5)
lines(x = XYdata.3$x[order(XYdata.3$x)],y = q4.med.upperCI[order(XYdata.3$x)],col="blue",lty=1,lwd=1.5)



# Mean Response function

post.par.est.mean <- c(mean(q1.res1.1$alp),mean(q1.res1.1$bet),mean(q1.res1.1$gam),mean(q1.res1.1$thet),mean(q1.res1.1$sig))
post.par.est.mean

mean.resp.est <- mu.i.func(X = XYdata.3$x,par = post.par.est.mean)
mean.resp.est


# Confidence Bands for Mean Response function

q4.mu.der <- der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])
q4.mu.cov <- cov(q1.res1.1[,-c(1:3)])[1:3,1:3]
q4.mu.vcov <- q4.mu.der%*%q4.mu.cov%*%t(q4.mu.der)
q4.mu.var <- diag(q4.mu.vcov)
q4.mu.se <- sqrt(q4.mu.var)
q4.mu.LowerCI <- mean.resp.est-1.96*q4.mu.se
q4.mu.upperCI <- mean.resp.est+1.96*q4.mu.se
cbind(q4.mu.LowerCI,mean.resp.est,q4.mu.upperCI)

plot(x = XYdata.3$x,y = XYdata.3$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(q4.mu.LowerCI),max(q4.mu.upperCI)),main="Mean Response function")
lines(x = XYdata.3$x[order(XYdata.3$x)],y = mean.resp.est[order(XYdata.3$x)],col="green",lwd=2)
lines(x = XYdata.3$x[order(XYdata.3$x)],y =  q4.mu.LowerCI[order(XYdata.3$x)],col="blue",lty=1,lwd=1.5)
lines(x = XYdata.3$x[order(XYdata.3$x)],y = q4.mu.upperCI[order(XYdata.3$x)],col="blue",lty=1,lwd=1.5)




# %*%t(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3]))
# diag(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])%*%cov(q1.res1.1[,-c(1:3)])[1:3,1:3]%*%t(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])))
# sqrt(diag(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])%*%cov(q1.res1.1[,-c(1:3)])[1:3,1:3]%*%t(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3]))))
# cbind(mean.resp.est-1.96*sqrt(diag(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])%*%cov(q1.res1.1[,-c(1:3)])[1:3,1:3]%*%t(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])))),mean.resp.est,mean.resp.est+1.96*sqrt(diag(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])%*%cov(q1.res1.1[,-c(1:3)])[1:3,1:3]%*%t(der.func(X = XYdata.3$x,par = post.par.est.mean[1:3])))))



plot(x = XYdata.3$x,y = XYdata.3$y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(XYdata.3$y),max(mean.resp.est)),main="Median versus Mean Response function")
lines(x = XYdata.3$x[order(XYdata.3$x)],y = med.resp.est[order(XYdata.3$x)],col="red",lwd=2)
lines(x = XYdata.3$x[order(XYdata.3$x)],y = mean.resp.est[order(XYdata.3$x)],col="green",lwd=2)
legend("bottomright",legend = c("Median Response","Mean Response"),col = c("red","green"),lwd = c(2,2))


"resp.func" <- function(X,Y,final.pars){
  
  r <- 0
  c <- 0
  
  if (is.data.frame(final.pars)){
    r <- dim(final.pars)[1]
    c <- dim(final.pars)[2]
  } else {
    cat("Parameter Input Should be a Dataframe from the Joint Posterior","\n")
    stop()
  }
  
  corr.col <- c()
  
  
  for (i in 1:c){
    if (any(final.pars[,i]!=0 & final.pars[,i]!=1)){
      corr.col <- c(corr.col,i)    
    }
  }
  final.pars <- final.pars[,corr.col]
  
  cat("First Few Rows of Dataframe of Parameters",fill=T,"\n")
  print(head(final.pars))
  
  med.par <- c()
  for (i in 1:ncol(final.pars)){
    med.par[i] <- median(final.pars[,i])
  }
  mean.par <- c()
  for (i in 1:ncol(final.pars)){
    mean.par[i] <- mean(final.pars[,i])
  }
  
  mean.resp <- mu.i.func(X = X,par = mean.par)
  
  mu.der <- der.func(X = X,par = mean.par[1:3])
  print(mu.der)
  #mu.cov <- cov(q1.res1.1[,-c(1:3)])[1:3,1:3]
  mu.cov <- cov(final.pars)[1:3,1:3]
  mu.vcov <- mu.der%*%mu.cov%*%t(mu.der)
  mu.var <- diag(mu.vcov)
  mu.se <- sqrt(mu.var)
  mu.LowerCI <- mean.resp-1.96*mu.se
  mu.upperCI <- mean.resp+1.96*mu.se
  #cbind(q4.mu.LowerCI,mean.resp.est,q4.mu.upperCI)
  
  med.resp <- mu.i.func(X=X,par = med.par)
  
  med.der <- der.func(X = X,par = med.par[1:3])
  med.cov <- cov(final.pars)[1:3,1:3]
  med.vcov <- med.der%*%med.cov%*%t(med.der)
  med.var <- diag(med.vcov)
  med.se <- sqrt(med.var)
  med.LowerCI <- mean.resp-1.96*med.se
  med.upperCI <- mean.resp+1.96*med.se
  #cbind(q4.med.LowerCI,mean.resp.est,q4.med.upperCI)
  
  
  
  
  
  windows(width = 10,height = 10)
  par(mfrow=c(3,1))
  plot(x = X,y = Y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(mu.LowerCI),max(mu.upperCI)),main="Mean Response function")
  lines(x = X[order(X)],y = mean.resp[order(X)],col="green",lwd=2)
  lines(x = X[order(X)],y =  mu.LowerCI[order(X)],col="blue",lty=1,lwd=1.5)
  lines(x = X[order(X)],y = mu.upperCI[order(X)],col="blue",lty=1,lwd=1.5)
  
  plot(x = X,y = Y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(med.LowerCI),max(med.upperCI)),main="Median Response function")
  lines(x = X[order(X)],y = med.resp[order(X)],col="red",lwd=2)
  lines(x = X[order(X)],y =  med.LowerCI[order(X)],col="blue",lty=1,lwd=1.5)
  lines(x = X[order(X)],y = med.upperCI[order(X)],col="blue",lty=1,lwd=1.5)
  
  plot(x = X,y = Y,xlab = "Covariate(X):Phosphorus content of Soil",ylab="Response(Y):Soil Sorption",ylim = c(min(Y,min(mean.resp,med.resp)),max(Y,max(mean.resp,med.resp))),main="Median versus Mean Response function")
  lines(x = X[order(X)],y = med.resp[order(X)],col="red",lwd=2)
  lines(x = X[order(X)],y = mean.resp[order(X)],col="green",lwd=2)
  legend("bottomright",legend = c("Median Response","Mean Response"),col = c("red","green"),lwd = c(2,2))
  
  par(mfrow=c(1,1))
  
  
  return(list(median.pars=med.par,median.resp=cbind(med.LowerCI,med.resp,med.upperCI),mean.pars=mean.par,mean.resp=cbind(mu.LowerCI,mean.resp,mu.upperCI)))
}

resp.func(X=XYdata.3$x,Y = XYdata.3$y,final.pars = q1.res1.1)
resp.func(X=XYdata.3$x,Y = XYdata.3$y,final.pars = q1.res2.1)
resp.func(X=XYdata.3$x,Y = XYdata.3$y,final.pars = q1.res3.1)


#### Subsequent Analysis of Question-4 based on Final Summaries of question-1 and Question-2 and Q.3.



q3.res2.1 # Here,sigma has squared Values Obtained From Q.3.

q4.res1 <- resp.func(X = XYdata.3$x,Y = XYdata.3$y,final.pars = q3.res2.1)
q4.res1  # Here, sigma has squared Values
q4.res1$median.pars
q4.res1$median.resp
q4.res1$mean.resp

### Export Inference Results into Latex --New Version Q4

print(xtable::xtable(x = head(q4.res1$median.resp), type = "latex", tabular.environment="longtable",digits=4), file = "Q4_Med_Resp_new_v1.tex")

print(xtable::xtable(x = head(q4.res1$mean.resp), type = "latex", tabular.environment="longtable",digits=4), file = "Q4_Mean_Resp_new_v1.tex")








## Question-5

### Comparison Results have been mentioned in Overleaf Latex Report of Project-3







graphics.off()


save.image("601_Project3.RData")

# },error=function(e) {message(e)}
# , warning=function(w) {message(w)
#   invokeRestart("muffleWarning")})
# },error = function(e){
#   message(e)
# },warning = function(w){
#   message(w)
# },finally =  {
#   cat("\n","Done with execution of the RScript File.","\n")
# }
# )

