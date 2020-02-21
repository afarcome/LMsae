estimation <- function(Y,X,Z,M,cN,k,nit=1000,taube = 0.01,tauga = 0.25,taude = 0.25,
                       regany = 50,Be=NULL,Ga=NULL,De=NULL,U=NULL){
	
# Augmented MCMC algorithm for estimation of the model proposed in
# "A hidden Markov space-time model for mapping world access to food dynamics"
#
# INPUT:
# Y        = matrix of responses (n x TT)
# X        = covariate at individual level (np1 x max(M) x n x TT)
# Z        = covariate at site level (np2 x n x TT)
# M        = matrix of number of units in each site (n x TT)
# k        = number of latent states
# cN       = matrix of neighborhood for each unit (n x cmax)
# taube    = sd of the proposal for beta
# tauga    = sd of the proposal for gamma
# taude    = sd of the proposal for delta
# regany   = register and iteration any
# Be       = initial value of beta parameters (k x 1+np1)
# Ga       = initial value of gamma parameters (k-1 x 1+np2+k)
# De       = initial value of beta parameters (k-1 x 1+np2+k x k)
# U        = initial value of beta parameters (n x TT)
#
# OUTPUT:
# BE       = parameters affecting the conditional response probabilities (k x 1+np1 x nit)
# UU       = matrix of latent state (n x TT x nit)
# GA       = parameters affecting the initial probabilities (k-1 x 1+np2+k x nit)
# DE       = parameters affecting the transition probabilities (k-1 x 1+np2+k x k x nit)
# lklprv   = vector of loglikelihood+logprio (nit)
# PR       = posterior expected value for any observation (mmax x n x TT)
# eldpwv   = vector of eldpw values
# seldpwv  = vector of eldpw s.e.
# ELDPW    = single components of eldpw (mmax x n x TT)

t1 = proc.time()

# hyperparameters
	sibe = siga = side = sqrt(1000)

# preliminaries
  n = as.integer(dim(Y)[2]); TT = as.integer(dim(Y)[3]); k = as.integer(k)
  MM = array(0,c(k,k-1,k))
  for(u in 1:k) MM[,,u] = diag(k)[,-u]
  X1 = aperm(X,c(2,1,3,4)); X1[is.na(X1)] = 0
  M = matrix(as.integer(M),n,TT)
  cmax = ncol(cN)
  cN = matrix(as.integer(cN),n,cmax)
  nnv = as.integer(rowSums(cN>0))
  C = matrix(as.integer(0),n,n)
  for(i in 1:n) if(nnv[i]>0) C[i,cN[1:nnv[i]]] = 1
  ncov1 = dim(X)[1]; ncov2 = dim(Z)[1]
  np1 = as.integer(1+ncov1); np2 = as.integer(ncov2+k)
  mmax = as.integer(max(M))
  nnv = as.integer(nnv)
  nT = as.integer(sum(M))
  X2 = matrix(0,nT,ncov1); Y2 = rep(0,nT)
  co = 0
  for(t in 1:TT) for(j in 1:n) if(M[j,t]>0){
    X2[co+(1:M[j,t]),] = X1[1:M[j,t],,j,t]
    Y2[co+(1:M[j,t])] = Y[1:M[j,t],j,t]
    co = co+M[j,t]
  }
  rm(X); rm(X1); rm(Y)

# MCMC algorithm
  if(is.null(Be)){
    est = glm(Y2 ~ X2,family=binomial())
    Be = matrix(rmvnorm(k,est$coefficients,k*vcov(est)),k,np1)
  }
  if(is.null(Ga)) Ga = matrix(rnorm(np2*(k-1),0,tauga),k-1,np2)
  if(is.null(De)) De = array(rnorm(np2*(k-1)*k,0,taude),c(k-1,np2,k))
  if(is.null(U)) if(k==1) U = matrix(as.integer(1),n,TT) else U = matrix(sample(1:k,n*TT,rep=TRUE),n,TT)

# iterate
  accbe = accU = accga = accde = 0
  nit1 = nit/regany
  BE = array(0,c(k,np1,nit1))
  dimnames(BE) = list(state=1:k,covariate=1:np1,iteraton=1:nit1)
  UU = array(0,c(n,TT,nit1))
  dimnames(UU) = list(cluster=1:n,time=1:TT,iteraton=1:nit1)
  GA = array(0,c(k-1,np2,nit1))
  if(k>1) dimnames(GA) = list(state=2:k,covariate=1:np2,iteraton=1:nit1)
  DE = array(0,c(k-1,np2,k,nit1))
  if(k>1) dimnames(DE) = list(state=2:k,covariate=1:np2,prev_state=1:k,iteration=1:nit1)
  lklprv = rep(0,nit1)
  eldpwv = seldpwv = rep(0,nit1)
  PR = array(NA,c(mmax,n,TT))
  for(t in 1:TT) for(j in 1:n) PR[1:M[j,t],j,t] = 0
  dimnames(PR) = list(unit=1:mmax,cluster=1:n,time=1:TT)
  ELDPW = PR
  LPR = LPR2 = rep(0,nT)
  eldpwaic = seldpwaic = 0
  t0 = proc.time()[3]
  for(it in -(nit/4):nit){

# level of log-likelihodd + log prior
    lklpr = 0

# Update beta parameters
    Bes = Be+rnorm(k*np1,0,taube)
    out = .Fortran("lk_be3",Be,Bes,U,M,X2,Y2,n,TT,k,np1,nT,mmax,lkv=rep(0,k),lksv=rep(0,k))
    for(u in 1:k){
      beu = Be[u,]; beus = Bes[u,]
      lk = out$lkv[u]; lks = out$lksv[u]
      lpr = sum(dnorm(beu,0,sibe,log=TRUE))
      lprs = sum(dnorm(beus,0,sibe,log=TRUE))
      al = min(1,exp(lprs+lks-lpr-lk))
      if(runif(1)<al){
        Be[u,] = beus
        accbe = accbe+1/k
        lk = lks; lpr = lprs
      }
      lklpr = lklpr+lk+lpr
    }
    if(is.nan(lklpr)){
      print(1)
      browser()
    }

# update latent variables
    if(k>1){
      MGa = MM[,,1]%*%Ga
      MDe = array(0,c(k,np2,k))
      for(u in 1:k) MDe[,,u] = matrix(MM[,,u],k,k-1)%*%matrix(De[,,u],k-1,np2)
      out = .Fortran("update_u4",Be,k,np1,MGa,np2,MDe,cN,cmax,U=U,M,X2,mmax,nT,Z,Y2,n,TT,
                     nnv,accU=accU,tlrit=0,RR=array(0,c(k,n,TT)))
      U = out$U; accU = out$accU; tlrit = out$tlrit
      lklpr = lklpr+tlrit
    }
    if(is.nan(lklpr)){
      print(2)
      browser()
    }

# update gamma parameters
    if(k>1) for(u in 1:(k-1)){
      gau = Ga[u,]
      gaus = gau+rnorm(np2,0,tauga)
      Gas = Ga; Gas[u,] = gaus
      MGa = MM[,,1]%*%Ga
      MGas = MM[,,1]%*%Gas
      out = .Fortran("lk_ga2",MGa,k,np2,MGas,cN,cmax,n,U,TT,Z,nnv,lk=0,lks=0)
      lk = out$lk; lks = out$lks
      lpr = sum(dnorm(gau,0,siga,log=TRUE))
      lprs = sum(dnorm(gaus,0,siga,log=TRUE))
      al = min(1,exp(lprs+lks-lpr-lk))
      if(runif(1)<al){
        Ga = Gas
        accga = accga+1/(k-1)
        lk = lks; lpr = lprs
      }
      lklpr = lklpr+lk+lpr
    }
    if(is.nan(lklpr)){
      print(3)
      browser()
    }

# update delta parameters
    if(k>1) for(up in 1:k){
      for(u in 1:(k-1)){
        deupu = De[u,,up]
        deupus = deupu+rnorm(np2,0,tauga)
        Des = De; Des[u,,up] = deupus
        MDe = MDes = array(0,c(k,np2,k))
        for(u1 in 1:k){
          MDe[,,u1] = matrix(MM[,,u1],k,k-1)%*%matrix(De[,,u1],k-1,np2)
          MDes[,,u1] = matrix(MM[,,u1],k,k-1)%*%matrix(Des[,,u1],k-1,np2)
        }
        out = .Fortran("lk_de2",MDe,k,np2,MDes,C,n,U,TT,Z,nnv,up=as.integer(up),lk=0,lks=0)
        lk = out$lk; lks = out$lks
        lpr = sum(dnorm(deupu,0,side,log=TRUE))
        lprs = sum(dnorm(deupus,0,side,log=TRUE))
        al = min(1,exp(lprs+lks-lpr-lk))
        if(runif(1)<al){
          De = Des
          accde = accde+1/(k*(k-1))
          lk = lks; lpr = lprs
        }
        lklpr = lklpr+lk+lpr
        if(is.nan(lklpr)){
          print(4)
          browser()
        }
      }
    }

# record iteration
    if(it>0 & it%%regany==0){
      it1 = it/regany
      BE[,,it1] = Be
      UU[,,it1] = U
      GA[,,it1] = Ga
      DE[,,,it1] = De
      lklprv[it1] = lklpr
      co = 0
      for(t in 1:TT) for(j in 1:n) if(M[j,t]>0){
        ind = co+(1:M[j,t])
        if(M[j,t]==1){
          lpr = c(1,X2[ind,])%*%Be[U[j,t],]
        }else{
          lpr = cbind(1,X2[ind,])%*%Be[U[j,t],]					
        }
        pr = exp(lpr)/(1+exp(lpr))
        PR[1:M[j,t],j,t] = (PR[1:M[j,t],j,t]*(it1-1)+pr)/it1
        lpr1 = Y2[co+(1:M[j,t])]*log(pr)+(1-Y2[co+(1:M[j,t])])*log(1-pr)
        LPR[ind] = (LPR[ind]*(it1-1)+lpr1)/it1
        LPR2[ind] = (LPR2[ind]*(it1-1)+lpr1^2)/it1
        co = co+M[j,t]
      }
      co = 0
      ELDPW = array(NA,c(mmax,n,TT))
      for(t in 1:TT) for(j in 1:n) if(M[j,t]>0){
        tmp1 = Y2[co+(1:M[j,t])]*log(PR[1:M[j,t],j,t])+
          (1-Y2[co+(1:M[j,t])])*log(1-PR[1:M[j,t],j,t])
        tmp2 = LPR2[co+(1:M[j,t])]-LPR[co+(1:M[j,t])]^2
        ELDPW[1:M[j,t],j,t] = tmp1-tmp2
        co = co+M[j,t]
      }
      eldpwaic = sum(ELDPW,na.rm=TRUE)
      seldpwaic = sqrt(nT*var(as.vector(ELDPW),na.rm=TRUE))
      eldpwv[it1] = eldpwaic
      seldpwv[it1] = seldpwaic
    }

# display acceptance rate
    if(it%%10==0){
      tt = proc.time()[3]-t0
      it2 = it+nit/4
      print(c(iteration=it,k=k,accbe=accbe/it2,accU=accU/it2,accga=accga/it2,accde=accde/it2,
              lklpr=lklpr,eldpwaic=eldpwaic,seldpwaic=seldpwaic,tt/it2))
      print(Be)
    }
  }

# final output
  out = list(BE=BE,UU=UU,GA=GA,DE=DE,lklprv=lklprv,accbe=accbe,accU=accU,accga=accga,
             accde=accde,eldpwv=eldpwv,seldpwv=seldpwv,ELDPW=ELDPW)
  return(out)

}