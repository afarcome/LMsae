simulate <- function(Be,Ga,De,X,Z,M,cN,nnv,nit = 100){
	
# simulate data from the model proposed in the paper
# "A hidden Markov space-time model for mapping world access to food dynamics"
#
# INPUT:
# Be = parameters affecting the conditional response probabilities (k x np1)
# Ga = parameters affecting the initial probabilities (k-1 x np2)
# De = parameters affecting the transition probabilities (k-1 x 1+np2+k x k)
# X  = covariate at individual level (np1 x max(M) x n x TT)
# Z  = covariate at site level (ncov2 x n x TT)
# M  = matrix of number of units in each site (n x TT)
# cN      = matrix of neighborhood for each unit (n x cmax)
#
# OUTPUT:
# U  = matrix of latent state (n x TT)
# Y  = matrix of responses (n x TT)

# preliminaries
  k = as.integer(dim(Be)[1])
  n = as.integer(dim(X)[3])
  TT = as.integer(dim(X)[4])
  cmax = ncol(cN)
  cN = matrix(as.integer(cN),n,cmax)
  nnv = as.integer(rowSums(C))
  np2 = as.integer(dim(Ga)[2])
  MM = array(0,c(k,k-1,k))
  for(u in 1:k) MM[,,u] = diag(k)[,-u]

# iterate to draw matrix of latent states
  U = matrix(sample(1:k,n*TT,rep=TRUE),n,TT)
  MGa = MM[,,1]%*%Ga
  MDe = array(0,c(k,np2,k))
  for(u in 1:k) MDe[,,u] = MM[,,u]%*%De[,,u]
  for(it in 1:nit){
    if(it%%100==0) print(it)
    U = .Fortran("simulate_u",k,MGa,np2,MDe,cN,cmax,U=U,Z,n,TT,nnv,0)$U
  }

# draw binary observations
  Y = array(NA,c(max(M),n,TT))
  for(t in 1:TT) for(j in 1:n){
    lpv = cbind(1,t(X[,1:M[j,t],j,t]))%*%Be[U[j,t],]
    pv = exp(lpv)/(1+exp(lpv))
    Y[1:M[j,t],j,t] = 1*(runif(M[j,t])<pv) 
  }

# output
  out = list(U=U,Y=Y)

}