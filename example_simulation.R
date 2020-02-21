# number of sites and time occasions
n = 1000
TT = 8

# matrix of the time-varying sample size of each site
la = 100
M = matrix(rpois(n*TT,la),n,TT)

# binary matrix of neighbors
C = matrix(0,n,n)
for(j in 1:n) C[j,] = runif(n)<1/n
C = pmin(C+t(C),1)
diag(C) = 0

# sets of neighboors of each unit
nnv = rowSums(C)
cN = matrix(0,n,max(nnv))
for(j in 1:n) if(nnv[j]>0) cN[j,1:nnv[j]] = which(C[j,]==1)

# number of latent states
k = 3

# covariates at unit and site levels
X = array(NA,c(2,max(M),n,TT))
for(j in 1:n) for(t in 1:TT) X[,1:M[j,t],j,t] = rnorm(2*M[j,t])/sqrt(2)
Z = array(0,c(3,n,TT))
for(j in 1:n) for(t in 1:TT) Z[,j,t] = rnorm(3)/sqrt(2)

# regression parameters (covariates and latent states) for the responses
np1 = 1+dim(X)[1]
Be = matrix(c(0.5,0,-0.5,-1,0,1,-1,0,1),k,np1)
dimnames(Be) = list(state=1:k,covariate=1:np1)

# regression parameters (covariates and neighbors) for initial time occasion
np2 = dim(Z)[1]+k
Ga = matrix(c(-1,-1,-1,-1,1,0,0,1,2,1,1,2),k-1,np2)
dimnames(Ga) = list(state=2:k,covariate=1:np2)
De = array(c(-3,-3,-1,-1,1,0,0,1,2,1,1,2,
             -1,-2,1,0,-1,-1,0,1,-2,-1,-1,1,
             -1,-2,1,0,0,1,-1,-1,-1,1,-2,-1),c(k-1,np2,k))
dimnames(De) = list(state=2:k,covariate=1:np2,prev_state=1:k)

# simulate data
out = simulate(Be,Ga,De,X,Z,M,cN,nnv,nit=10000)
Y = out$Y; U = out$U
save(file="sample_simulated.RData",Y,X,Z,M,cN)