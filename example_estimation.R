# load data
load("sample_simulated.RData")

# number of iteration of the MCMC estimation algorithm
nit = 10000

# estimate with different numbers of classes
est1 = estimation(Y,X,Z,M,cN,k=1,nit=nit,taube=0.001,tauga=0.01,taude=0.01)
est2 = estimation(Y,X,Z,M,cN,k=2,nit=nit,taube=0.001*2,tauga=0.01,taude=0.01)
est3 = estimation(Y,X,Z,M,cN,k=3,nit=nit,taube=0.001*3,tauga=0.01,taude=0.01)
est4 = estimation(Y,X,Z,M,cN,k=4,nit=nit,taube=0.001*4,tauga=0.01,taude=0.01)
est5 = estimation(Y,X,Z,M,cN,k=5,nit=nit,taube=0.001*5,tauga=0.01,taude=0.01)

# save results
save.image("example_estimation.RData")