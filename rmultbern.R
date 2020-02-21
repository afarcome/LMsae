rmultbern <- function(piv){

# draw a category from a multivariate Bernoulli distribution with
# parameter vector piv
  k = piv
  cpiv = cumsum(piv)
  x = 1+sum(runif(1)>cpiv)
  return(x)

}