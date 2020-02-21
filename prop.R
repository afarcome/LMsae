prop <- function(xv,k){

# compute the proportion of each of the k categories
  n = length(xv)
  pv = rep(0,k)
  for(h in 1:k) pv[h] = sum(xv==h)/n
  return(pv)

}