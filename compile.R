# compile all R functions
require("mvtnorm")
source("prop.R")
source("rmultbern.R")
source("simulate.R")
source("estimation.R")

# compile and link Fortran codes (written for Mac users)
system("R CMD SHLIB simulate_u.f --preclean")
system("R CMD SHLIB lk_be3.f --preclean")
system("R CMD SHLIB update_u4.f --preclean")
system("R CMD SHLIB lk_ga2.f --preclean")
system("R CMD SHLIB lk_de2.f --preclean")
dyn.load("simulate_u.so")
dyn.load("lk_be3.so")
dyn.load("update_u4.so")
dyn.load("lk_ga2.so")
dyn.load("lk_de2.so")

