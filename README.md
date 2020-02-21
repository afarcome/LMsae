R and Fortran scripts and functions implementing the approach described in the paper
"A hidden Markov space-time model for mapping world access to food dynamics"

The main script files are:

- example_estimation.R: to run the Bayesian estimation algorithm on example data that have been simulated;

- example_simulation.R: to simulate a dataset;

- compile.R: compile all R and Fortran functions (to be run before the examples).

Suitable comments are directly included in the three scripts in order to illustrate their use in R and intput/output arguments.

Note that the commands to compile and link the Fortran functions are written for Mac users. Simple adjustments are required to compile and link them in other platforms, following the manual available at https://cran.r-project.org/doc/manuals/r-release/R-exts.html




