## R code for fitting Stan model described in "Relationships between
## \beta-amyloid and tau in an elderly population: An accelerated
## failure time model" by Terry Therneau et al.
## Note: For local Mayo run we need `scl enable devtoolset-8 bash`
##       prior to `R CMD BATCH aft-stan-fit.R`.
library(rstan)

# This loads aft.list.stan (data is not provided), a list containing the
#  data elements N= number of PET scans, M= number of subjects,
#  id = vector of subject id values (1,1, 2,2,2, 3, ...), x= matrix with
#  columns for apoe, sex, and education, adrc = vector of 0=MCSA 1=ADRC,
#  y = log(amyloid or tau), outcome= 0/1 indicates amyloid or tau, 
#  dft = degrees of freedom for the t distribution
load("stan.rda")


mod <- stan("aft-model.stan",
            data = aft.list.stan,
	    pars = c("alpha", "B", "beta", "sigma", "rho", "adrc_shift"),
            chains = 4,
            iter = 5000 * 2,
            init = "random",
            thin = 1,
            seed = 39392,
	    include = TRUE,
            cores = 4)

save(mod, file = "aft-model-20210415.rda")
