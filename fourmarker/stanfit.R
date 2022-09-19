#
# Fit the Stan model, using the rstan library
#
library(rstan)
simdata <- readRDS("simdata.rds") # load the simulated data

# Create a list object from the data, which is what the call needs.
# id needs to be 1, 2, ... with no holes
# outcome needs to be 1, 2, 3, or 4
# center education at 16
uid <- sort(unique(simdata$id))
if (any(simdata$id != match(simdata$id, uid))) stop("id is not sequential")

stan.list <- list(id = simdata$id, y= simdata$y, 
               outcome= as.numeric(simdata$outcome),
               age = simdata$age,
               x = cbind(simdata$male, simdata$apoe, simdata$apoemis,
                         simdata$educ-16),
               adrc = simdata$adrc,
               N = nrow(simdata),  # number of obs
               M = length(uid),    # number of unique subjects
               P = 3,              # number of covariates
               dft = 10            # df for the t distribution
)

# Set up initial values.  This should not be strictly necessary, but we 
#  get a shorter burn in with rational choices for the shape
set.seed(39392)  # for reproducable initial values
init.fun <- function(chain) {
    target <- cbind(c(.4, .05, .4, 7.5, .1), # amyloid
                    c(.3, .05, .4, 8.5, .2), # tau
                    c(.1, .01, .4, 8.5, .5), # wmh
                    c(.4, .02, .4, 8.5, .15)) # 1-FA
    list(B= target * runif(length(target), .8, 1.2),
         sigma= abs(rnorm(4, 0, .2)))
}

nchains <- 6
init.list <- lapply(1:nchains, init.fun)

fit4<- stan("aft4.stan",
            data= stan.list,
            init= init.list,
            chains= nchains,
            iter = 10000,
            cores = nchains,
            pars = c("alpha", "B", "beta", "sigma", "Omega", "atau", 
                     "adrc_shift")
            )

# save(fit4, "fit4save.rda")

# Example output
sumfit <- summary(fit4)[[1]]
rname <- rownames(sumfit)

B <- matrix(sumfit[grepl("B", rname), "50%"], byrow=TRUE, ncol=4,
            dimnames= list(paste0("B", 1:5), c("amyloid", "tau", "WMH", "FA")))
round(B, 4)

beta <- matrix(sumfit[grepl("beta", rname), "50%"], byrow=TRUE, ncol=4,
            dimnames= list(c("male", "APOE+", "APOE-missing", "educ"), 
                           c("amyloid", "tau", "WMH", "FA")))
round(beta, 4)

 
                 
