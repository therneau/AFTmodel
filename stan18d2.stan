// Using STAN for a PIB fit
// B[1]- B[5] are the parameters of the target function f
// y ~ Gaussian(f(age + offset), sigma)
// Add covariates
functions{
  real ghyper(real b1, real b2, real b3, real k, real x) {
      return(b1 + b2*(x+2) + .5*(b3-b2)*(x + sqrt(x^2 + k^2)));
  }
}
data {
  int<lower=0> N;  // number of observations
  int<lower=0> M;  //number of persons (subjects)
  int<lower=0> P;  // number of variables
  real y[N];       // observed values
  real age[N];     // age values
  matrix [N,P] x;  // covariates
  int  id[N];     // the id index, 1,2 etc
}
parameters {
  real<lower=0> B[5];     // parameters of the logistic
  real<lower=0> tau;      // varariance of the subject intercepts
  real<lower=0> sigma;    // variance of the data about the means
  real alpha[M];          // per subject intercepts
  real beta[P];
}
transformed parameters{
   vector[N] tage;
   vector[N] yhat;
   for (i in 1:N) {
     tage[i] = (age[i] + alpha[id[i]]  + beta[1]*x[i,1] + beta[2]*x[i,2]
                       + beta[3] * x[i,3])/10 - B[4];
     yhat[i] = ghyper(B[1], B[2], B[3], B[5], tage[i]);
    }
}
model {
  B ~ normal(0, 50);   //vague prior, sd of 50, all our values are < 10
  beta ~ normal(0, 20);
  tau ~ normal(0, 5); // stronger constraint on tau
  alpha ~ student_t(3, 0, tau);
  y  ~ normal(yhat, sigma);
}
