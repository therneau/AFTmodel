// Using STAN for a PIB fit
// B[1]- B[5] are the parameters of the target function, two columns
// y ~ Gaussian(f(age + offset), sigma)
// Add covariates
//  The only difference from 18d5 is that we have 4 covariates
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
  int  outcome[N]; // is each Y a pointer to pib or tau
  real age[N];     // age values
  matrix [N,P] x;  // covariates
  int  id[N];     // the id index, 1,2 etc
  real tau;
  real dft;
}
parameters {
  real<lower=0> B[5,2];     // parameters of the logistic
  real<lower=0> sigma[2];   // residual std
  matrix[M,2] alpha;          // per subject intercepts
  matrix[P,2] beta;
  real rho;
}
transformed parameters {
   matrix[N,2] lin;   // linear predictors 
   vector[N] tage;
   vector[N] yhat;
   matrix[2,2] Sigma;
   row_vector[2] zero;
   Sigma[1,1] = tau*tau;
   Sigma[2,2] = tau*tau;
   Sigma[1,2] = rho*(tau*tau);
   Sigma[2,1] = rho*(tau*tau);

   zero = rep_row_vector(0,2);
   lin = x*beta;
   for (i in 1:N) {
     tage[i] = (age[i] + alpha[id[i],outcome[i]] + lin[i, outcome[i]])/10 - B[4,outcome[i]];
     yhat[i] = ghyper(B[1,outcome[i]], B[2, outcome[i]], B[3,outcome[i]], 
         B[5, outcome[i]], tage[i]);
    }
}
model {
  for (k in 1:2){
    B[,k] ~ normal(0, 50);   //vague prior, sd of 50, all our values are < 10
    beta[,k] ~ normal(0,50); //vague again, all values are < 20
    }
  sigma ~ normal(0, 1); // the value is near .04, so this is no constraint
  rho  ~  beta(2,1);
  for (j in 1:M) {
     alpha[j,] ~ multi_student_t(dft, zero, Sigma);
     }
   for (i in 1:N) {
      y[i]  ~ normal(yhat[i], sigma[outcome[i]]);
    }
}
