// Using STAN for a four parameter fit
// B[1]- B[5] are the parameters of the target function, two columns
// y ~ Gaussian(f(age + offset), sigma)
functions{
  real ghyper(real b1, real b2, real b3, real k, real x) {
      return(b1 + b2*x + .5*(b3-b2)*(x + sqrt(x^2 + k^2)));
  }
}
data { 
  int<lower=0> N;  // number of observations
  int<lower=0> M;  // number of persons (subjects)
  int<lower=0> P;  // number of variables
  real y[N];       // observed values
  int  outcome[N]; // is each Y a pointer to amyloid, tau, WMH or FA
  real age[N];     // age values
  matrix [N,P] x;  // covariates
  int  id[N];      // the id index, 1,2 etc
  real adrc[N]; 
  real dft;        // degrees of freedom of the bivariate t
}

parameters {
  real<lower=0> B[5,4];     // parameters of the logistic
  real<lower=0> sigma[4];   // residual std
  matrix[M,4] alpha;        // per subject intercepts
  matrix[P,4] beta;
  corr_matrix[4] Omega;     // correlation between random effects
  vector<lower=0>[4] atau;  // scale for random effect
  real adrc_shift;          // referral effect assumed to be the same for all 4
}
transformed parameters {
   matrix[N,4] lin;   // linear predictors 
   vector[N] tage;
   vector[N] yhat;
   row_vector[4] zero;

   zero = rep_row_vector(0,4);
   lin = x*beta;  
   for (i in 1:N) { 
       tage[i] = (age[i] + alpha[id[i],outcome[i]] + lin[i, outcome[i]] + 
             adrc_shift * adrc[i]);
       if (outcome[i] ==4) yhat[i] = B[1,4] + B[2,4]*tage[i];
       else {
          tage[i] = tage[i]/10 - B[4, outcome[i]]; 
          yhat[i] = ghyper(B[1,outcome[i]], B[2, outcome[i]], B[3,outcome[i]], 
                B[5, outcome[i]], tage[i]);
       }
   }    
}
model {
  atau ~ cauchy(0, 2.5); 
  Omega ~ lkj_corr(2);   // correlation matrix
  for (k in 1:4){
    B[,k] ~ normal(0, 15);   //vague prior
    beta[,k] ~ normal(0,30); //vague again, all values are < 20
    }
  {
    matrix[4, 4] Sigma_beta;    // not worth inspecting, review Omega and atau
    Sigma_beta = quad_form_diag(Omega, atau);
    for (j in 1:M) {
       alpha[j,] ~ multi_student_t(dft, zero, Sigma_beta);
       }
  }
  for (i in 1:N) {
      y[i]  ~ normal(yhat[i], sigma[outcome[i]]);
    }
  adrc_shift ~ normal(0, 50); // <---- Model the common ADRC shift
}

