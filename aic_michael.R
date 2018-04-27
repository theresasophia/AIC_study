## Taken from v69i12.R
library(pomp)

##' ## More complex models.
##' ### Simple SIR.
##' C snippets expressing the two faces of the measurement model.

rmeas <- "
  cases = rnbinom_mu(theta, H);
"
dmeas <- "
  lik = dnbinom_mu(cases, theta, H, give_log);
"

##' The process model simulator.
##' This takes one step from time t -> t+dt
##' The per-capita rates of the elementary transitions are stored in 'rate'.
##' The numbers of individuals making each transition is stored in 'trans'.
##' Births are Poisson, transitions are Euler-multinomial.
##' 'H' accumulates the recoveries (and will be zeroed after each observation).

sir.step <- "
  double rate[6];
  double dN[6];
  double P;
  P = S + I + R;

  rate[0] = mu * P;       // birth
  rate[1] = beta * I / P; // transmission
  rate[2] = mu;           // death from S
  rate[3] = gamma;        // recovery
  rate[4] = mu;           // death from I
  rate[5] = mu;           // death from R

  dN[0] = rpois(rate[0] * dt);
  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
  reulermultinom(1, R, &rate[5], dt, &dN[5]);

  //State development
  S += dN[0] - dN[1] - dN[2];
  I += dN[1] - dN[3] - dN[4];
  R += dN[3] - dN[5];
  H += dN[1];
"

##' Construct the pomp object and fill with simulated data.
fromEstScale <- Csnippet("
 Tbeta = exp(beta);
 Tgamma = exp(gamma);
 Tmu = exp(mu);
 Ttheta = exp(theta);
")

toEstScale <- Csnippet("
 Tbeta = log(beta);
 Tgamma = log(gamma);
 Tmu = log(mu);
 Ttheta = log(theta);
")

param_vec1 <- c(popsize = 5e+05, beta = 400, gamma = 26, mu = 1/50,
               theta = 100, S.0 = 26/400, I.0 = 0.002, R.0 = 1)

sir1 <- pomp(data = data.frame(cases = NA,
             time = seq(0, 10, by = 1/52)),
             times = "time",
             t0 = -1/52,
             dmeasure = Csnippet(dmeas),
             rmeasure = Csnippet(rmeas),
             rprocess = euler.sim(step.fun = Csnippet(sir.step),
                                  delta.t = 1/52/20),
             statenames = c("S", "I", "R", "H"),
             paramnames = c("gamma", "mu", "theta", "beta", "popsize", "S.0", "I.0", "R.0"),
             zeronames = c("H"),
             initializer = function(params, t0, ...) {
               fracs <- params[c("S.0", "I.0", "R.0")]
               setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R",
                                                                             "H"))
             },
             toEstimationScale = toEstScale,
             fromEstimationScale = fromEstScale,
             params = param_vec1)

##Simulate a trajectory
sir1 <- simulate(sir1, seed = 1914679908L)

ops <- options(scipen = -10)
plot(sir1, mar = c(0, 5, 2, 0))
options(ops)

######################################################################
##Fit model!
######################################################################

##Fitting using iterated filtering
m1 <-  mif2(
  sir1,
  start=param_vec1, ##start at true values
  Np=2000,
  Nmif=50,
  cooling.type="geometric",
  cooling.fraction.50=0.5,
  transform=TRUE,
  rw.sd=rw.sd(beta=0.02,gamma=0.02,mu=0.02, theta=0.02))

summary(m1)

plot(m1)

cbind(fit=coef(m1), true=param_vec1)

##Loglikelihood of the filter approximation (not the same as the model!)
logLik(m1)

##Model likelihood at MLE found by MIF
fit1 <- pfilter(sir1,params=coef(m1),Np=5000)
ll1 <- logLik(fit1)


######################################################################
## Model with waiting time being gamma instead of exponential
######################################################################

sir.step2 <- "
  double rate[10];
  double dN[10];
  double P;

  //Total
  P = S + I + (R1 + R2 + R3);

  rate[0] = mu * P;       // birth
  rate[1] = beta * I / P; // transmission
  rate[2] = mu;           // death from S

  //From I->R (now R1, R2, R3) to mimic gamma
  rate[3] = gamma/3 ;        // recovery I->R1
  rate[4] = mu;           // death from I

  rate[5] = gamma/3 ;        // recovery R1->R2
  rate[6] = mu;           // death from R1

  rate[7] = gamma/3 ;        // recovery R2->R3
  rate[8] = mu;           // death from R2

  rate[9] = mu;           // death from R3

  dN[0] = rpois(rate[0] * dt);

  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
//New R states
  reulermultinom(2, R1, &rate[5], dt, &dN[5]);
  reulermultinom(2, R2, &rate[7], dt, &dN[7]);
  reulermultinom(1, R3, &rate[9], dt, &dN[9]);

  //State development: N[0] = birth, N[1] = new infected, N[2] = Dead from S, N[3] = I->R1, N[4] = Dead from I, N[5] = R1->R2, N[6]: dead from R1, ...
  S += dN[0] - dN[1] - dN[2];
  I += dN[1] - dN[3] - dN[4];
  R1 += dN[3] -  dN[5] - dN[6];
  R2 += dN[5] -  dN[7] - dN[8];
  R3 += dN[7] - dN[9];
  H += dN[1];
"

param_vec2 <- c(popsize = 5e+05, beta = 400, gamma = 26, mu = 1/50,
               theta = 100, S.0 = 26/400, I.0 = 0.002, R1.0 = 1, R2.0 = 0, R3.0 = 0)

sir2 <- pomp(data = data.frame(cases = NA,
             time = seq(0, 10, by = 1/52)),
             times = "time",
             t0 = -1/52,
             dmeasure = Csnippet(dmeas),
             rmeasure = Csnippet(rmeas),
             rprocess = euler.sim(step.fun = Csnippet(sir.step2),
                                  delta.t = 1/52/20),
             statenames = c("S", "I", "R1", "R2", "R3", "H"),
             paramnames = c("gamma", "mu", "theta", "beta", "popsize", "S.0", "I.0", "R1.0", "R2.0", "R3.0"),
             zeronames = c("H"),
             initializer = function(params, t0, ...) {
               fracs <- params[c("S.0", "I.0", "R1.0", "R2.0", "R3.0")]
               setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R1", "R2", "R3", "H"))
             },
             toEstimationScale = toEstScale,
             fromEstimationScale = fromEstScale,
             params = param_vec2)

##Simulate a trajectory
sir2 <- simulate(sir2, seed = 1914679908L)
plot(sir2)

##Fitting using iterated filtering
m2 <-  mif2(
  sir2,
  start=param_vec2, ##start at true values
  Np=2000,
  Nmif=50,
  cooling.type="geometric",
  cooling.fraction.50=0.5,
  transform=TRUE,
  rw.sd=rw.sd(beta=0.02,gamma=0.02,mu=0.02, theta=0.02))

summary(m2)

plot(m2)

cbind(fit=coef(m2), true=param_vec2)

##Loglikelihood of the filter approximation (not the same as the model!)
logLik(m2)

##Model likelihood at MLE found by MIF
fit2 <- pfilter(sir2,params=coef(m2),Np=5000)
ll2 <- logLik(fit2)
ll2

######################################################################
##Fit model 1 to the data of model 2
######################################################################

sir1 <- pomp(data = data.frame(cases = as.numeric(sir2@data),
             time = seq(0, 10, by = 1/52)),
             times = "time",
             t0 = -1/52,
             dmeasure = Csnippet(dmeas),
             rmeasure = Csnippet(rmeas),
             rprocess = euler.sim(step.fun = Csnippet(sir.step),
                                  delta.t = 1/52/20),
             statenames = c("S", "I", "R", "H"),
             paramnames = c("gamma", "mu", "theta", "beta", "popsize", "S.0", "I.0", "R.0"),
             zeronames = c("H"),
             initializer = function(params, t0, ...) {
               fracs <- params[c("S.0", "I.0", "R.0")]
               setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R",
                                                                             "H"))
             },
             toEstimationScale = toEstScale,
             fromEstimationScale = fromEstScale,
             params = param_vec1)

m1 <-  mif2(
  sir1,
  start=param_vec1, ##start at true values
  Np=2000,
  Nmif=50,
  cooling.type="geometric",
  cooling.fraction.50=0.5,
  transform=TRUE,
  rw.sd=rw.sd(beta=0.02,gamma=0.02,mu=0.02, theta=0.02))

summary(m1)

plot(m1)

cbind(fit=coef(m1), true=param_vec1)

##Loglikelihood of the filter approximation (not the same as the model!)
logLik(m1)

##Model likelihood at MLE found by MIF
fit1 <- pfilter(sir1,params=coef(m1),Np=5000)
ll1 <- logLik(fit1)

##Compare the two likelihoods. Did the MIF2 run long enough in each case?
c(ll1, ll2)

##Which likelihood is bigger? Expectation is that 2 is better, but
##this could also be a convergence issue...
which.max(c(ll1, ll2))
