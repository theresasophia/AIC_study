rm(list = ls())     # clear objects
graphics.off()      #close graphics windows
#' Taken from v69i12.R
library(pomp)
require(doParallel)
library(magrittr)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(foreach)
cores <- 8
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

N <- 10  # number of simulations to generate
coef_list <- vector("list", N)

for(i in 1:N){

##' Simple SIRS.
##' C snippets expressing the two faces of the measurement model.

#'hoehle 2018-05-03: whenwhere is give_log defined?
dmeas <- "if (ISNA(cases)) {
                  lik = (give_log) ? 0 : 1;
                  } else {
                  lik =  dnbinom_mu(cases, theta, H, 1);
                  lik = (give_log) ? lik : exp(lik);}
"

rmeas <- "
 cases = rnbinom_mu(theta,H);
"

##' The process model simulator.
##' This takes one step from time t -> t+dt
##' The per-capita rates of the elementary transitions are stored in 'rate'.
##' The numbers of individuals making each transition is stored in 'trans'.
##' Births are Poisson, transitions are Euler-multinomial.
##' 'H' accumulates the recoveries (and will be zeroed after each observation).

sir.step <- "
  double rate[7];
  double dN[7];
  double P;
  P = S + I + R;

  rate[0] = mu * P;       // birth
  rate[1] = Beta * I / P; // transmission
  rate[2] = mu;           // death from S
  rate[3] = gamma;        // recovery
  rate[4] = mu;           // death from I
  rate[5] = mu;           // death from R
  rate[6] = omega;        // waning of immunity


  dN[0] = rpois(rate[0] * dt);
  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
  reulermultinom(2, R, &rate[5], dt, &dN[5]);


  //State development
  S += dN[0] - dN[1] - dN[2] + dN[6];
  I += dN[1] - dN[3] - dN[4];
  R += dN[3] - dN[5] - dN[6];
  H += dN[1];
"


fromEstScale <- Csnippet("
 TBeta = exp(Beta);
 Tgamma = exp(gamma);
 Tmu = exp(mu);
 Ttheta = exp(theta);
")

toEstScale <- Csnippet("
 TBeta = log(Beta);
 Tgamma = log(gamma);
 Tmu = log(mu);
 Ttheta = log(theta);
")

param_vec <- c(popsize = 10000, Beta = 4, gamma = 1, mu = 1/(80*52),
               theta = 100, S.0 = 90000, I.0 = 10, R.0 = 0, omega=1/4)


#'load the simulated data set - they were simulated form sir1 and sir 2 as explained below
read.table(paste("sir1_data",i,".txt", sep="")) %>%
  rbind(data.frame(time=0,cases=NA)) %>%
  arrange(time) -> sir1_dat
head(sir1_dat)

read.table(paste("sir2_data",i,".txt", sep="")) %>%
  rbind(data.frame(time=0,cases=NA)) %>%
  arrange(time) -> sir2_dat
head(sir2_dat)

read.table(paste("sir_seas_data",i,".txt", sep="")) %>%
  rbind(data.frame(time=0,cases=NA)) %>%
  arrange(time) -> sir_seas_dat
head(sir_seas_dat)

#' Construct the pomp object and fill with simulated data.
#' we start at t0=-10 so the system has time to equilibrate and reach the endemic level
sir1 <- pomp(data =sir1_dat,
             times = "time",
             t0 = -10,
             dmeasure = Csnippet(dmeas),
             rmeasure = Csnippet(rmeas),
             rprocess = euler.sim(step.fun = Csnippet(sir.step),
                                  delta.t = 1/40),
             statenames = c("S", "I", "R", "H"),
             paramnames = c("gamma", "mu", "theta", "Beta","omega", "popsize", "S.0", "I.0", "R.0"),
             zeronames = c("H"),
             initializer = function(params, t0, ...) {
               fracs <- params[c("S.0", "I.0", "R.0")]
               setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R",
                                                                             "H"))
             },
             toEstimationScale = toEstScale,
             fromEstimationScale = fromEstScale,
             params = param_vec)

 plot(sir1)
 #'This is how the data was created
 #'  we create 10 simulations, seed starts from 0, 01, etc last 01112345678= 8 and 01112345670 =9
 #' hoehle 2018-05-03: why not seed = 01112345670 + N?
 #'        a<- simulate(sir1,seed = 01112345670, obs=TRUE, as.data.frame=TRUE)
 #'        sir1_data<- data.frame(time=a$time[-1],cases=a$cases[-1])
 #'        plot(sir1_data$cases, type="l")
 #'        write.table(sir1_data, file = "~/Dropbox/AIC_study/AIC_study_git/sir1_data9.txt" )

######################################################################
##Fit model 1
######################################################################


sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- param_vec[c("gamma","omega","theta","mu","popsize","S.0","I.0","R.0")]


#' Fitting using iterated filtering with drawing 8 times from the sir_box
stew(file=paste("sir1_",i,".rda", sep=""),{

  t_global <- system.time({
    mifs_global <- foreach(i=1:cores,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
      sir1,
      start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
      Np=1000,
      Nmif=20,
      cooling.type="geometric",
      cooling.fraction.50=0.5,
      transform=TRUE,
      rw.sd=rw.sd(Beta=0.02))
    }
  })
},seed=1270401374,kind="L'Ecuyer")


##visualizing the diagnosistics of iterated filtering
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>%   mutate(variable = factor(variable)) %>% rename(f=L1) -> t
##hoehle 2018-04-28: easier and more consistent to use the dplyr rename
##function colnames(t)[4]<- "f".
##Not sure why it's necessary to rename though...

##Make the plot
ggplot(t,aes(x=iteration,y=value,color=variable,group=f))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()


## approximating the liklihood of every mif search outcome
stew(file=paste("sir1_lik_",i,".rda", sep=""),{

  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir1,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")



aic <- function(k,loglik){
  2*k-2*loglik
}



best <- which.max(liks_global[,1])
sir1_loglik <- round(liks_global[best,],3)
sir1_aic <- aic(1,sir1_loglik[1])


##Compare true parameters with fitted ones
cbind(fit=coef(mifs_global[[best]]), true=param_vec)


######################################################################
## Model with waiting time being gamma instead of exponential;
## the infectious period is now gamma(3,3*gamma) distributed with
## mean 1/gamma
######################################################################
sir.step2 <- "
double rate[11];
double dN[11];
double P;

//Total
P = S + (I1 + I2 + I3) + R;

rate[0] = mu * P;       // birth
rate[1] = Beta * (I1 + I2 + I3)/ P; // transmission
rate[2] = mu;           // death from S

//From I1->I2 (now I1, I2, I3) to mimic gamma
rate[3] = gamma*3 ;        // recovery I1->I2
rate[4] = mu;           // death from I1

rate[5] = gamma*3 ;        // recovery I2->I3
rate[6] = mu;           // death from I2

rate[7] = gamma*3 ;        // recovery I3->R
rate[8] = mu;           // death from I3

rate[9] = mu;           // death from R
rate[10] = omega;        // waning of immunity

dN[0] = rpois(rate[0] * dt);

reulermultinom(2, S, &rate[1], dt, &dN[1]);
reulermultinom(2, I1, &rate[3], dt, &dN[3]);
//New R states
reulermultinom(2, I2, &rate[5], dt, &dN[5]);
reulermultinom(2, I3, &rate[7], dt, &dN[7]);
reulermultinom(2, R, &rate[9], dt, &dN[9]);

//State development: N[0] = birth, N[1] = new infected, N[2] = Dead from S, N[3] = I->R1, N[4] = Dead from I, N[5] = R1->R2, N[6]: dead from R1, ...
S  += dN[0] - dN[1] - dN[2] + dN[10];
I1 += dN[1] - dN[3] - dN[4];
I2 += dN[3] - dN[5] - dN[6];
I3 += dN[5] - dN[7] - dN[8];
R  += dN[7] - dN[9] - dN[10];
H  += dN[1];
"

param_vec2 <- c(popsize = 10000, Beta = 4, gamma = 1, mu = 1/(80*52), omega = 1/4,
               theta = 100, S.0 = 90000, I1.0 = 10, I2.0 = 0, I3.0 = 0, R.0 = 0)

sir2 <- pomp(data = sir2_dat,
             times = "time",
             t0 = -10,
             dmeasure = Csnippet(dmeas),
             rmeasure = Csnippet(rmeas),
             rprocess = euler.sim(step.fun = Csnippet(sir.step2),
                                  delta.t = 1/40),
             statenames = c("S", "I1", "I2", "I3", "R", "H"),
             paramnames = c("omega","gamma", "mu", "theta", "Beta", "popsize", "S.0", "I1.0", "I2.0", "I3.0", "R.0"),
             zeronames = c("H"),
             initializer = function(params, t0, ...) {
               fracs <- params[c("S.0", "I1.0", "I2.0", "I3.0", "R.0")]
               setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S",  "I1", "I2", "I3", "R","H"))
             },
             toEstimationScale = toEstScale,
             fromEstimationScale = fromEstScale,
             params = param_vec2)

plot(sir2)
   #  a<- simulate(sir2,seed = 01234556789, obs=TRUE, as.data.frame=TRUE)
   #   sir2_data<- data.frame(time=a$time[-1],cases=a$cases[-1])
   #   plot(sir2_data$cases, type="l")
   # write.table(sir2_data, file = "~/Dropbox/AIC_study/AIC_study_git/sir2_data9.txt" )
 # #



sir_box <- rbind(
  Beta = c(0.5,10)
)

##better version IMHO -- less prone to copy paste errors:
sir_fixed_params <- param_vec2[c("gamma","theta","mu","popsize","S.0","I1.0","I2.0","I3.0","R.0","omega")]

##Fitting using iterated filtering
stew(file=paste("sir2_",i,".rda", sep=""),{

  t_global <- system.time({
    mifs_global <- foreach(i=1:cores,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir2,
        start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=1000,
        Nmif=20,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02))
    }
  })
},seed=1270401374,kind="L'Ecuyer")

#
#used dplyr to rename the column instead of separate call outside the pipeline
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>%   mutate(variable = factor(variable)) %>% rename(f="L1") -> t
##colnames(t)[4]<- "f"

ggplot(t,aes(x=iteration,y=value,color=variable,group=f))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()

stew(file=paste("sir2_lik_",i,".rda", sep=""),{

  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir2,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")



best <- which.max(liks_global[,1])
sir2_loglik <- round(liks_global[best,],3)
sir2_aic <- aic(1,sir2_loglik[1])


######################################################################
##Fit model 1 to the data of model 2
######################################################################

sir1_cross <- pomp(data =sir2_dat,
             times = "time",
             t0 = -10,
             dmeasure = Csnippet(dmeas),
             rmeasure = Csnippet(rmeas),
             rprocess = euler.sim(step.fun = Csnippet(sir.step),
                                  delta.t = 1/40),
             statenames = c("S", "I", "R", "H"),
             paramnames = c("gamma", "mu","omega", "theta", "Beta", "popsize", "S.0", "I.0", "R.0"),
             zeronames = c("H"),
             initializer = function(params, t0, ...) {
               fracs <- params[c("S.0", "I.0", "R.0")]
               setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R",
                                                                             "H"))
             },
             toEstimationScale = toEstScale,
             fromEstimationScale = fromEstScale,
             params = param_vec)


sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- param_vec[c("gamma","theta","mu","popsize","S.0","I.0","R.0","omega")]

##Fitting using iterated filtering
stew(file=paste("sir1_cross_",i,".rda", sep=""),{

  t_global <- system.time({
    mifs_global <- foreach(i=1:cores,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir1_cross,
        start=c(apply(sir_box,1,function(x) runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=1000,
        Nmif=20,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02))
    }
  })
},seed=1270401374,kind="L'Ecuyer")


mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>%
##  subset(iteration>17) %>%
  mutate(variable = factor(variable)) %>%
  rename(f="L1") -> t


ggplot(t,aes(x=iteration,y=value,color=variable,group=f)) +
  geom_line() +
  guides(color=FALSE) +
  labs(x="MIF2 Iteration",y="") +
  facet_wrap(~variable,scales="free_y",ncol=2) +
  theme_bw()


stew(file=paste("sir1_cross_lik_",i,".rda", sep=""),{

  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir1_cross,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")



best <- which.max(liks_global[,1])
sir1_cross_loglik <- round(liks_global[best,],3)
sir1_cross_aic <- aic(1,sir1_cross_loglik[1])

######################################################################
##Fit model 2 to the data of model 1
######################################################################

sir2_cross <- pomp(data = sir1_dat,
                   times = "time",
                   t0 = -10,
                   dmeasure = Csnippet(dmeas),
                   rmeasure = Csnippet(rmeas),
                   rprocess = euler.sim(step.fun = Csnippet(sir.step2),
                                        delta.t = 1/40),
                   statenames = c("S", "I1", "I2", "I3", "R", "H"),
                   paramnames = c("gamma", "mu", "theta", "Beta", "popsize", "S.0", "I1.0", "I2.0", "I3.0", "R.0", "omega"),
                   zeronames = c("H"),
                   initializer = function(params, t0, ...) {
                     fracs <- params[c("S.0", "I1.0", "I2.0", "I3.0", "R.0")]
                     setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I1", "I2", "I3", "R", "H"))
                   },
                   toEstimationScale = toEstScale,
                   fromEstimationScale = fromEstScale,
                   params = param_vec2)


sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- param_vec2[c("gamma","theta","mu","popsize","S.0","I1.0","I2.0","I3.0","R.0","omega")]


##Fitting using iterated filtering
stew(file=paste("sir2_cross_",i,".rda", sep=""),{

  t_global <- system.time({
    mifs_global <- foreach(i=1:cores,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir2_cross,
        start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=1000,
        Nmif=20,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02))
    }
  })
},seed=1270401374,kind="L'Ecuyer")

#
#
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>%   mutate(variable = factor(variable)) %>%
  rename(f="L1") -> t

ggplot(t,aes(x=iteration,y=value,color=variable,group=f))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()

stew(file=paste("sir2_cross_lik_",i,".rda", sep=""),{

  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir2_cross,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")



best <- which.max(liks_global[,1])
sir2_cross_loglik<- round(liks_global[best,],3)
sir2_cross_aic <- aic(1,sir2_cross_loglik[1])


#logliks of true model and other model
#is the difference only due to Monte Carlo noise?
sir2_loglik
sir1_cross_loglik

##Maybe better to check directly, if the intervals overlap
i2_on2 <- sir2_loglik[1] + c(-1,1) * qnorm(0.975) * sir2_loglik[2]
i1_on2 <- sir1_cross_loglik[1] + c(-1,1) * qnorm(0.975) * sir1_cross_loglik[2]

##Small helper function to assess if two intervals overlap. Assuming
##i1 and i2 are both sorted so i[1] < i[2].
overlap <- function(i1, i2) {
  (i1[2] >= i2[1]) & (i1[1] <= i2[2])
}
##They don't overlap, so AIC would not be the same? But we don't find
##the right model.. (AIC/loglik of cross fitted model is better?)
overlap(i2_on2, i1_on2)
sort(c(sir2_loglik=sir2_loglik[1], sir1_cross_loglik=sir1_cross_loglik[1]),decreasing=TRUE)


sir1_loglik
sir2_cross_loglik

sort(c(sir1_loglik=sir1_loglik[1], sir2_cross_loglik=sir2_cross_loglik[1]),decreasing=TRUE)



###########################################################################
##  SIRS with seasonaility                                                 
###########################################################################

sir.step_seas <- "
double rate[7];
double dN[7];
double P;
P = S + I + R;

rate[0] = mu * P;       // birth
rate[1] = Beta *(1 + rho * cos(M_2PI/5*t + phi))* I / P; // transmission
rate[2] = mu;           // death from S
rate[3] = gamma;        // recovery
rate[4] = mu;           // death from I
rate[5] = mu;           // death from R
rate[6] = omega;        // waning of immunity


dN[0] = rpois(rate[0] * dt);
reulermultinom(2, S, &rate[1], dt, &dN[1]);
reulermultinom(2, I, &rate[3], dt, &dN[3]);
reulermultinom(2, R, &rate[5], dt, &dN[5]);


//State development
S += dN[0] - dN[1] - dN[2] + dN[6];
I += dN[1] - dN[3] - dN[4];
R += dN[3] - dN[5] - dN[6];
H += dN[1];
"

param_vec_seas <- c(popsize = 10000, Beta = 4, gamma = 1, mu = 1/(80*52),
               theta = 100, S.0 = 90000, I.0 = 10, R.0 = 0, omega=1/4, phi=0, rho=0.6)


fromEstScale_seas <- Csnippet("
 TBeta = exp(Beta);
                         Tgamma = exp(gamma);
                         Tmu = exp(mu);
                         Ttheta = exp(theta);
                         Trho = expit(rho);
                         Tphi    = M_2PI*expit(phi);
                         ")

toEstScale_seas <- Csnippet("
                       TBeta = log(Beta);
                       Tgamma = log(gamma);
                       Tmu = log(mu);
                       Ttheta = log(theta);
                         Trho = logit(rho);
                         Tphi    = M_2PI*logit(phi);
                       ")



# we start at t0=-10 so the system has time to equilibrate and reach the endemic level
sir_seas <- pomp(data = sir_seas_dat,
             times = "time",
             t0 = -10,
             dmeasure = Csnippet(dmeas),
             rmeasure = Csnippet(rmeas),
             rprocess = euler.sim(step.fun = Csnippet(sir.step_seas),
                                  delta.t = 1/40),
             statenames = c("S", "I", "R", "H"),
             paramnames = c("gamma", "mu", "theta", "Beta","omega", "phi","rho","popsize", "S.0", "I.0", "R.0"),
             zeronames = c("H"),
             initializer = function(params, t0, ...) {
               fracs <- params[c("S.0", "I.0", "R.0")]
               setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R",
                                                                             "H"))
             },
             toEstimationScale = toEstScale_seas,
             fromEstimationScale = fromEstScale_seas,
             params = param_vec_seas)

# until 01122344566= 6, then 77 etc
# a<- simulate(sir_seas, seed=7789,obs=TRUE, as.data.frame=TRUE)
# sir1_data<- data.frame(time=a$time[-1],cases=a$cases[-1])
# plot(sir1_data$cases, type="l")
# write.table(sir1_data, file = "~/Dropbox/AIC_study/AIC_study_git/sir_seas_data9.txt" )


sir_box <- rbind(
  Beta = c(0.5,5),
  phi = c(0.001,0.1),
  rho = c(0.05,0.2)
)

sir_fixed_params <- param_vec_seas[c("gamma","omega","theta","mu","popsize","S.0","I.0","R.0")]


##Fitting using iterated filtering with drawing 8 times from the sir_box
stew(file=paste("sir_seas_",i,".rda", sep=""),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:8,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir_seas,
        start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=1000,
        Nmif=20,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02,rho=0.02, phi=0.02))
    }
  })
},seed=1270401374,kind="L'Ecuyer")


##visualizing the diagnosistics of iterated filtering
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","phi", "rho", "theta")) %>%
  melt() %>%  subset(iteration>0) %>%  mutate(variable = factor(variable)) %>% rename(f=L1) -> t

##Make the plot
ggplot(t,aes(x=iteration,y=value,color=variable,group=f)) +
  geom_line() +
  guides(color=FALSE) +
  labs(x="MIF2 Iteration",y="") +
  facet_wrap(~variable,scales="free_y",ncol=2) +
  theme_bw()


## approximating the liklihood of every mif search outcome
stew(file=paste("sir_seas_lik_",i,".rda", sep=""),{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir_seas,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")



best <- which.max(liks_global[,1])
sir_seas_loglik <- round(liks_global[best,],3)
sir_seas_aic <- aic(3,sir_seas_loglik[1])


##########################################
## Fit model 1 to data from model seas
##########################################

sir1_cross_seas <- pomp(data =sir_seas_dat,
                   times = "time",
                   t0 = -10,
                   dmeasure = Csnippet(dmeas),
                   rmeasure = Csnippet(rmeas),
                   rprocess = euler.sim(step.fun = Csnippet(sir.step),
                                        delta.t = 1/40),
                   statenames = c("S", "I", "R", "H"),
                   paramnames = c("gamma", "mu","omega", "theta", "Beta", "popsize", "S.0", "I.0", "R.0"),
                   zeronames = c("H"),
                   initializer = function(params, t0, ...) {
                     fracs <- params[c("S.0", "I.0", "R.0")]
                     setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R",
                                                                                   "H"))
                   },
                   toEstimationScale = toEstScale,
                   fromEstimationScale = fromEstScale,
                   params = param_vec)


sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- param_vec[c("gamma","theta","mu","popsize","S.0","I.0","R.0","omega")]

##Fitting using iterated filtering
stew(file=paste("sir1_cross_seas_",i,".rda", sep=""),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:8,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir1_cross_seas,
        start=c(apply(sir_box,1,function(x) runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=1000,
        Nmif=20,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02))
    }
  })
},seed=1270401374,kind="L'Ecuyer")


##[Q]hoehle 2018-04-28: why subset iteration>17?
##Also optimized the renaming
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>%
  subset(iteration>0) %>%
  mutate(variable = factor(variable)) %>%
  rename(f="L1") -> t
##colnames(t)[4]<- "f"

ggplot(t, aes(x=iteration,y=value,color=variable,group=f)) +
  geom_line() +
  guides(color=FALSE) +
  labs(x="MIF2 Iteration",y="") +
  facet_wrap(~variable,scales="free_y",ncol=2) +
  theme_bw()


stew(file=paste("sir1_cross_seas_lik_",i,".rda", sep=""),{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir1_cross_seas,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")


best <- which.max(liks_global[,1])
sir1_seas_cross_loglik <- round(liks_global[best,],3)
sir1_cross_seas_aic <- aic(1,sir1_seas_cross_loglik[1])




##################################################################
## Fit model seas to data from Model 1
##################################################################


sir_seas_cross <- pomp(data = sir1_dat,
                 times = "time",
                 t0 = -10,
                 dmeasure = Csnippet(dmeas),
                 rmeasure = Csnippet(rmeas),
                 rprocess = euler.sim(step.fun = Csnippet(sir.step_seas),
                                      delta.t = 1/40),
                 statenames = c("S", "I", "R", "H"),
                 paramnames = c("gamma", "mu", "theta", "Beta","omega", "phi","rho","popsize", "S.0", "I.0", "R.0"),
                 zeronames = c("H"),
                 initializer = function(params, t0, ...) {
                   fracs <- params[c("S.0", "I.0", "R.0")]
                   setNames(c(round(params["popsize"] * fracs/sum(fracs)), 0), c("S", "I", "R",
                                                                                 "H"))
                 },
                 toEstimationScale = toEstScale_seas,
                 fromEstimationScale = fromEstScale_seas,
                 params = param_vec_seas)



#a<- simulate(sir_seas, obs=TRUE, as.data.frame=TRUE)
#sir1_data<- data.frame(time=a$time[-1],cases=a$cases[-1])
#plot(sir1_data$cases, type="l")
#write.table(sir1_data, file = "~/Dropbox/AIC_study/AIC_study_git/sir_seas_data.txt" )


sir_box <- rbind(
  Beta = c(0.5,5),
  phi = c(0.001,0.1),
  rho = c(0.05,0.2)
)


sir_fixed_params <- param_vec[c("gamma","omega","theta","mu","popsize","S.0","I.0","R.0")]


##Fitting using iterated filtering with drawing 8 times from the sir_box
stew(file=paste("sir_seas_cross_",i,".rda", sep=""),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:8,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir_seas_cross,
        start=c(apply(sir_box,1,function(x)runif(1,min=x[1],max=x[2])),sir_fixed_params),
        Np=1000,
        Nmif=20,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(Beta=0.02,rho=0.02, phi=0.02))
    }
  })
},seed=1270401374,kind="L'Ecuyer")


##visualizing the diagnosistics of iterated filtering
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","phi", "rho", "theta")) %>%
  melt() %>%  subset(iteration>0) %>%  mutate(variable = factor(variable)) %>% rename(f=L1) -> t


##Make the plot
ggplot(t,aes(x=iteration,y=value,color=variable,group=f)) +
  geom_line() +
  guides(color=FALSE) +
  labs(x="MIF2 Iteration",y="") +
  facet_wrap(~variable,scales="free_y",ncol=2) +
  theme_bw()


## approximating the liklihood of every mif search outcome
stew(file=paste("sir_seas_cross_lik_",i,".rda", sep=""),{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:8,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir_seas_cross,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

best <- which.max(liks_global[,1])
sir_seas_cross_loglik <- round(liks_global[best,],3)
sir_seas_cross_aic<- aic(3,sir_seas_cross_loglik[1])


res<- c(sir1_aic=sir1_aic,
        sir2_aic=sir2_aic,
        sir1_cross_aic=sir1_cross_aic,
        sir2_cross_aic=sir2_cross_aic,
        sir_seas_aic=sir_seas_aic,
        sir_seas_cross_aic=sir_seas_cross_aic,
        sir1_cross_seas_aic=sir1_cross_seas_aic)

##Store AIC for each
coef_list[[i]] <-  res
}


names(coef_list) <- 1:10

##hoehle 2018-05-03: better imho - do it tibble style
tbl <- as.data.frame(t(bind_rows(coef_list)))
names(tbl) <- names(coef_list[[1]])

tbl %>% select(sir1_aic,sir2_cross_aic)

tbl %>% summarise( `1betterThan2_on1` = mean(sir1_aic < sir2_cross_aic),
                   `2betterThan1_on2` = mean(sir2_aic < sir1_cross_aic),
                   `seas_1betterThan2_on1` = mean(sir_seas_aic < sir1_cross_seas_aic),
                   `seas_2betterThan1_on2` = mean(sir1_aic < sir_seas_cross_aic))

