rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows
## Taken from v69i12.R
library(pomp)
require(doParallel)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
library(foreach)
cores <- 8
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)


##' ## More complex models.
##' ### Simple SIRS.
##' C snippets expressing the two faces of the measurement model.

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

##' Construct the pomp object and fill with simulated data.
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


#load the simulated data set - they were simulated form sir1 and sir 2 as explained below
read.table("sir1_data.txt") %>%
  rbind(data.frame(time=0,cases=NA)) %>%
  arrange(time) -> sir1_dat
head(sir1_dat)

read.table("sir2_data.txt") %>%
  rbind(data.frame(time=0,cases=NA)) %>%
  arrange(time) -> sir2_dat
head(sir2_dat)


# we start at t0=-10 so the system has time to equilibrate and reach the endemic level
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
 #This is how the data was created
     # a<- simulate(sir1,seed = 1916908L, obs=TRUE, as.data.frame=TRUE)
     # sir1_data<- data.frame(time=a$time[-1],cases=a$cases[-1])
     # plot(sir1_data$cases, type="l")
     # write.table(sir1_data, file = "~/Dropbox/AIC_study_michael/sir1_data.txt" )

######################################################################
##Fit model!
######################################################################

##One evaluation of the likelihood?
fit1 <- pfilter(sir1,params=param_vec,Np=1000)
logLik(fit1)


sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- c(param_vec["gamma"],param_vec["omega"],param_vec["theta"],param_vec["mu"],param_vec["popsize"], param_vec["S.0"], param_vec["I.0"], param_vec["R.0"])


##Fitting using iterated filtering with drawing 8 times from the sir_box
stew(file="sir1.rda",{
  
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


#visualizing the diagnosistics of iterated filtering
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>%   mutate(variable = factor(variable))->t
colnames(t)[4]<- "f"

  ggplot(t,aes(x=iteration,y=value,color=variable,group=f))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()


  # approximating the liklihood of every mif search outcome
stew(file="sir1_lik-%d.rda",{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir1,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")


liks_global
best <- which.max(liks_global[,1])
sir1_loglik <- round(liks_global[best,],3)
round(coef(mifs_global[[best]]),3)


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
S += dN[0] - dN[1] - dN[2] + dN[10];
I1 += dN[1] - dN[3] - dN[4];
I2 += dN[3] -  dN[5] - dN[6];
I3 += dN[5] -  dN[7] - dN[8];
R += dN[7] - dN[9] - dN[10];
H += dN[1];
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
  # a<- simulate(sir2,seed = 191469908L, obs=TRUE, as.data.frame=TRUE)
  #  sir2_data<- data.frame(time=a$time[-1],cases=a$cases[-1])
  #  plot(sir2_data$cases, type="l")
  # write.table(sir2_data, file = "~/Dropbox/AIC_study_michael/sir2_data.txt" )
 # # 

##One evaluation of the likelihood?
fit2 <- pfilter(sir2,params=param_vec2,Np=1000)
logLik(fit2)


sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- c(param_vec2["gamma"],param_vec2["theta"],param_vec2["mu"],param_vec2["popsize"],
                      param_vec2["S.0"], param_vec2["I1.0"], param_vec2["I2.0"], 
                      param_vec2["I3.0"], param_vec2["R.0"], param_vec2["omega"])


##Fitting using iterated filtering
stew(file="sir2.rda",{
  
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
# 
mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>%   mutate(variable = factor(variable))->t
colnames(t)[4]<- "f"

ggplot(t,aes(x=iteration,y=value,color=variable,group=f))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()

stew(file="sir2_lik-%d.rda",{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir2,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")



# 
liks_global
best <- which.max(liks_global[,1])
sir2_loglik <- round(liks_global[best,],3) 
round(coef(mifs_global[[best]]),3)

# 
# cbind(fit=coef(mifs_global[[best]]), true=param_vec2)





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
plot(sir1_cross)




sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- c(param_vec["gamma"],param_vec["theta"],param_vec["mu"],param_vec["popsize"], 
                      param_vec["S.0"], param_vec["I.0"], param_vec["R.0"],param_vec["omega"])


##Fitting using iterated filtering
stew(file="sir1_cross.rda",{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:cores,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar% {
      mif2(
        sir1_cross,
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



mifs_global %>%
  conv.rec(c("loglik","nfail","Beta","gamma", "mu", "theta")) %>%
  melt() %>% subset(iteration>17)  %>%mutate(variable = factor(variable))->t
colnames(t)[4]<- "f"

ggplot(t,aes(x=iteration,y=value,color=variable,group=f))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()


stew(file="sir1_lik_cross-%d.rda",{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir1_cross,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")


liks_global
best <- which.max(liks_global[,1])
sir1_cross_loglik <- round(liks_global[best,],3)
round(coef(mifs_global[[best]]),3)


cbind(fit=coef(mifs_global[[best]]), true=param_vec)



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


plot(sir2_cross)

sir_box <- rbind(
  Beta = c(0.5,10)
)

sir_fixed_params <- c(param_vec2["gamma"],param_vec2["theta"],param_vec2["mu"],param_vec2["popsize"],
                      param_vec2["S.0"], param_vec2["I1.0"], param_vec2["I2.0"], 
                      param_vec2["I3.0"], param_vec2["R.0"],param_vec2["omega"])


##Fitting using iterated filtering
stew(file="sir2_cross.rda",{
  
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
  melt() %>%   mutate(variable = factor(variable))->t
colnames(t)[4]<- "f"

ggplot(t,aes(x=iteration,y=value,color=variable,group=f))+
  geom_line()+
  guides(color=FALSE)+
  labs(x="MIF2 Iteration",y="")+
  facet_wrap(~variable,scales="free_y",ncol=2)+
  theme_bw()

stew(file="sir2_lik_cross-%d.rda",{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:cores,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(10, logLik(pfilter(sir2_cross,params=coef(mifs_global[[i]]),Np=1000)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")


# 
liks_global
best <- which.max(liks_global[,1])
sir2_cross_loglik<- round(liks_global[best,],3)
round(coef(mifs_global[[best]]),3)


#logliks of true model and other model 
#is the difference only due to Onte Carlo noise?
sir2_loglik
sir1_cross_loglik

sir1_loglik
sir2_cross_loglik