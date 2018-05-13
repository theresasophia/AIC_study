# AIC_study
Short simulation study to answer question 1.2 of  Reviewer 1 concerning wether the AIC can detect model miss-specifications in the transmission model or not. For this we simulate from a simple SIRS model with birth and death and exponentially distributed infectious period with mean 1/gamma. We let the system equilibrate so that it reaches the endemic state. Moreover we simulate from a similar  SIRS model with the only difference that now the infectious period is Gamma(3,gamma*3) distributed. Both models have the same mean of infectious period 1/gamma. We treat the simulations (sir1_data.txt and sir2_data.txt) as data and fit them to the respective models. We only fit the parameter \beta and keep all other parameters fixed at their true value. We do this because this is most similar to what we did in the manuscript.

