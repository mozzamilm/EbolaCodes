Theory of Infectious Disease Spillover at an Ecological Boundary: Impacts of Seasonality and Cross-boundary Movement
(Kaniz Fatema Nipa, Mozzamil Mohammed, Patrick R. Stephens, John M. Drake)

The numerical implementation of a two-patch deterministic and stochastic model and the results are presented in the above-mentioned paper.
The deterministic ODE system is coded in **Two-patch_Ebola_ODE_Model.R**.
The Continuous-Time Markov Chain (CTMC) Model is coded in **CTMC_Ebola_V02.R**.
The main script **MainDeterministicStochastic.Rmd** is where these functions are called and the manuscript results are produced.
**CTMC_SeasonPersistence.R** codes the effect of seasonality strength on the persistence probability of the disease.
**CTMC_Movement_BA_Persistence.R** codes the effect of movement strength of human settlement (Patch B) to the reservoir region (Patch A) on the persistence probability of the disease.
