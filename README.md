Pseudo-dynamic rupture generator (PDRG) based on 1-point and 2-point statistics of earthquake source parameters (Hernández-Aguirre et al., 2026).

The adopted simulation method is similar to that of Song et al. (2014), with the difference that Cholesky factorization is done on the covariance matrix conditioned on the slip, defined previously. It includes components from various previous models, such as Mai et al. (2002), Mai et al. (2005), Aquib et al. (2015).

The two main codes are:

RuptureGeneration_Scenario.m: To create a new source model, starting from the slip.
RuptureGeneration_Modify.m: Merges a deterministic slip distribution with a stochastic part at short wavelengths. 

Supporting functions are inside the folder SongModified and LibraryFunctions.

