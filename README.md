Hernández-Aguirre, V.M., Rupakhety, R., Paolucci, R., Smerzini, C., Bessason, B., Erlingsson, S. (2026) 
A kinematic rupture generator for ground-motion simulations: Validation and scenarios in South Iceland.
Manuscript submitted for publication.

This repository includes a kinematic rupture generation framework useful for ensemble PBS. The generator integrates key elements from 
Song et al. (2014), Savran and Olsen (2020), and Aquib et al. (2025) to produce rupture realizations whose kinematic parameters are
not treated as independent ad hoc inputs but are generated with physically motivated structure and cross-parameter consistency. 

Slip distribution is simulated first, and then distributions of Vmax and Vr are simulated conditioned on slip, using a multivariate covariance model.

Function RuptureGeneration_Scenario.m applies the rupture generator step by step. All input parameters are defined in the first section: Input and preprocessing.
Function RuptureGeneration_Modify.m has prescribes slip distribution as a combination of a deterministic part at large wavelengths and a stochastic part at short wavelengths.

Supporting functions are provided in subfolder LibraryFunctions. Some functions are based on the work of Seok G. Song, and Martin Mai and collaborators.
