# DataDriven-FlowMaps
Data-Driven Flow-Map Models for Data-Efficient Discovery of Dynamics and Fast Uncertainty Quantification of Biological and Biochemical Systems

## Abstract
Computational models are increasingly used to investigate and predict the complex dynamics of biological and biochemical systems. Nevertheless, governing equations of a biochemical system may not be (fully) known, which would necessitate learning the system dynamics directly from, often limited and noisy, observed data. On the other hand, when expensive models are available, systematic and efficient quantification of the effects of model uncertainties on quantities of interest can be an arduous task. This paper leverages the notion of flow-map (de)compositions to present a framework that can address both of these challenges via learning data-driven models useful for capturing the dynamical behavior of biochemical systems. Data-driven flow-map models seek to directly learn the integration operators of the governing differential equations in a black-box manner, irrespective of structure of the underlying equations. As such, they can serve as a flexible approach for deriving fast-to-evaluate surrogates for expensive computational models of system dynamics, or, alternatively, for reconstructing the long-term system dynamics via experimental observations. We present a data-efficient approach to data-driven flow-map modeling based on polynomial chaos Kriging. The approach is demonstrated for discovery of the dynamics of various benchmark systems and a co-culture bioreactor subject to external forcing, as well as for uncertainty quantification of a microbial electrosynthesis reactor. Such data-driven models and analyses of dynamical systems can be paramount in the design and optimization of bioprocesses and integrated biomanufacturing systems.

## Case Study 1
The first benchmark is the Morris-Lecar system which describes neuronal excitability. Using limited data, we predict the long term dynamics of the system under varying values of injected current.

## Case Study 2
In this case we attempt to reconstruct the dynamics of a chaotic system described by the well-known Lorenz equations.  The dynamics vastly differ depending on the initial conditions and model parameters, yet, using limited training samples, the limit cycles and other main characteristics are adequately captured.

## Case Study 3
We demonstrate the ability of PCK-based data-driven flow-map models to learn the transient behavior of a co-culture system with variable inputs.  Using a non-linear bioreactor model, we obtain observations that are corrupted with noise to emulate real data collection. In that case, the PCK model yields predictive error/confidence bounds on one-step ahead predictions which can be utilized in various ways to quantify uncertainty for trajectory generation

## Case Study 4
We examine the utility of data-driven flow-maps for the UQ of a Microbial Electrosynthesis (MES) bioreactor using a high-fidelity computational model that is subject to uncertainty in model parameters and initial conditions.  The data-driven flow map can be used to predict trajectories under various conditions and parameters much faster compared to the high-fidelity integrator. This renders the problem of forward and inverse UQ computationally tractable
