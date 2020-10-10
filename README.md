# Kelvin-Helmholtz Simulations

Project description:
> This paper discusses some aspects of 2D simulations of the Kelvin-Helmholtz instability. Note that the simulations mentioned here are incompressible, while you should run compressible cases. Use the paper to initiate your own simulations, and to make contact with the paper results.

Summary of the [reference paper](./schroeder2019reference.pdf):

# On reference solutions and the sensitivity of the 2D Kelvin-Helmholtz instability problem
## 1 Introduction
## 2 The Kelvin-Helmholtz instability problem
## 3 Self-organization and 2D turbulence
## 4 High-order divergence-free IMEX HDG methods
## 5 Computational studies
## 6 Computational studies of some possible sources for perturbations
## 7 Summary and conclusions
Contents of each section:

2. Problem setting, overview of quantities of interest (of which most important one is palinstrophy)
3. Van Groesen theory of self-organization applied to the Kelvin-Helmholtz instability problem (important: lemma 3.5)
4. Spatial (div-free HDG) and temporal (BDF2 IMEX) discretization schemes
5. Study quantities of interest for different Reynolds numbers; reliable reference results up to t = 200
6. Common sources of perturbations: kind of mesh used, linear system solver, numerical integration, compiler settings

## 8 Appendix
[Reference results](https://gitlab.gwdg.de/KHdata/KelvinHelmholtz) are available, can be useful to compare against.