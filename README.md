# Kelvin-Helmholtz Simulations

Project description:
> This paper discusses some aspects of 2D simulations of the Kelvin-Helmholtz instability. Note that the simulations mentioned here are incompressible, while you should run compressible cases. Use the paper to initiate your own simulations, and to make contact with the paper results.

How to run the stuff in `SimCode` (where you can replace "4" by the number of processors you want to use):
```bash
# One time
cd SimCode
mkdir SimData
$AMRVAC_DIR/setup.pl -d=2

# Each time mod_usr.t is updated
make -j 4

# Each simulation
mpirun -np 4 ./amrvac -i test.par
```

## TODO
- [ ] Reproduce some figures from the paper (assuming incompressible is easier than compressible)
    - [x] Translate the problem statement in [2.1](#2.1-Problem-statement) into a `.par` parameter file (see the base file, [settings.par](./SimCode/settings.par))
    - [x] Translate initial & boundary conditions into the `mod_usr.t` file (see the base file, [mod_usr.t](./SimCode/mod_usr.t))
    - [x] We also need to do some extra stuff with the vorticity, might want to check `amrvac/tests/hd/Karman_Vortex_2D` (there, the vorticity is also already being calculated, so might also be interesting for the following point)
    - [ ] Calculate the quantities of interest from simulation output 
        - See [here](http://amrvac.org/md_doc_mpiamrvac_nw.html) how to integrate this in the code
        - Palinstrophy seems to be the most important one
        - Maybe also consider numerical dissipation in section 5.4 from the paper
    - [ ] Read up on `&methodlist`, maybe make an overview below
- [ ] Run compressible cases of the simulations
- [ ] Make a presentation (for 9 / 16 December)

## Notes
### Overview of `&methodlist`
- `time_stepper`
- `flux_scheme`
- `limiter`


---

Summary of the [reference paper](./schroeder2019reference.pdf):

## On reference solutions and the sensitivity of the 2D Kelvin-Helmholtz instability problem
### 1 Introduction
Contributions of this paper:
- Reference results for 2D Kelvin-Helmholtz up to t = 200
- Even with state-of-the-art methods, the time of the final vortex pairing remains very difficult to predict
- Application of the theory for self-organization to Kelvin-Helmholtz

### 2 The Kelvin-Helmholtz instability problem
#### 2.1 Problem statement
- Equation (1) is what we need to simulate (but *compressible cases*, somehow)
- Boundary and initial conditions are also defined here
- Definition of Reynolds number

#### 2.2 Quantities of interest
This section gives 8 computable quantities which can be used to evaluate a given numerical scheme, check the paper for precise definitions. Kinetic energy (spectra), vorticity (thickness), enstrophy, palinstrophy, and the time intervals of the pairings are of interest.

### 3 Self-organization and 2D turbulence
The Van Groesen theory behind self-organization in 2D Navier-Stokes (with homogenous Dirichlet BC) is applied to the incompressible 2D Kelvin-Helmholtz problem (with mixed BC). Rather technical, uses some optimization and fancy analysis (didn't really read this section). The most important results are summarized in 2 remarks:
- The theory does not make any prediction about the position (in space nor time) of the final vortex
- 2D flows with high Reynolds numbers are very sensitive to perturbations

### 4 High-order divergence-free IMEX HDG methods
This section is less relevant to us than it seems. They use some finite element library in Matlab. As far as I know, MPI-AMRVAC only has finite differences and finite volumes.
#### 4.1 Space discretization
- Divergence-free HDC (Hybrid Discontinuous Galerkin) method based on H(div) finite elements (rectangular RT8 or triangular BDM8 elements, both using polynomials of order 8). 
- Some additional fancy stuff (projected jumps etc)
- For computing vorticity, enstrophy and palinstrophy: element-wise derivatives
- Careful with stabilization mechanisms, they might introduce a lot of numerical dissipation (resulting in lower accuracy) due to the sensitivity of the KH problem.

#### 4.2 Time discretization
- Linear multistep IMEX (BDF) method (in MPI-AMRVAC we will probably use IMEX-RK)
- Equation (13): result of a spatial discretization of incompressible Navier-Stokes
    - $Au$ is the Stokes part, this will be handled implicitly
    - $C(u)u$ is the convection part, this will be handled explicitly
- For linear system solving: sparse Cholesky with iterative refinement (not sure if we can choose the solver though)

### 5 Computational studies
- Sequence of square meshes with 16² to 256² elements
- Some vorticity figures (Figure 3) and energy spectra (Figure 4) are shown and discussed
- Then the quantities of interest from 2.2 are shown and discussed for Reynolds numbers `1E2`, `1E3` and `1E4`; for different mesh sizes:
    1. *Kinetic energy*
    2. *Enstrophy and vorticity thickness*
    3. *Palinstrophy*
    4. *Numerical dissipation* (not a problem quantity, but a general measure for numerical schemes)

### 6 Computational studies of some possible sources for perturbations
This section tries to find some explanations for why the results in literature vary so wildly, by a sensitivity analysis (i.e. identifying the reasons for small perturbations).
1. *Structured and unstructured triangular meshes*: rectangular, structured triangular or unstructured triangular meshes all seem to give different results
2. *Inaccurate solution of linear systems and numerical integration*: 
    - The iterative refinement for solving linear systems more accurately definitely seemed to have an impact on the end result
    - Impact of under-integration (because of a too low polynomial order used in the elements)
3. *Different compiler settings (FMA)*: even rounding modes seemed to have a non-negligible impact, the example of fused multiply-add is given.

### 7 Summary and conclusions
Contents of each section:

2. Problem setting, overview of quantities of interest (of which most important one is palinstrophy)
3. Van Groesen theory of self-organization applied to the Kelvin-Helmholtz instability problem (important: lemma 3.5)
4. Spatial (div-free HDG) and temporal (BDF2 IMEX) discretization schemes
5. Study quantities of interest for different Reynolds numbers; reliable reference results up to t = 200
6. Common sources of perturbations: kind of mesh used, linear system solver, numerical integration, compiler settings

### 8 Appendix
[Reference results](https://gitlab.gwdg.de/KHdata/KelvinHelmholtz) are available, can be useful to compare against.
