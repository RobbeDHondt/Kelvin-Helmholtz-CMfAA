# Report of simulations

## Comparison of space discretization

### Finite differences

Investigated with `typeaverage = {"default", "roe", "arithmetic"}`. No difference between the three types of average visible in the simulations.

### TVDLF

+ limiter `ppm`: yields 4 vortices


## Meeting Jack 11/11/2020

### Preparation: some questions
* Error with `typeaverage`
* How to choose best discretization method, slope limiter?
  => suggestion: look at behaviour of enstrophy
* Show results for different Reynolds 
* Show attempt to reproduce simulation lecture yesterday (slide 23/52; density top 1.5 bottom 1, using `hll` and `mp5` and `ssprk5`, low Mach number flow)
* doublecheck: do we have the right settings? (compressible, not solving for energy...) 

### On the setting of our problem

Since the problem is nonlinear, the numerical behaviour can differ from simulation to simulation.
So it's hard to track what the benchmark solution is. One way to do so, is to run different schemes
and check which one runs the longest (and some other variables?) to set one scheme as a benchmark to
compare against. 

**(Robbe)** The way I understand it, is that we search for the "least artificial" flux scheme (on the AMRVAC site) and combine it with a high order timestepper (probably `fivestep ssprk5`) and a highly accurate slope limiter, and take the results of this simulation (possibly with a small timestep & gridsize) as benchmark result.
However, we don't *need* to make a benchmark result, we can also just compare a bunch of methods with one another.

**(Daniela)** From what I understood, there is no right results, therefore there is not a unique and right way to handle this. All depends on how we want to set our presentation, namely we can choose to have a benchmark and to compare the other simulation with, or we can show different simulations showing the differences between them. In case we want to have a benchmark, the *safest* way is probably to choose the most stable scheme. 

It's normal that at lower Reynolds number, less vortices are formed than the setup predicts. 

### On numerical artefacts

When we solve the incompressible case, some strange lines pop up in our visualisation of the vorticity.
This is a numerical artefact: it has something to do with sound speed and energy dissipation, but I did
not get it entirely so I hope Jack can give us some resource on that. 

**(Robbe)** What I have written about these strange lines:
- They are shocks (resulting from speeds higher than the sound speed). This corresponds to energy leaking from the system (?), and is therefore non-physical behavior. 
- In the incompressible case (no energy equation), you normally can't go above the sound speed (i.e. you can't add pressure to the system).
- Making the problem compressible (i.e. adding the energy equation) removes the issue of the sound speed.

**(Daniela)** I wrote down more or less the same thing. In the paper they have an approach which is able to better handle velocity and better resolve shocks. The difference between compressible and incompressible lies in the fact that in compressible cases, when a shock is formed, the system responds with and increase of density and pressure. This is physically not possible for incompressible cases and therefore we end up leaking energy from the system. 

The reason why the paper does not show these artefacts, has to do with the fact that they solve the problem
with a program specifically written for this problem. So they have a specific way to handle velocity which
makes the solving more accurate. AMRVAC on the other hand is a toolbox with predefined equations and 
modules to handle them, which makes it hard for the incompressible setting to be solved taking into
account all that may behave bad.

**(Daniela)** If I dod not misinterpret Jack's explanation, I think that the difference in the approach involves also the initial conditions. Basically, AMRVAC first solve the equation in a general form that can be applied to different system, and then it solve the specific problem taking into account all the parameters given in the .par file. In the paper, instead, they solved the equation taking into account all the parameters since the beginning, modifying the scheme step by step. 

Anyway, if we want to solve the compressible case, the only thing we should do is turn on the energy equation
in settings.par. When doing so, we should set the pressure in the mod_usr.t file. It can be done pretty
similar as to the test-problem in the amrvac-files (Dirichlet boundary conditions, uniform (?) initial 
condition).
**(Robbe)** Yeah, I also understood uniform pressure (using the `pint` variable as before) and Dirichlet BC.

### On discretization methods

Jack showed us a pretty nice simulation, which was made with `tvdlf` and `woodward` slope limiter. This
combination would be the most stable for computations, yet not necessarily the most accurate or correct
one. Yet as said before, it's hard to know what's accurate, except when you compare with physical experiments.

`minmod` is the most diffusive slope limiter. 
`ppm` (piecewise power method) on the other hand is 5th order in space, that's very sharp. 
Make sure that the order of your method for space discretization is not lower than the order of 
the method for time discretization. Apart from that, you can mess around as you like by combining different
methods. Hooray!

When using the `tvdmu` flux-scheme, we should also set `typeentropy`.
He also said something about an `entropyfix`.

### Possible settings to play with

* spin the other way around
* change pressure
* change velocity top, bottom
* sine pattern

### Grid refinement

It is possible to do the refinement according to the vorticity. Now it is only based on the four (or
three if we turn off energy equations) default variables, which is 'rho', 'v' (two dimensions) and 'p'.

### Visualisation

When using pytools, the output variables are interpolated op to second order between the grid points.
This may smoothen out sharp variations, which in some cases is not desirable. With yt, this can be 
done better.

### Good news

We do not need to worry about deadlines for having our simulations running on supercomputers: Jack is
so kind to offer us computing time on the desktop on his office, which has four cores. So this may 
be interesting in case we'd like to do deeper mesh refinement simulations.

### TO DO

For our next meeting, on wednesday the 18th of november, all of us think of a possible lay-out for the 
presentation. What are we going to investigate/analyse/discuss, in which order, with which emphasises.

Each of us also has a personal task: 

* Judith can catch up with how amrvac works and how our problem is set
* Daniela will look how we can switch to the compressible case
* Robbe will try to set the grid refinement for the vorticity
* Marie will try to visualize the enstrophy and palintrophy in Python using pytools and/or yt

**The end**
