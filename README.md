# StokesMFS2D
Method-of-fundamental-solutions (MFS) solvers for 2D Stokes and Laplace boundary value problems on geometries of close-to-touching disks in free space. This research code demonstrates pair-based preconditioning. In the _2B solvers, this is implemented using a fine set of sources for each pair of disks. Solvers with the _peanut extension utilise compression of this two-body fine discretisation to an equivalent coarse set of proxy nodes within each body. Greater gains are achieved from pair-based preconditioning for Stokes BVPs than for Laplace, owing to the stronger singularities in close-to-touching geometries.

## Repository overview
- `mobility/`
  Stokes mobility solvers (input force/torque, output rigid-body velocities), including:
  `solve_mob_2B_enhanced`, `solve_mob_peanut_enhanced`.
- `resistance/`
  Stokes resistance solvers (input rigid-body velocities, output force/torque), including:
  `solve_res_2B_enhanced`, `solve_res_peanut_enhanced`.
- `laplace/`
  Scalar Laplace solvers for capacitance and elastance problems, where the capacitance formulation corresponds to a modified BVP.[^stein-barnett-2022]
  - Capacitance (`solve_cap_1B`, `solve_cap_2B`, `solve_cap_peanut`):
    prescribe per-body voltages `v_body`, solve for per-body net charges `Q_body`.
  - Elastance (`solve_elast_1B`, `solve_elast_2B`, `solve_elast_peanut`):
    prescribe per-body net charges `Q_body`, solve for per-body voltages `v_body`.
  - In the Laplace solvers, `opt.rad ~= 1` is used to avoid the unit
    logarithmic-capacity case in 2D.
- `geometry/`
  Grid generation and close-pair enhancement tools.
- `common/`
  Shared kernels for resistance and mobility: basis assembly, transformations, and helper utilities.
- `gmres/`
  GMRES utilities.
- `demo/`
  Practical comparison scripts and usage examples, with two primary entry points:
  - `test_mob_res.m`: compares Stokes mobility/resistance solver families (1B, 2B, peanut), including two-way checks.
  - `test_cap_elast.m`: compares Laplace capacitance/elastance solver families (1B, 2B, peanut), including two-way checks.
- `experiments/`
  Some research experiments long-range preconditioning.
- `tools/`
  Plotting and helper routines.

## Self tests
Most core solver/operator functions contain a built-in self test.  
Run the function with no input arguments, for example:
- `solve_mob_2B_enhanced()`
- `solve_res_peanut_enhanced()`
- `solve_cap_2B()`
- `solve_elast_peanut()`

## Highlighted demos
- `demo/test_mob_res.m`
  - Stokeslet-only enhanced comparisons (`solve_mob_2B_enhanced`, `solve_mob_peanut_enhanced`, `solve_res_2B_enhanced`, `solve_res_peanut_enhanced`).
  - Optional 1B image-based comparisons (`solve_mob_1B`, `solve_res_1B`) that still use mixed source types.
  - Demos support line, dumbbell, and cluster particle layouts.
  - Two-way check means: solve mobility from `(F,T)` to get `(U,W)`, then solve resistance using that `(U,W)` and compare recovered `(F,T)` to the original input (and vice versa).
- `demo/test_cap_elast.m`
  - Capacitance: prescribed `v_body`, solve for net charge `Q_body`.
  - Elastance: prescribed net charge `Q_body`, solve for `v_body`.
  - Two-way check means: solve capacitance from `v_body` to get `Q_body`, then solve elastance using that `Q_body` and compare recovered `v_body` to the original input (and vice versa).

## Source-type note
- The `dev_image_based/` folders under `mobility/`, `resistance/`, and `common/` contain older and not as polished helper functions for Stokes close interaction resolution based on multiple source types.
- Outside these `dev_image_based/` paths, the Stokes solver paths used in the main comparisons are Stokeslet-only.


## Dependencies
- [fmm2d](https://github.com/flatironinstitute/fmm2d) for fast 2D kernel summation.


[^stein-barnett-2022]:
  Stein, D. B., and Barnett, A. H. (2022).  
  *Quadrature by fundamental solutions: kernel-independent layer potential evaluation for large collections of simple objects.*  
  Advances in Computational Mathematics (ACOM). The capacitance formulation here follows the lemma introducing the modified Laplace exterior 2D BVP.


<!--
## Previous README (commented for comparison)
# StokesMFS2D
Solves the Stokes resistance and mobility BVPs in 2D using the method of fundamental solutions.

## Close-to-touching disks
We investigate preconditioning ideas based on pair corrections and peanut compression.
Test examples in demo/ and by running solver functions without arguments.


***To do***
- Introduce effective mapping \mu -> \lambda
- Compare timings w/o peanut compression
- Speedup of transformation \mu -> \lambda, beta: determine in parallel
- Write direct Stokeslet eval in fortran / C. Same for pot. dipole and stresslet.
- Later: Make sure Krylov precond can be switched on to compare #GMRES iterations, most interesting with time-stepping (mobility).

## Other geometries
The `experiments/` folder contains test cases for a variety of non-circular geometries:
- Stars
- Ellipses
- Objects with corners

## Dependencies

- [fmm2d](https://github.com/flatironinstitute/fmm2d) — Fast Multipole Method library for 2D kernels
- [Chebfun](https://www.chebfun.org/) — Needed for testing with non-circular geometries
-->
