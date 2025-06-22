# StokesMFS2D
Solves the Stokes resistance and mobility BVPs in 2D using the method of fundamental solutions.

## Close-to-touching discs
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