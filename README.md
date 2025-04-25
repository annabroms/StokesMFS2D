# StokesMFS2D
Solves the Stokes resistance and mobility BVPs in 2D for close to touching disks. Investigation of preconditioning ideas based on pair corrections and peanut compression (see document on overleaf).
Test example in demo/ 

Dependencies: fmm2d https://github.com/flatironinstitute/fmm2d

**To do**
- Compare timings w/o peanut compression
- Speedup of precompute step: determine pair corrections in parallel
- Speedup of transformation \mu -> \lambda, beta: determine in parallel
- Write direct Stokeslet eval in fortran / C. Same for pot. dipole and stresslet.
- Make sure Krylov precond can be switched on to compare #GMRES iterations, most interesting with time-stepping (mobility). 
