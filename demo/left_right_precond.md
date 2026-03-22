## Left vs Right 1-Body Preconditioning (Resistance+mobility)

The results of the demo left_right_precond.m:

The 2D Stokes resistance and mobility problems are solved on randomized clusters of 
P disks with separation parameter delta. Left and right 1-body preconditioning, using the same discretisation, are compared in terms of GMRES iterations, relative boundary residuals, and force/torque mismatch.

### Resistance
<p align="center">
  <img src="figs/left_right_iterations.png" width="100%" />  
</p>

<p align="center">
  <img src="figs/left_right_residuals.png" width="100%" />  
</p>

<p align="center">
  <img src="figs/left_right_ft_error.png" width="100%" />  
</p>

### Mobility

Larger difference in the number of iterations for small delta.

<p align="center">
  <img src="figs/mob_left_right_iterations.png" width="100%" />  
</p>

<p align="center">
  <img src="figs/mob_left_right_residuals.png" width="100%" />  
</p>

