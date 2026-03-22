## Aligned particle experiments
Investigates iteration counts and boundary residuals for 
P particles in a line.

Reproduce using particle_line_mob.m and particle_line_elast.m.
### Mobility / Resistance

<p align="center">
  <img src="figs/mobility_iterations.png" width="45%" />
  <img src="figs/mobility_residual.png" width="45%" />
</p>

Mobility. Left: iterations. Right: relative residual.

<p align="center">
  <img src="figs/resistance_iterations.png" width="45%" />
  <img src="figs/resistance_residual.png" width="45%" />
</p>

Resistance. Left: iterations. Right: relative residual.
Long-range effects are dominant.

---


### Elastance / Capacitance

<p align="center">
  <img src="figs/elastance_iterations.png" width="45%" />
  <img src="figs/elastance_residual.png" width="45%" />
</p>

Elastance. Left: iterations. Right: relative residual.

<p align="center">
  <img src="figs/capacitance_iterations.png" width="45%" />
  <img src="figs/capacitance_residual.png" width="45%" />
</p>

Capacitance. Left: iterations. Right: relative residual.
Compare with resistance. Long-range effects are dominant.

**Note:**
- A higher resolution is used for Stokes than for Laplace, leading to the lower residuals.
- A larger pair-detection parameter `delta_pair` is needed for capacitance. With this choice, each particle is considered in a pair not only with its nearest neighbours, but also with particles beyond them. This does not appear to be necessary for capacitance problems on a hexagonal packing. Why?