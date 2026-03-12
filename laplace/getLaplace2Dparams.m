function opt = getLaplace2Dparams(P,R)
%GETLAPLACE2DPARAMS Default parameters for scalar Laplace MFS in 2D.
%
% Input: P - number of particles.
%        R - radius
%
% See also: solve_cap_1B, solve_cap_2B, solve_cap_peanut.
%
% Anna Broms, Mar 2026

opt.rad = R;
opt.P = P; 

N_c = 80;
opt.N_c = N_c;

opt.a_c = 1.2;
opt.a_f = 1.2;

tol = 1e-12;
sep = (1/N_c)*log(1/tol);
opt.Rp_c = opt.rad*max([1-sep,0.01]);

N_f = 150;
opt.N_f = N_f;
sep = (1/N_f)*log(1/tol);
opt.Rp_f = opt.rad*max([1-sep,0.01]);

opt.pc = 1; %Do pair corrections? Is this field still active?
%opt.delta_pair = (opt.rad-opt.Rp_c)^2/opt.Rp_c;
opt.delta_pair = 0.2; 
opt.beta = 0.3;
opt.Nclust = 100;

% Only relevant with pair corrections
opt.N_peanut = 400; %Number of nodes on peanut separation surface
opt.precomp = 1;
opt.cmap = 1;
opt.show_counter = 1; % show progress for pair compression

opt.use_fmm = true;
opt.debug = 0;
opt.gmres_tol = 1e-7;
opt.project_charge = false;
opt.gmres_verbose = 0; % 0=silent, 1=final summary, 2=per-iteration

opt.visualise_grid = 1; %draw source and collocation points at setup stage.

end
