function opt = getLaplace2Dparams()
%GETLAPLACE2DPARAMS Default parameters for scalar Laplace MFS in 2D.
%
% See also: solve_cap_1B, solve_cap_2B, solve_cap_peanut.
%
% Anna Broms, Mar 2026

% Geometry scale: physical particle radius.
opt.rad = 2;

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

opt.pc = 1;
opt.delta_pair = (opt.rad-opt.Rp_c)^2/opt.Rp_c;
opt.beta = 0.3;
opt.Nclust = 100;

opt.N_peanut = 0;
opt.precomp = 1;
opt.cmap = 0;

opt.P = [];
opt.use_fmm = true;
opt.project_charge = false;

end
