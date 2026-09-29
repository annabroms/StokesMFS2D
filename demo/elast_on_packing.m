%LARGE_ELAST_ON_PACKING Reproduce the capacitance / elastance demo 
% on a close random packing but run it only for elastance so that we can
% crank up the number of disks
%
%
% Anna Broms, Sept 2026

clear;
close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;

fprintf('=== elastance_on_close_random_packing ===\n\n');

%% Interpolated peanut-compressed elastance solve on a close random packing
rng(8);

R = 2;
delta = 1e-3*R;

% Solve capacitance on dense cluster
P = 1000;
opt_mc.phi = 0.65;
opt_mc.rad = R;
q = random_discs_mc(P,opt_mc);
P = length(q);
Q_body = buildAlternatingVoltages(q,R);

% Set parameters and settings
N_c = 80;
N_f = 150; %used in the interpolation
opt = getLaplace2Dparams(P,R,N_c,N_f);
opt.use_interpolation = 'full';
opt.N_cmap = 80;
opt.delta_pair = 0.4;
opt.N_peanut = 400;
opt.visualise_grid = 0;
opt.gmres_tol = 1e-10;
opt.debug = 0;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.compress_cmap = 0; % use low rank approximation of coarse-coarse map
opt.reuse_pair_basis_by_sep = 1;
opt.use_big_sparse = 0; 
opt.get_bndry_field = 1;
opt.refit = 1; 
opt.ellipse_constant = 1;
opt.Nclust = 150;
opt.volt_charge_interp_tol = 1e-6;
opt.interpolation_tol = 1e-5;
opt.tikhonov_tol = 1e-11; %1e-11;

opt_tw = opt;
opt_tw.visualise_sol = 0;
opt_tw.gmres_tol = 1e-10;
opt_tw.refit = 0; % Fourier rotation on equal coarse grids
opt_tw.N_cmap = opt.N_c;
opt_tw.interpolation_model_file = '';
opt_tw = prepareLaplaceCmapInterpolation(opt_tw,'elastance');
Q_body = Q_body-mean(Q_body); 

tic;
[v_back_peanut,sol_elast_peanut] = solve_elast_peanut(q,Q_body,opt_tw);
t_tw = toc;
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n', ...
    sol_elast_peanut.it,sol_elast_peanut.maxres,t_tw);


fprintf('Done.\n');

function e = relerr(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
