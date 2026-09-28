%CAPACITANCE_ON_CLUSTER Reproduce the hexagonal-pack capacitance demo but on 
% a randomized cluster with fixed interparticle distances.
%
% Running this script performs:
%   1) A peanut-compressed capacitance solve 
%   2) A peanut cap->elast two-way check on the same  geometry.
%
% Anna Broms, Sept 2026

clear;
close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;

fprintf('=== capacitance_on_hexagonal_pack ===\n\n');

%% Peanut-compressed capacitance solve on a large hexagonal packing
rng(8);

R = 2;
delta = 1e-3*R;
P=100;

% Solve capacitance on dense cluster
q = grow_cluster(P,delta,2,R);
P = length(q);
v_body = buildAlternatingVoltages(q,R);

% Set parameters and settings
N_c = 80;
N_f = 150; %used in the interpolation
opt = getLaplace2Dparams(P,R,N_c,N_f);
opt.N_cmap = 80;
opt.delta_pair = 0.4;
opt.Nclust = 100;
opt.N_peanut = 400;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_tol = 1e-8;
opt.debug = 0;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.compress_cmap = 0; % use low rank approximation of coarse-coarse map
opt.reuse_pair_basis_by_sep = 1;
opt.use_big_sparse = 0; 
opt.get_bndry_field = 1;
opt.refit = 1; 
opt.ellipse_constant = 1;

% Use the preferred reference-free reduced interpolation model.  The
% preparation call loads the compatible parameter-keyed MAT file; it does
% not silently retrain a missing model in batch mode.
opt.use_interpolation = 'reduced_noconst';
opt = prepareLaplaceCmapInterpolation(opt);

fprintf('Peanut-compressed capacitance solve\n');
fprintf('P=%d, delta=%.1e\n',P,delta);
tic;
[Q_peanut,sol_peanut] = solve_cap_peanut(q,v_body,opt);
t_peanut = toc;
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n\n', ...
    sol_peanut.it,sol_peanut.maxres,t_peanut);

opt.use_interpolation = 'none'; %not implemented for elastance
tic;
[Q_peanut_pre,sol_peanut_pre] = solve_cap_peanut(q,v_body,opt);
t_peanut = toc;
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n\n', ...
    sol_peanut_pre.it,sol_peanut_pre.maxres,t_peanut);

opt_tw = opt;
opt_tw.visualise_sol = 0;
opt_tw.gmres_tol = 1e-10;
opt_tw.refit = 0; %non-fft based
opt_tw.N_cmap = opt.N_c; 

fprintf('Peanut cap->elast two-way check\n');
tic;
[v_back_peanut,sol_elast_peanut] = solve_elast_peanut(q,Q_peanut,opt_tw);
t_tw = toc;
two_way_cap_elast_peanut = relerr(v_back_peanut,v_body);
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n', ...
    sol_elast_peanut.it,sol_elast_peanut.maxres,t_tw);
fprintf('  two-way error ||v_back-v_body||_inf / max(1,||v_body||_inf) = %.3e\n\n', ...
    two_way_cap_elast_peanut);

fprintf('Done.\n');

function e = relerr(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
