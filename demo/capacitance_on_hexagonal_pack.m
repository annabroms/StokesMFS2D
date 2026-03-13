%CAPACITANCE_ON_HEXAGONAL_PACK Reproduce the hexagonal-pack capacitance demo.
%
% Running this script performs:
%   1) A peanut-compressed capacitance solve on a large hexagonal packing.
%   2) A 1B capacitance solve on a smaller hexagonal packing.
%
% Anna Broms, Mar 2026

clear;
close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;

fprintf('=== capacitance_on_hexagonal_pack ===\n\n');

%% Peanut-compressed capacitance solve on a large hexagonal packing
rng(8);

R = 2;
delta = 1e-3;

% Solve capacitance for hexagonal lattice
mode = 16;
q = hexagonal_lattice(delta,mode,R);
P = length(q);
v_body = buildAlternatingVoltages(q,R);

% Set parameters and settings
N_c = 60;
opt = getLaplace2Dparams(P,R,N_c);
opt.delta_pair = 0.2;
opt.Nclust = 100;
opt.N_peanut = 400;
opt.visualise_sol = 1;
opt.visualise_grid = 0;
opt.gmres_tol = 1e-8;
opt.debug = 0;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.compress_cmap = 0; % use low rank approximation of coarse-coarse map
opt.cmap_tol = 1e-8;
opt.reuse_pair_basis_by_sep = 1;

fprintf('Peanut-compressed capacitance solve\n');
fprintf('  hex mode=%d, P=%d, R=%.3f, delta=%.1e\n',mode,P,R,delta);
tic;
[Q_peanut,sol_peanut] = solve_cap_peanut(q,v_body,opt);
t_peanut = toc;
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n\n', ...
    sol_peanut.it,sol_peanut.maxres,t_peanut);

%% One-body capacitance solve on a smaller hexagonal packing
R = 2;

mode = 8;
delta = 1e-3;
q = hexagonal_lattice(delta,mode,R);
P = length(q);
v_body = buildAlternatingVoltages(q,R);

opt = getLaplace2Dparams(P,R);
opt.visualise_sol = 1;
opt.visualise_grid = 1;
opt.gmres_tol = 1e-8;
opt.debug = 0;
opt.use_fmm = true;
opt.N_c = 150;
opt.gmres_verbose = 0;

fprintf('1B capacitance solve\n');
fprintf('  hex mode=%d, P=%d, R=%.3f, delta=%.1e\n',mode,P,R,delta);
tic;
[Q_1B,sol_1B] = solve_cap_1B(q,v_body,opt);
t_1B = toc;
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n\n', ...
    sol_1B.it,sol_1B.maxres,t_1B);

fprintf('Done.\n');
