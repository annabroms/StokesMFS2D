%CAPACITANCE_ON_HEXAGONAL_PACK Reproduce the hexagonal-pack capacitance demo.
%
% Running this script performs:
%   1) A peanut-compressed capacitance solve on a large hexagonal packing.
%   2) A peanut cap->elast two-way check on the same large geometry.
%   3) A 1B capacitance solve on a smaller hexagonal packing.
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
rings = 9; %8 rings around the center disk 
q = hexagonal_lattice(delta,rings,R);
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
opt.reuse_pair_basis_by_sep = 1;

fprintf('Peanut-compressed capacitance solve\n');
fprintf('  hex rings=%d, P=%d, R=%.3f, delta=%.1e\n',rings,P,R,delta);
tic;
[Q_peanut,sol_peanut] = solve_cap_peanut(q,v_body,opt);
t_peanut = toc;
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n\n', ...
    sol_peanut.it,sol_peanut.maxres,t_peanut);

opt_tw = opt;
opt_tw.visualise_sol = 0;

fprintf('Peanut cap->elast two-way check\n');
tic;
[v_back_peanut,sol_elast_peanut] = solve_elast_peanut(q,Q_peanut,opt_tw);
t_tw = toc;
two_way_cap_elast_peanut = relerr(v_back_peanut,v_body);
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n', ...
    sol_elast_peanut.it,sol_elast_peanut.maxres,t_tw);
fprintf('  two-way error ||v_back-v_body||_inf / max(1,||v_body||_inf) = %.3e\n\n', ...
    two_way_cap_elast_peanut);

%% One-body capacitance solve on a smaller hexagonal packing
R = 2;

rings = 4; %four rings around the center disk 
delta = 1e-3;
q = hexagonal_lattice(delta,rings,R);
P = length(q);
v_body = buildAlternatingVoltages(q,R);

N_c = 150; %use fine grid for more fair accuracy comparison
opt = getLaplace2Dparams(P,R,N_c);
opt.visualise_sol = 1;
opt.visualise_grid = 1;
opt.gmres_tol = 1e-8;
opt.debug = 0;
opt.use_fmm = true;
opt.gmres_verbose = 2; %0: silent

fprintf('1B capacitance solve\n');
fprintf('  hex rings=%d, P=%d, R=%.3f, delta=%.1e\n',rings,P,R,delta);
tic;
[Q_1B,sol_1B] = solve_cap_1B(q,v_body,opt);
t_1B = toc;
fprintf('  GMRES it=%d, maxres=%.3e, time=%.2f s\n\n', ...
    sol_1B.it,sol_1B.maxres,t_1B);

fprintf('Done.\n');

function e = relerr(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
