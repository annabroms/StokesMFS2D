clear;
close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Peanut Weighting Comparison On 7-Disk Hex Lattice (Mar 26, 2026) ===\n');

rng(8);

rad = 1;
delta = 1e-3;
q = hexagonal_lattice(delta,1,rad);
P = numel(q);
assert(P == 7,'Expected 7 disks from hexagonal_lattice(delta,1,rad).');

F = randn(P,2);
F = F - mean(F,1);
T = randn(P,1);
U = randn(P,2);
W = randn(P,1);

opt = get2Dparams(P,150,150);
opt.rad = rad;
opt.delta_pair = 0.2;
opt.N_peanut = 400;
opt.gmres_tol = 1e-8;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.surface_error_mode = 'rel';
opt.reuse_pair_basis_by_sep = true;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.get_bndry_field = 1;

fprintf('Geometry: P=%d, delta=%.2e, delta_pair=%.2e, N_c=%d, N_f=%d, N_peanut=%d\n', ...
    P,delta,opt.delta_pair,opt.N_c,opt.N_f,opt.N_peanut);

modes = struct( ...
    'name',{'none','column','left','both'}, ...
    'column_weight',{false,true,false,true}, ...
    'left_weight',{false,false,true,true});

results = repmat(struct( ...
    'name','', ...
    'column_weight',false, ...
    'left_weight',false, ...
    'UW',[], ...
    'FT',[], ...
    'mob_sol',struct(), ...
    'res_sol',struct(), ...
    'mob_time',NaN, ...
    'res_time',NaN),numel(modes),1);

for it = 1:numel(modes)
    opt_run = opt;
    opt_run.column_weight = modes(it).column_weight;
    opt_run.left_weight = modes(it).left_weight;

    results(it).name = modes(it).name;
    results(it).column_weight = opt_run.column_weight;
    results(it).left_weight = opt_run.left_weight;

    tic;
    [results(it).UW,results(it).mob_sol] = solve_mob_peanut_enhanced(q,F,T,opt_run);
    results(it).mob_time = toc;

    %tic;
    %[results(it).FT,results(it).res_sol] = solve_res_peanut_enhanced(q,U,W,opt_run);
    %results(it).res_time = toc;

   % clear mex;
end

fprintf('\nMobility peanut solver:\n');
fprintf('  %-8s %6s %6s %10s %8s %14s %14s\n', ...
    'mode','col','left','time(s)','it','rel_res','abs_res');
for it = 1:numel(results)
    fprintf('  %-8s %6d %6d %10.2f %8d %14.3e %14.3e\n', ...
        results(it).name,results(it).column_weight,results(it).left_weight, ...
        results(it).mob_time,results(it).mob_sol.it, ...
        results(it).mob_sol.rel_res,results(it).mob_sol.abs_res);
end

fprintf('\nResistance peanut solver:\n');
fprintf('  %-8s %6s %6s %10s %8s %14s\n', ...
    'mode','col','left','time(s)','it','rel_res');
for it = 1:numel(results)
    fprintf('  %-8s %6d %6d %10.2f %8d %14.3e\n', ...
        results(it).name,results(it).column_weight,results(it).left_weight, ...
        results(it).res_time,results(it).res_sol.it, ...
        results(it).res_sol.rel_res);
end

fprintf('\nComparisons vs none:\n');
for it = 2:numel(results)
    fprintf('  %-8s relerr(UW,%s)=%.3e   relerr(FT,%s)=%.3e\n', ...
        results(it).name,results(1).name,relerr(results(it).UW,results(1).UW), ...
        results(1).name,relerr(results(it).FT,results(1).FT));
end

fprintf('\nCross comparisons:\n');
fprintf('  both   vs column: relerr(UW)=%.3e   relerr(FT)=%.3e\n', ...
    relerr(results(4).UW,results(2).UW),relerr(results(4).FT,results(2).FT));
fprintf('  both   vs left  : relerr(UW)=%.3e   relerr(FT)=%.3e\n', ...
    relerr(results(4).UW,results(3).UW),relerr(results(4).FT,results(3).FT));

function e = relerr(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
