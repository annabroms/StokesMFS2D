% Compare mobility big-sparse solves with sparse_map_coarse off and on.

close all;
clear

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (Apr 24, 2026) ===\n',mfilename);

compare_case('Case 1',50,0.65,240424,80,60,200,true);
fprintf('\n');
compare_case('Case 2',100,0.65,240425,60,150,200,false);

function compare_case(label,P,phi,rng_seed,N_c,N_f,maxit, ...
    report_output_comparison)
q = build_geometry(P,phi,rng_seed);
opt = make_options(P,N_c,N_f);
opt.maxit = maxit;

F = randn(P,2);
F = F - mean(F,1);
T = randn(P,1);

[~,~,~,~,~,pairs] = getEnhancedGrid(q,opt);
N_check = ceil(opt.a_c*opt.N_c);

fprintf('%s: P=%d, phi=%.2f, close pairs=%d\n', ...
    label,P,phi,size(pairs,1));
runs = repmat(struct('UW',[],'sol',struct()),2,1);
for sparse_map_coarse = [false true]
    opt_run = opt;
    opt_run.sparse_map_coarse = sparse_map_coarse;
    est = estimateMobPeanutBigSparseRamStokes( ...
        P,opt_run.N_c,N_check,size(pairs,1),opt_run);
    wall_timer = tic;
    [UW,sol] = solve_mob_peanut_enhanced(q,F,T,opt_run);
    wall_time = toc(wall_timer);
    runs(1 + sparse_map_coarse) = struct('UW',UW,'sol',sol);
    fprintf(['  sparse_map_coarse=%d  wall=%.3fs  solve=%.3fs  ', ...
        'fmm/total=%.3f  matrix=%.2f MB  peak=%.2f MB  ', ...
        'it=%d  final=%.3e\n'], ...
        sparse_map_coarse,wall_time,sol.solve_time.total, ...
        sol.solve_time.fmm/sol.solve_time.total, ...
        est.big_sparse_matrix_bytes/1024^2, ...
        est.big_sparse_peak_bytes/1024^2,sol.it,sol.resvec(end));
end

if report_output_comparison
    fprintf(['  relerr(UW on/off)=%.3e  resvec diff=%.3e  ', ...
        'it=%d/%d\n'], ...
        relerr(runs(2).UW,runs(1).UW), ...
        norm(runs(2).sol.resvec-runs(1).sol.resvec,inf), ...
        runs(1).sol.it,runs(2).sol.it);
end
end

function q = build_geometry(P,phi,rng_seed)
geom_opt = struct('domain','boxed','phi',phi,'rad',1, ...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',rng_seed, ...
    'visualise',false);
[q,~] = random_discs_mc(P,geom_opt);
end

function opt = make_options(P,N_c,N_f)
opt = get2Dparams(P,N_c,N_f);
opt.delta_pair = 0.2;
opt.N_peanut = 400;
opt.gmres_tol = 1e-8;
opt.get_bndry_field = 0;
opt.get_precomp_time = true;
opt.get_solve_time = true;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.use_big_sparse = true;
opt.reuse_pair_basis_by_sep = false;
opt.parallel_precomp = false;
opt.big_sparse_direct_u_corr = true;
opt.mob_big_sparse_build_mode = 'streaming';
opt.mob_big_sparse_chunk_pairs = 8;
opt.single_threaded = 1; 
end

function e = relerr(a,b)
e = norm(a(:)-b(:),inf)/max(1,norm(b(:),inf));
end
