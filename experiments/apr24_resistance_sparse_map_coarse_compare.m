% Compare resistance big-sparse solves with res_sparse_map_coarse and
% res_smaller_ucorr toggled on and off.

close all;
clear

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (Apr 24, 2026) ===\n',mfilename);

compare_case('Case 1',30,0.65,240424,120,60,250,true);
fprintf('\n');
compare_case('Case 2',60,0.65,240425,60,150,1000,false);

function compare_case(label,P,phi,rng_seed,N_c,N_f,maxit, ...
    report_output_comparison)
q = build_geometry(P,phi,rng_seed);
opt = make_options(P,N_c,N_f);
opt.maxit = maxit;

U = randn(P,2);
W = randn(P,1);

[~,~,~,~,~,pairs] = getEnhancedGrid(q,opt);
N_check = ceil(opt.a_c*opt.N_c);

fprintf('%s: P=%d, phi=%.2f, close pairs=%d\n', ...
    label,P,phi,size(pairs,1));
runs = repmat(struct('FT',[],'sol',struct()),2,2);
for sparse_map_coarse = [false true]
    for res_smaller_ucorr = [false true]
        opt_run = opt;
        opt_run.res_sparse_map_coarse = sparse_map_coarse;
        opt_run.res_smaller_ucorr = res_smaller_ucorr;
        est = estimateResPeanutBigSparseRam( ...
            P,opt_run.N_c,N_check,size(pairs,1),opt_run);
        wall_timer = tic;
        [FT,sol] = solve_res_peanut_enhanced(q,U,W,opt_run);
        wall_time = toc(wall_timer);
        runs(1 + sparse_map_coarse,1 + res_smaller_ucorr) = struct( ...
            'FT',FT,'sol',sol);
        fprintf(['  res_sparse_map_coarse=%d res_smaller_ucorr=%d  wall=%.3fs  ', ...
            'solve=%.3fs  fmm/total=%.3f  matrix=%.2f MB  ', ...
            'peak=%.2f MB  it=%d  final=%.3e\n'], ...
            sparse_map_coarse,res_smaller_ucorr,wall_time,sol.solve_time.total, ...
            sol.solve_time.fmm/sol.solve_time.total, ...
            est.big_sparse_matrix_bytes/1024^2, ...
            est.big_sparse_peak_bytes/1024^2,sol.it,sol.resvec(end));
    end
end

if report_output_comparison
    for sparse_map_coarse = [false true]
        baseline = runs(1 + sparse_map_coarse,1);
        smaller = runs(1 + sparse_map_coarse,2);
        fprintf(['  res_smaller_ucorr relerr(FT) at res_sparse_map_coarse=%d: ', ...
            '%.3e  resvec diff=%.3e  it=%d/%d\n'], ...
            sparse_map_coarse,relerr(smaller.FT,baseline.FT), ...
            norm(smaller.sol.resvec-baseline.sol.resvec,inf), ...
            baseline.sol.it,smaller.sol.it);
    end
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
opt.res_big_sparse_u_corr_mode = 'combined';
opt.res_smaller_ucorr = false;
opt.res_big_sparse_chunk_pairs = 8;
opt.single_threaded = 1; 
end

function e = relerr(a,b)
e = norm(a(:)-b(:),inf)/max(1,norm(b(:),inf));
end
