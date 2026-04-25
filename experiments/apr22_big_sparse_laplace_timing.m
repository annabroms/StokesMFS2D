% Compare Laplace peanut solve time with serial, precomputed, and streaming
% big-sparse close-pair matvecs.

close all;
clear

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (Apr 22, 2026) ===\n',mfilename);

P = 50; 
phi = 0.65;
rad = 2;
rng_seed = 220422;

geom_opt = struct('domain','boxed','phi',phi,'rad',rad, ...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',rng_seed, ...
    'visualise',false);
[q,meta] = random_discs_mc(P,geom_opt);

Q_body = randn(P,1);
Q_body = Q_body - mean(Q_body);

opt = getLaplace2Dparams(P,rad,40,80);
opt.delta_pair = 0.2*rad;
opt.N_peanut = 200;
opt.gmres_tol = 1e-7;
opt.get_bndry_field = 0;
opt.get_precomp_time = true;
opt.get_solve_time = true;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 1; 
opt.cmap = 1;
opt.reuse_pair_basis_by_sep = false; % change here if running with reused rotations. 
%                                       Then, the gain is not as large.
opt.use_fmm = true;
opt.rotation_mode = 'fft';
opt.single_threaded = 1; 

[~,~,~,~,~,pairs] = getEnhancedGrid(q,opt);
N_check = ceil(opt.a_c*opt.N_c);
fprintf('Geometry: P=%d, phi=%.3f, close pairs=%d\n', ...
    P,meta.phi,size(pairs,1));

variants = struct( ...
    'label',{'serial loop','big sparse precomputed','big sparse streaming'}, ...
    'use_big_sparse',{false,true,true}, ...
    'build_mode',{'auto','precomputed','streaming'});

for k = 1:3 %numel(variants)
    opt_run = opt;
    opt_run.use_big_sparse = variants(k).use_big_sparse;
    opt_run.lap_big_sparse_build_mode = variants(k).build_mode;

    if opt_run.use_big_sparse
        est = estimateLaplacePeanutBigSparseRam( ...
            numel(q),opt_run.N_c,N_check,size(pairs,1),opt_run);
        fprintf('%-22s sparse %.2f MB  build %.2f MB  peak %.2f MB\n', ...
            variants(k).label,est.estimated_sparse_MB, ...
            est.estimated_build_MB,est.estimated_peak_MB);
    end

    wall_timer = tic;
    [v_body(:,k),sol] = solve_elast_peanut(q,Q_body,opt_run); 
    wall_time = toc(wall_timer);

    fmm_fraction = get_fmm_fraction(sol);
    final_res = sol.resvec(end);
    build_time = get_big_sparse_build_time(sol);
    fprintf(['%-22s wall %.3fs  GMRES %.3fs  fmm/total %.3f  ', ...
        'big-build %.3fs  it %d  final %.3e\n'], ...
        variants(k).label,wall_time,sol.solve_time.total,fmm_fraction, ...
        build_time,sol.it,final_res);
end

function frac = get_fmm_fraction(sol)
frac = nan;
if isfield(sol,'solve_time') && isstruct(sol.solve_time) && ...
        isfield(sol.solve_time,'total') && sol.solve_time.total > 0 && ...
        isfield(sol.solve_time,'fmm')
    frac = sol.solve_time.fmm/sol.solve_time.total;
end
end

function t = get_big_sparse_build_time(sol)
t = 0;
if isfield(sol,'big_sparse_stats') && isstruct(sol.big_sparse_stats) && ...
        isfield(sol.big_sparse_stats,'build_time')
    t = sol.big_sparse_stats.build_time;
end
end
