% Compare peanut mobility solve time with and without opt.use_big_sparse.


close all;
clear

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

script_date = 'Apr 21, 2026';
fprintf('=== %s (%s) ===\n',mfilename,script_date);

if ~exist('P','var') || isempty(P)
    P = 100;
end
phi = 0.65;
rad = 1;
rng_seed = 210421;

geom_opt = struct('domain','boxed','phi',phi,'rad',rad, ...
    'min_gap',1e-3,'n_sweeps',200,'rng_seed',rng_seed, ...
    'visualise',false);
[q,meta] = random_discs_mc(P,geom_opt);

F = randn(P,2);
F = F - mean(F,1);
T = randn(P,1);

opt = get2Dparams(P,60,150);
opt.delta_pair = 0.2;
opt.N_peanut = 400;
opt.gmres_tol = 1e-8;
opt.maxit = 1000;
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
opt.reuse_pair_basis_by_sep = false;
opt.single_threaded = 0; 
opt.big_sparse_direct_u_corr = 1;  % 1 is a little less stable than 0 here

fprintf('Geometry: P=%d, phi=%.3f, close pairs=%d\n', ...
    P,meta.phi,count_close_pairs(q,opt.delta_pair,rad));

[~,~,~,~,~,pairs] = getEnhancedGrid(q,opt);
N_check = ceil(opt.a_c*opt.N_c);
labels = {'serial pair loop','big sparse'};
for k = 1:2
    opt.use_big_sparse = (k==2);

    if opt.use_big_sparse
        est = estimateMobPeanutBigSparseRamStokes( ...
            numel(q),opt.N_c,N_check,size(pairs,1),opt);
        fprintf('%-16s  sparse %.2f MB  build %.2f MB  peak %.2f MB\n', ...
            labels{k},est.estimated_sparse_MB,est.estimated_build_MB, ...
            est.estimated_peak_MB);
    end

    wall_timer = tic;
    [UW(:,k),sol] = solve_mob_peanut_enhanced(q,F,T,opt);
    wall_time = toc(wall_timer);
    
    fprintf(['%-16s  wall %.3fs  GMRES %.3fs  fmm/total %.3f  ', ...
        'it %d  final %.3e\n'], labels{k},wall_time,sol.solve_time.total, ...
        sol.solve_time.fmm/sol.solve_time.total,sol.it, ...
        sol.resvec(end));
end

