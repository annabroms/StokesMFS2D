close all; 
clear;

% Random seeds
geom_seed = 1;
load_seed = 11;

% Geometry
P = 20;
rad = 1;
domain = 'boxed';
phi = 0.65;
delta = 1e-3; 
n_sweeps = 30;
visualise_geometry = false;
N_peanut = 400; 
N_c = 100; 
N_f = 40; %activates ellipse segment earlier

visualise = 1; 
debug = 0;

delta_pair = 0.2;
%delta_pair = 0.01;



opt = get2Dparams(1,N_c,N_f);
opt.delta_pair = delta_pair;
opt.N_peanut = N_peanut;
opt.visualise_sol = visualise;
opt.debug = debug;
opt.cmap = 1; 
opt.beta = 0.3;
opt.reuse_pair_basis_by_sep = 0;
opt.rotation_mode = 'oversampled_fft';
opt.rotation_oversample = 8;
opt.visualise_sol = visualise;
opt.gmres_tol = 1e-8;
opt.debug = debug;
%opt.Nclust = 200; 
%opt.beta = 0.5;

%% Check that parameters make sense
report = test_pair_corrections_stokes(opt,@solve_mob_peanut_enhanced,@solve_res_peanut_enhanced);

%% Setup test and solve



geom_opt = struct();
geom_opt.domain = domain;
geom_opt.phi = 0.65;
geom_opt.rad = 1;
geom_opt.min_gap = delta;
geom_opt.n_sweeps = n_sweeps;
geom_opt.rng_seed = geom_seed;
geom_opt.visualise = visualise_geometry;

[q,geom_meta] = random_discs_mc(P,geom_opt);

q = grow_cluster(P,delta,2); 

U = rand(P,2); W = rand(P,1); 
%W = zeros(P,1); 

opt.P = P; 

[FT2p,sol2p] = solve_res_peanut_enhanced(q,U,W,opt);
[FT2,sol2] = solve_res_2B_enhanced(q,U,W,opt);

str = sprintf('Relative residual with peanut compression: %1.2e vs 2B preconditioner without compression: %1.2e\n Converging in %u resp %u iterations', ...
    sol2p.rel_res,sol2.rel_res,sol2p.it,sol2.it);
disp(str)

alignfigs;
