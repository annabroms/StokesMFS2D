function opt = get2Dparams(P,N_c,N_f)
%get2Dparams returns struct with default parameters for a Stokes MFS 2D problem

if nargin<2
    N_c = 60; % sources on coarse proxy grid
    N_c = 100;
    N_c = 150;
    N_f = 150; % sources on fine proxy grid
elseif nargin<3
    N_f = 150;
end

opt.P = P; %number of particles
opt.N_c = N_c; 

%set separation between true boundary and proxy surface
tol = 1e-12;
sep = (1/N_c)*log(1/tol); %separation between surfaces based on rule of thumb in Stein & Barnett QFS paper (2021) 
opt.Rp_c = max([1-sep,0.01]); %radius of proxy surface

opt.a_c = 1.2; % upsampling factor for collocation points so that M_c = a_c*N_c; 
opt.image = 1;%use images

% Image based enhancement
opt.M_image = 35;  % Not always in use - sets extra collocation points for close to touching region
opt.s = [0 0 1 1 0 0 0]; %source types at clustered nodes
opt.alpha = 0; %use two lines?
opt.proj_all = 0; %only for the mobility problem

% Enhancement with Stokeslets only, using ellipse segments
opt.beta = 0.3;
opt.Nclust = 100; 

% Pair corrections
opt.pc = 0; %use pair corrections?
opt.delta_pair = 0.2; %activate pair_corrections
opt.N_f = N_f;
tol = 1e-12;
sep = (1/N_f)*log(1/tol); %separation between surfaces based on rule of thumb in Stein & Barnett QFS paper (2021) 
opt.Rp_f = max([1-sep,0.01]);
opt.a_f = 1.2; % upsampling factor for collocation points so that M_f = a_f*N_f; 
 
opt.beta = 0.3; %determines distance between focus and tip of ellipse for enhancing nodes
opt.cmap = 1; %use coarse-to-coarse mapping. Only applicable with peanut compression.
opt.precomp = 1; % Store evaluation of onebody field from precomputation
opt.gmres_verbose = 0; % 0=silent, 1=final summary, 2=per-iteration
opt.visualise_sol = 0; % plot postprocessing diagnostics
opt.visualise_grid = 0; % plot source/collocation grids during setup
opt.get_bndry_field = 1; % determine velocity on the boundary in post-processing?
opt.project_force = false; % project out net force/torque-producing modes in Stokes mobility solves
opt.reuse_pair_basis_by_sep = true; % reuse canonical pair bases for repeated separations
opt.parallel_precomp = false; % parallelise pair-basis builds when a parallel pool is available
opt.shared_sep_tol = 1e-4; % separation tolerance for grouping repeated pairs
opt.rotation_mode = 'oversampled_fft'; % 'fft' | 'oversampled_fft' for cached pair rotations
opt.rotation_oversample = 8; % oversampling factor used when rotation_mode = 'oversampled_fft'
opt.use_matrix_free_Lc_pair = true; % matrix-free pair rigid projection; set false for dense comparison
opt.show_counter = 1; % show progress for pair compressions
opt.self_correct = 1; % enforce identity diagonal matrices in system matrix
opt.use_dense = 1; % use stored matrices for evaluation of Stokeslet on single body / pair
opt.pair_basis_debug = 0; % check accuracy in fine and peanut least squares problem
opt.column_weight = false; % scale LS operator columns before SVD in peanut solvers
opt.left_weight = false; % scale LS operator rows by local arclength weights in peanut solvers

end
