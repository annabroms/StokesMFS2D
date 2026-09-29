function opt = getLaplace2Dparams(P,R,N_c,N_f)
%GETLAPLACE2DPARAMS Default parameters for scalar Laplace MFS in 2D.
%
% Input: P - number of particles.
%        R - physical particle radius
%        N_c - Number of coarse proxy points per particle used by the
%              global solve and interaction representation
%        N_f - Number of fine proxy points per particle
%
% See also: solve_cap_1B, solve_cap_2B, solve_cap_peanut.
%
% Anna Broms, Mar 2026

% Discretisation parameters, basic discretisation:
opt.rad = R;
opt.P = P; 
if nargin<4
    N_f = 150;
    N_c = 80;
elseif nargin<3
    N_c = 80;
end
opt.N_c = N_c;
% Number of coarse sources in the canonical Cmap representation. Keeping
% this equal to N_c recovers the original square-map implementation. The
% capacitance peanut solver permits N_cmap ~= N_c only with canonical
% pair-basis reuse and charge-preserving field refits around Cmap.
opt.N_cmap = N_c;
opt.refit = false; % false: equal-grid Fourier rotation; true: per-particle MFS field refit
opt.a_c = 1.2;
tol = 1e-12;
sep = (1/N_c)*log(1/tol);
opt.Rp_c = opt.rad*max([1-sep,0.01]);

% Enhancing discretisation for close interations
%opt.delta_pair = (opt.rad-opt.Rp_c)^2/opt.Rp_c;
opt.delta_pair = 0.2*R; % largest distance for which pair corrections are applied. 
opt.beta = 0.3; % beta is a parameter determining the shape of the enhancing ellipse 
% segments for close pairs. Smaller beta means tip of ellipse closer to image accumulation points.
opt.Nclust = 150; % Chebyshev nodes on each ellipse segment for close pairs, a portion of which are used as enhancing sources.
opt.ellipse_constant = false; % if true, freeze the ellipse-segment discretisation at opt.smallest_delta
opt.smallest_delta = 1e-3*R; % smallest gap ever requested; required when opt.ellipse_constant=true
opt.a_f = 1.2;
opt.N_f = N_f;
sep = (1/N_f)*log(1/tol);
opt.Rp_f = opt.rad*max([1-sep,0.01]);

% Only relevant with pair corrections (2B/peanut)
opt.N_peanut = 400; %Number of nodes on peanut separation surface
opt.cmap = 1; % use compressed coarse to coarse map for pair compression? Only relevant with peanut compression.
opt.show_counter = 1; % show progress for pre-computation step for all pairs
opt.pc = 1; %Do pair corrections? %% Is this field still active?
opt.compress_cmap  = 0; % low rank compression of cmap
opt.cmap_tol = 1e-8; %tolerance in the low rank compression
opt.use_tikhonov = false; % smoothly regularize pair/peanut pseudoinverses in two-body setup
% Relative parameter lambda/sigma_max for those Tikhonov filters. Empty
% uses the legacy 1e-14 pair and peanut pseudoinverse tolerance;
% interpolation instead selects its 1e-11 default.
opt.tikhonov_tol = [];
% Pair-map construction for Laplace capacitance:
%   'none'            exact map for every distinct separation (default off)
%   'full'            interpolate the full canonical Cmap in alpha
%   'reduced'         interpolate Cref+U*B(alpha)*V' panelwise
%   'reduced_noconst' interpolate U*B(alpha)*V' panelwise (preferred mode)
% The interpolation modes use canonical rotations for both refit=false
% (Fourier rotation) and refit=true (charge-preserving field refits).
% Elastance supports only 'full'; capacitance supports all listed modes.
% Interpolation remains opt-in and requires a compatible saved model.
opt.use_interpolation = 'none';
% Accuracy targets for the canonical alpha interpolator.  The first is
% applied to the final full/reconstructed coarse correction map C; the
% second is applied independently to the 2-by-(2*N_cmap) voltage/charge map.
opt.interpolation_tol = 1e-6;
opt.volt_charge_interp_tol = 1e-8;
% Empty selects a deterministic, parameter-keyed file under data/.
% Pass 'capacitance' or 'elastance' to prepareLaplaceCmapInterpolation.
opt.interpolation_model_file = '';
% Adaptive training search.  Defaults reproduce the documented current
% search scale; most runs only need to change the two tolerances above.
opt.interpolation_panel_count_candidates = [1 2 4 8];
opt.interpolation_node_candidates = 3:2:17;
opt.interpolation_validation_nodes = 33;
opt.reuse_pair_basis_by_sep = true; % build one canonical x-axis pair basis per repeated separation
opt.parallel_precomp = false; % parallelise pair-basis builds when a parallel pool is available
opt.check_rotations = false; % store per-pair pair-basis data alongside the canonical cache for debugging
opt.shared_sep_tol = 1e-6*max(1,opt.rad); % separation matching tolerance used when grouping close pairs
opt.rotation_mode = 'oversampled_fft'; % 'fft' | 'oversampled_fft' for cached pair rotations
opt.rotation_oversample = 8; % oversampling factor used when rotation_mode = 'oversampled_fft'
opt.use_big_sparse = false; % use global sparse close-pair correction matrices in peanut GMRES
opt.lap_big_sparse_build_mode = 'auto'; % 'auto' | 'precomputed' | 'streaming'
opt.lap_big_sparse_chunk_pairs = 8; % pairs per sparse-triplet assembly chunk
opt.lap_big_sparse_max_build_bytes = inf; % guardrail for estimated sparse-build peak memory
opt.lap_big_sparse_max_build_ram_fraction = 1; % optional fraction of current MemAvailable for sparse-build peak guardrail

% Solver and postprocessing control fields:
opt.use_fmm = true;
opt.debug = 0; % build/plot/investigate system matrix corresponding to matvec
opt.gmres_tol = 1e-7;
opt.project_charge = false; % false for capacitance problems, true for elastance problems
opt.gmres_verbose = 0; % 0=silent, 1=final summary, 2=per-iteration
opt.visualise_sol = 0; % draw solution/postprocessing quantities after solve
opt.visualise_grid = 0; % draw source and collocation points at setup stage
opt.get_bndry_field = 1; % reconstruct/evaluate boundary fields in postprocessing?
opt.get_solve_time = true; % measure GMRES wall-clock solve time and the FMM subset
opt.RAM_check = false; % estimate/report RAM usage for precomp, solve, and postprocessing
opt.single_threaded = 0;

end
