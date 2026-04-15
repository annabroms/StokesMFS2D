function opt = getLaplace2Dparams(P,R,N_c,N_f)
%GETLAPLACE2DPARAMS Default parameters for scalar Laplace MFS in 2D.
%
% Input: P - number of particles.
%        R - physical particle radius
%        N_c - Number of coarse proxy points per particle
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
opt.a_c = 1.2;
tol = 1e-12;
sep = (1/N_c)*log(1/tol);
opt.Rp_c = opt.rad*max([1-sep,0.01]);

% Enhancing discretisation for close interations
%opt.delta_pair = (opt.rad-opt.Rp_c)^2/opt.Rp_c;
opt.delta_pair = 0.2*R; % largest distance for which pair corrections are applied. 
opt.beta = 0.3; % beta is a parameter determining the shape of the enhancing ellipse 
% segments for close pairs. Smaller beta means tip of ellipse closer to image accumulation points.
opt.Nclust = 100; % Chebyshev nodes on each ellipse segment for close pairs, a portion of which are used as enhancing sources.
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
opt.reuse_pair_basis_by_sep = true; % build one canonical x-axis pair basis per repeated separation
opt.parallel_precomp = false; % parallelise pair-basis builds when a parallel pool is available
opt.check_rotations = false; % store per-pair pair-basis data alongside the canonical cache for debugging
opt.shared_sep_tol = 1e-2*max(1,opt.rad); % separation matching tolerance used when grouping close pairs
opt.rotation_mode = 'oversampled_fft'; % 'fft' | 'oversampled_fft' for cached pair rotations
opt.rotation_oversample = 8; % oversampling factor used when rotation_mode = 'oversampled_fft'

% Solver and postprocessing control fields:
opt.use_fmm = true;
opt.debug = 0; % build/plot/investigate system matrix corresponding to matvec
opt.gmres_tol = 1e-7;
opt.project_charge = false; % false for capacitance problems, true for elastance problems
opt.gmres_verbose = 0; % 0=silent, 1=final summary, 2=per-iteration
opt.visualise_sol = 0; % draw solution/postprocessing quantities after solve
opt.visualise_grid = 0; % draw source and collocation points at setup stage
opt.get_bndry_field = 1; % reconstruct/evaluate boundary fields in postprocessing?

end
