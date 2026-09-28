% apr26_reuse_sep_rotation_compare.m
%
% Test whether transform_lap_peanut gives matching results for
% reuse_pair_basis_by_sep=0 (per-pair basis, no FFT rotation) vs =1
% (canonical group basis, FFT rotation) across a sweep of inter-row shift
% angles.
%
% Geometry: a rectangular multi-row packing where all nearest-neighbour
% pairs have the same surface gap delta, but the offset between consecutive
% rows is parametrised by alpha in [0,1]:
%
%   Row ir: q = ic*d + ir * (alpha*d + 1i*sqrt(1-alpha^2)*d)
%
% With alpha = 0.5 this is the standard equilateral-triangle / hexagonal
% packing (all inter-row pairs at angle 60 degrees). For other alpha the
% inter-row angle differs and no equilateral triangles are formed.  All
% nearest-neighbour pairs still share the same centre-to-centre separation
% d, so reuse_pair_basis_by_sep=1 will group them all into one canonical
% basis—relying on the FFT rotation to map between the canonical frame
% (real axis) and the actual pair angle.
%
% For each alpha the script:
%   1. Builds geom/basis for reuse=0 and reuse=1.
%   2. Calls transform_lap_peanut with both and a random tau.
%   3. Reports relative errors in lam_c, lam_c_nonp, u_corr, pair_qv_nonp.
%
% If the FFT rotation is exact the errors should be at machine-precision
% level for all alpha. Any systematic growth with alpha reveals a rotation
% inaccuracy and its nature (e.g., 'oversampled_fft' vs 'fft' modes).
%
% Anna Broms, Apr 26, 2026

script_name = mfilename;
script_date = 'Apr 26, 2026';
repo_root   = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (%s) ===\n\n', script_name, script_date);

% -------------------------------------------------------------------------
%  Parameters (can be overridden before running the script)
% -------------------------------------------------------------------------
if ~exist('n_cols','var')     || isempty(n_cols),     n_cols     = 4;    end
if ~exist('n_rows','var')     || isempty(n_rows),     n_rows     = 3;    end
if ~exist('R','var')          || isempty(R),          R          = 2;    end
if ~exist('delta','var')      || isempty(delta),      delta      = 1e-3; end  % relative surface gap
if ~exist('N_c','var')        || isempty(N_c),        N_c        = 40;   end
if ~exist('N_f','var')        || isempty(N_f),        N_f        = 80;   end
if ~exist('rng_seed','var')   || isempty(rng_seed),   rng_seed   = 260426; end
% Alpha controls the horizontal shift of each next row as a fraction of d.
% alpha = 0.5  ->  equilateral triangle / standard hex packing (angle 60 deg)
% alpha != 0.5 ->  non-equilateral, inter-row angle varies
if ~exist('alpha_vals','var') || isempty(alpha_vals)
    alpha_vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
end
% Rotation mode for the FFT-based pair rotation
% NOTE: In the current implementation, 'oversampled_fft' zero-pads the
% spectrum of the same N samples before shifting and then samples back on
% the N-grid. It therefore agrees with 'fft' to roundoff and does not add
% information or resolve modes that were unresolved on the original grid.
if ~exist('rot_modes','var') || isempty(rot_modes)
    rot_modes = {'fft', 'oversampled_fft'};
end

rng(rng_seed, 'twister');
d = 2*R*(1 + delta);   % nearest-neighbour centre-to-centre separation

% -------------------------------------------------------------------------
%  Base options (shared by both reuse variants)
% -------------------------------------------------------------------------
P_max = n_rows * n_cols;
opt_base = getLaplace2Dparams(P_max, R, N_c, N_f);
% This threshold includes all intended nearest neighbours, but it does not
% select only those neighbours for every alpha. In particular, for the
% default geometry alpha=0.4 and 0.6 also have diagonal pairs with surface
% gap about 0.386 when R=2, below 0.2*R=0.4. Those cases consequently form
% two separation groups and do not isolate a single shared separation. Use
% a smaller value such as 0.05*R when that isolation is desired.
opt_base.delta_pair    = 0.2 * R;
opt_base.N_peanut      = 100;
opt_base.cmap          = 1;
opt_base.compress_cmap = 0;
opt_base.get_bndry_field = 0;
opt_base.use_fmm       = false;
opt_base.show_counter  = 0;
opt_base.visualise_sol = 0;
opt_base.visualise_grid = 0;
opt_base.gmres_verbose = 0;
opt_base.single_threaded = 1;

% -------------------------------------------------------------------------
%  Run sweep
% -------------------------------------------------------------------------
for im = 1:numel(rot_modes)
    rot_mode = rot_modes{im};
    fprintf('Rotation mode: %s\n', rot_mode);
    fprintf('%-6s  %-9s  %-8s  %-8s  %-12s  %-12s  %-12s  %-12s\n', ...
        'alpha', 'angle_deg', 'P', 'n_pairs', ...
        'err_lam_c', 'err_lam_cnp', 'err_u_corr', 'err_pqv');
    fprintf('%s\n', repmat('-', 1, 95));

    for ia = 1:numel(alpha_vals)
        alpha = alpha_vals(ia);

        % -----------------------------------------------------------------
        %  Build geometry: n_rows x n_cols discs, all NN pairs at dist d
        %
        %  Row offset: row i starts at ic*d + ir*(row_dx + 1i*row_dy).
        %  For alpha < 0.5 the close inter-row neighbour is the same-column
        %  particle (at distance d).  For alpha > 0.5 it is the column-1
        %  neighbour, also at distance d.  Both cases require
        %    row_dy = d * sqrt(1 - min(alpha, 1-alpha)^2)
        %  so that the nearest inter-row distance equals d exactly and no
        %  discs overlap.
        % -----------------------------------------------------------------
        row_dx = alpha * d;
        row_dy = d * sqrt(max(1 - min(alpha, 1-alpha)^2, 0));
        q = zeros(n_rows * n_cols, 1);
        for ir = 0:n_rows-1
            for ic = 0:n_cols-1
                q(ir*n_cols + ic + 1) = ic*d + ir*(row_dx + 1i*row_dy);
            end
        end
        P = numel(q);
        inter_row_angle_deg = atan2d(row_dy, row_dx);

        opt = opt_base;
        opt.P = P;
        opt.rotation_mode = rot_mode;

        % -----------------------------------------------------------------
        %  Build grids
        % -----------------------------------------------------------------
        nout = ceil(opt.a_c * N_c);
        t_in  = linspace(0,2*pi,N_c+1)';     t_in  = t_in(1:end-1);
        t_out = linspace(0,2*pi,nout+1)';    t_out = t_out(1:end-1);
        t_f   = linspace(0,2*pi,N_f+1)';     t_f   = t_f(1:end-1);
        t_of  = linspace(0,2*pi,ceil(opt.a_f*N_f)+1)'; t_of = t_of(1:end-1);

        rbase_in_c  = opt.Rp_c * (cos(t_in)  + 1i*sin(t_in));
        rbase_out_c = R        * (cos(t_out) + 1i*sin(t_out));
        rbase_in_f  = opt.Rp_f * (cos(t_f)   + 1i*sin(t_f));
        rout_base_f = R        * (cos(t_of)  + 1i*sin(t_of));

        rvec_in_c = zeros(P*N_c, 1);
        rout      = zeros(P*nout, 1);
        for k = 1:P
            rvec_in_c((k-1)*N_c+1:k*N_c) = q(k) + rbase_in_c;
            rout((k-1)*nout+1:k*nout)     = q(k) + rbase_out_c;
        end

        [~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q, opt);
        [U, Y] = getSelfPseudoLaplace(1, rbase_in_c, rbase_out_c, [0 nout]);

        % -----------------------------------------------------------------
        %  Variant 0: reuse_pair_basis_by_sep = false (no rotation)
        % -----------------------------------------------------------------
        opt0 = opt;
        opt0.reuse_pair_basis_by_sep = false;
        [Uf0,Yf0,Up0,Yp0,Cm0,CmQV0,pc0] = getPairBasisLaplace( ...
            q, rbase_in_c, rbase_in_f, rout_base_f, rbase_out_c, ...
            rimage_vec, refine, pairs, opt0);

        geom0 = makeGeom(q, rbase_in_c, rbase_in_f, refine, opt0, ...
            rout, pairs, rimage_vec, rvec_in_c, pc0);
        basis0 = makeBasis(U, Y, Uf0, Yf0, Up0, Yp0, Cm0, CmQV0, pc0, ...
            rbase_in_c, rbase_out_c);

        % -----------------------------------------------------------------
        %  Variant 1: reuse_pair_basis_by_sep = true (FFT rotation)
        % -----------------------------------------------------------------
        opt1 = opt;
        opt1.reuse_pair_basis_by_sep = true;
        [Uf1,Yf1,Up1,Yp1,Cm1,CmQV1,pc1] = getPairBasisLaplace( ...
            q, rbase_in_c, rbase_in_f, rout_base_f, rbase_out_c, ...
            rimage_vec, refine, pairs, opt1);

        geom1 = makeGeom(q, rbase_in_c, rbase_in_f, refine, opt1, ...
            rout, pairs, rimage_vec, rvec_in_c, pc1);
        basis1 = makeBasis(U, Y, Uf1, Yf1, Up1, Yp1, Cm1, CmQV1, pc1, ...
            rbase_in_c, rbase_out_c);

        % -----------------------------------------------------------------
        %  Compare transform_lap_peanut outputs
        % -----------------------------------------------------------------
        % randn is deliberately an operator stress test, not smooth
        % boundary data and not a convergence test in N_c. After applying
        % the one-body pseudoinverse, these samples can produce coarse proxy
        % coefficients dominated by the highest Fourier modes (including
        % the ambiguous even-N Nyquist mode). Increasing N_c while drawing
        % a new length-dependent random vector continues to populate the
        % new Nyquist range. For a resolution/convergence study, instead
        % sample one fixed smooth function, or a random Fourier series with
        % a fixed maximum mode, on every boundary grid. See also
        % sep23_laplace_offgrid_pair_rotation.m.
        rng(rng_seed, 'twister');
        tau = randn(P * nout, 1);

        [lam_c0, ~, ~, ~, u0, pqv0, lam_cn0] = transform_lap_peanut(tau, geom0, basis0);
        [lam_c1, ~, ~, ~, u1, pqv1, lam_cn1] = transform_lap_peanut(tau, geom1, basis1);

        ref_lam  = max(1, norm(lam_c0));
        ref_lcnp = max(1, norm(lam_cn0));
        % u_corr is a difference of the fine and compressed fields. Its
        % norm can be small through cancellation, so this relative metric
        % can look large even when both component fields agree well. A
        % component-scaled field metric is used in the Sep 23 diagnostic.
        ref_u    = max(1, norm(u0));
        ref_pqv  = max(1, norm(pqv0));

        err_lam  = norm(lam_c1  - lam_c0)  / ref_lam;
        err_lcnp = norm(lam_cn1 - lam_cn0) / ref_lcnp;
        err_u    = norm(u1      - u0)       / ref_u;
        err_pqv  = norm(pqv1    - pqv0)     / ref_pqv;

        n_groups = pc1.n_groups;
        fprintf('%-6.2f  %-9.2f  %-8d  %-8d  %-12.3e  %-12.3e  %-12.3e  %-12.3e   [%d grp]\n', ...
            alpha, inter_row_angle_deg, P, size(pairs,1), ...
            err_lam, err_lcnp, err_u, err_pqv, n_groups);
    end
    fprintf('\n');
end

% -------------------------------------------------------------------------
%  Helper functions
% -------------------------------------------------------------------------

function geom = makeGeom(q, rbase_in_c, rbase_in_f, refine, opt, ...
        rout, pairs, rimage_vec, rvec_in_c, pair_cache)
geom = struct();
geom.rbase_in_c  = rbase_in_c;
geom.rbase_in_f  = rbase_in_f;
geom.refine      = refine;
geom.opt         = opt;
geom.rvec_out    = rout;
geom.rcheck      = rout;
geom.q           = q;
geom.pairs       = pairs;
geom.rimage_vec  = rimage_vec;
geom.rvec_in     = rvec_in_c;
geom.pair_cache  = pair_cache;
end

function basis = makeBasis(U, Y, Upf, Ypf, DC_all, YC_all, Cmap, Cmap_QV, ...
        pair_cache, rbase_in_c, rbase_out_c)
basis = struct();
basis.U         = U;
basis.Y         = Y;
basis.Upf       = Upf;
basis.Ypf       = Ypf;
basis.DC_all    = DC_all;
basis.YC_all    = YC_all;
basis.Cmap      = Cmap;
basis.Cmap_QV   = Cmap_QV;
basis.pair_cache = pair_cache;
basis.Nii       = lapSLPmat(rbase_in_c, rbase_out_c);
end
