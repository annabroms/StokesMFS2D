function [Q,sol] = solve_cap_peanut(q,v_body,opt)
%SOLVE_CAP_PEANUT Solve exterior Dirichlet Laplace (capacitance:
%known voltages, unknown charges) with peanut-compressed 2B preconditioner.
%
% Syntax:
%   [Q,sol] = solve_cap_peanut(q,v_body,opt)
% Inputs:
%   q          - Complex particle centers (P x 1).
%   v_body     - Constant boundary values per body (P x 1).
%   opt        - Options struct (see getLaplace2Dparams.m).
%     Required fields:
%       rad           physical particle radius
%       N_c,N_f       coarse/fine proxy point counts
%       a_c,a_f       coarse/fine collocation upsampling factors
%       Rp_c,Rp_f     coarse/fine proxy radii
%       delta_pair    pair-detection threshold
%       N_peanut      peanut boundary node count
%       Nclust        total Chebyshev nodes on each enclosing ellipse used
%                     to extract the shielding arc of enhancing sources for
%                     each close pair
%     Solver-control fields:
%       gmres_tol     GMRES tolerance
%       gmres_verbose GMRES print level:
%                     0 = silent, 1 = final summary only,
%                     2 = per-iteration estimated residuals + final summary
%       debug         build/plot/investigate system matrix corresponding to
%                     matvec.
%       visualise_sol plot postprocessing diagnostics
%       use_fmm       use fmm2d (of flatiron) for Laplace field evals
%       cmap          use compressed coarse to coarse map
%       get_bndry_field
%                     if true, reconstruct boundary fields/residuals in
%                     postprocessing
%       body_plot_font_size
%                     base font size used in the voltage/charge subplot figure
%
% Outputs:
%   Q           - Per-body unweighted sums of source strengths.
%   sol         - Struct with fields:
%                 lambda_proxy : compressed coarse source strengths
%                 it           : GMRES iteration count
%                 gmres_tol    : GMRES tolerance used
%                 maxres       : max relative residual on independent boundary points
%                                (NaN if opt.get_bndry_field = 0)
%                 resvec       : GMRES convergence history
%
% Notes:
%   The radius parameter should be chosen with rad ~= 1 to avoid unit logarithmic
%   capacity in 2D.
%
% To test: call without inputs.
%
% See also: solve_cap_1B, solve_cap_2B, solve_elast_peanut, ...
%   transform_lap_peanut, matvec_lap_peanut_enhanced.
%
% Anna Broms, Mar 2026

if nargin==0
    test_solve_cap_peanut;
    return
end

if nargin < 3 || ~isstruct(opt)
    error('solve_cap_peanut requires q, v_body, and an options struct opt.');
end

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-7);
debug = logical(getOptField(opt,'debug',false));
use_fmm = logical(getOptField(opt,'use_fmm',true));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
get_bndry_field = logical(getOptField(opt,'get_bndry_field',true));
body_plot_font_size = getOptField(opt,'body_plot_font_size',14);
opt_solve = opt;
opt_solve.get_bndry_field = false;

q = q(:);
v_body = v_body(:);
P = numel(q);
assert(numel(v_body)==P,'v_body must have one entry per particle.');

%% Parameters
maxit = 800;

if ~exist('solver_name','var') || isempty(solver_name)
    solver_name = mfilename;
end
fprintf('==== START: %s ====\n', solver_name);

R = getOptField(opt,'rad',2);
opt.gmres_verbose = gmres_verbose;

N_c = getOptField(opt,'N_c',80);
N_f = getOptField(opt,'N_f',150);
a_c = getOptField(opt,'a_c',1.2);
a_f = getOptField(opt,'a_f',1.2);

tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',R*max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',R*max([1-sep_f,0.01]));

%% Build grids
nout = ceil(a_c*N_c);
tout = linspace(0,2*pi,nout+1)';
tout = tout(1:end-1);
rbase_out_c = R*(cos(tout)+1i*sin(tout));

tin = linspace(0,2*pi,N_c+1)';
tin = tin(1:end-1);
rbase_in_c = Rp_c*(cos(tin)+1i*sin(tin));

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*(cos(tin_f)+1i*sin(tin_f));

tout_f = linspace(0,2*pi,ceil(a_f*N_f)+1)';
tout_f = tout_f(1:end-1);
rout_base_f = R*(cos(tout_f)+1i*sin(tout_f));

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*nout,1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k)+rbase_in_c;
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

%% Basis factors
[UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_QV,pair_cache] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rimage_vec,refine,pairs,opt);
[UU,YY] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout]);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.opt = opt_solve;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rimage_vec = rimage_vec;
geom.rvec_in = rvec_in_c;
geom.pair_cache = pair_cache;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_QV = Cmap_QV;
basis.pair_cache = pair_cache;
basis.Nii = lapSLPmat(rbase_in_c,rbase_out_c);

%% RHS
fout = zeros(P*nout,1);
for k = 1:P
    fout((k-1)*nout+1:k*nout) = v_body(k);
end

%% Solve
if debug
    x = zeros(length(rout),1);
    CC = zeros(length(rout));
    ncols = length(rout);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_lap_peanut_enhanced(x,geom,basis);
    end
    figure(); imagesc(log10(abs(CC))); colorbar
    title([solver_name ': log_{10}|matvec system matrix|'],'interpreter','none')
    [V,D] = eig(CC);
    D = diag(D);
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title([solver_name ': eigenvalues of matvec system matrix'],'interpreter','none')
end

disp(' == Solving... == ');
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_lap_peanut_enhanced(x,geom,basis), ...
    fout,length(rout),maxit,gmres_tol,opt.gmres_verbose,rout);

figure(); semilogy(resvec)
title('GMRES convergence capacitance peanut','interpreter','latex')

disp(' == Postprocessing == ');
%% Postprocess
if get_bndry_field
    n_bound = 803;
    tb = linspace(0,2*pi,n_bound+1)';
    tb = tb(1:end-1);
    rcheck_b = zeros(P*n_bound,1);
    for k = 1:P
        rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+R*(cos(tb)+1i*sin(tb));
    end

    geom_eval = geom;
    geom_eval.opt = opt;
    geom_eval.rcheck = rcheck_b;
else
    geom_eval = geom;
end

[lam_c,~,~,~,u_corr,pair_qv_nonp,~,lam_self_nonp,lam_f_nonp,lam_e_nonp] = ...
    transform_lap_peanut(tau,geom_eval,basis);

lambda_proxy = lam_c;

if get_bndry_field
    u_b = lapSLPfield(rvec_in_c,rcheck_b,lam_c,use_fmm) + u_corr;

    g_true = zeros(P*n_bound,1);
    for k = 1:P
        g_true((k-1)*n_bound+1:k*n_bound) = v_body(k);
    end

    maxres = max(abs(u_b-g_true))/max(1,max(abs(g_true)));
    fprintf('Max surface relative residual at new nodes %.3e\n',maxres);
else
    u_b = [];
    g_true = [];
    maxres = nan;
    fprintf('Boundary field evaluation skipped (opt.get_bndry_field=0)\n');
end

Q = zeros(P,1);
for k = 1:P
    idx = (k-1)*N_c+1:k*N_c;
    if opt.cmap
        Q(k) = sum(lam_self_nonp(idx)) + pair_qv_nonp(k);
    else
        Q(k) = sum(lam_self_nonp(idx)) + sum(lam_f_nonp{k}) + sum(lam_e_nonp{k});
    end
end

if visualise_sol
    if get_bndry_field
        figure();
        plot(u_b); hold on; plot(g_true)
        title('Boundary values: lhs vs rhs (capacitance peanut)', ...
            'Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
    end

    figure();
    semilogy(abs(lambda_proxy))
    title('Compressed source strengths (capacitance peanut)', ...
        'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')
    axis tight

    plotBodyScalars(q,R,v_body,Q,body_plot_font_size);
end

sol = struct();
sol.lambda_proxy = lambda_proxy;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.maxres = maxres;
sol.resvec = resvec;

end

function test_solve_cap_peanut
fprintf('--- solve_cap_peanut self-test ---\n');

close all; 
rng(8);

% Set geometry and data
R = 2;
P = 10;
delta = 1e-3;
q = grow_cluster(P,delta,2,R);

% Solve capacitance for hexagonal lattice
%rings = 8; 
%rings = 2; 
%q = hexagonal_lattice(delta,rings,R);
P = length(q); 
v_body = buildAlternatingVoltages(q,R);
check_multi_compress = 0; 

% Set parameters and settings
N_c = 60; 
opt = getLaplace2Dparams(P,R,N_c);
opt.delta_pair = 0.2;
opt.Nclust = 100;
opt.N_peanut = 400;
opt.visualise_sol = 1;
opt.visualise_grid =0; 
opt.gmres_tol = 1e-8;
opt.debug = 0;
opt.use_fmm = true;
opt.gmres_verbose = 0;
opt.compress_cmap = 0; %use low rank approximation of coarse-coarse map
opt.cmap_tol = 1e-8; % tolerance used in the low-rank compression
opt.reuse_pair_basis_by_sep = 0; 
%tic
[Qp,solp] = solve_cap_peanut(q,v_body,opt);
%t_one = toc;

opt.visualise_sol = 0;
opt.reuse_pair_basis_by_sep = 0; 
if check_multi_compress
  %  tic 
    [Qp2,solp2] = solve_cap_peanut(q,v_body,opt);
   % t_recomp = toc;
    fprintf('Peanut solution times: one pair=%1.2d s, multipair=%1.2d s, iters: %u vs %u, soldiff: %1.2d\n',t_one,...
        t_recomp,solp.it,solp2.it,norm(norm(Qp-Qp2)/max(1,norm(Qp2))));
else
%    fprintf('Solved and postprocessed in %1.2d s\n',t_one);
end
[Q2,sol2] = solve_cap_2B(q,v_body,opt);
itp = solp.it;
resp = solp.maxres;
it2 = sol2.it;
res2 = sol2.maxres;

fprintf('Peanut: it=%d, maxres=%.3e\n',itp,resp);
fprintf('2B    : it=%d, maxres=%.3e\n',it2,res2);
fprintf('Rel diff in Q (peanut vs 2B): %.3e\n',norm(Qp-Q2)/max(1,norm(Q2)));

end

function plotBodyScalars(q,R,v_body,Q,font_size)
q = q(:);
v_body = v_body(:);
Q = Q(:);
P = numel(q);

if nargin < 5 || isempty(font_size)
    font_size = 14;
end

if isscalar(R)
    rad = repmat(R,P,1);
else
    rad = R(:);
end

theta = linspace(0,2*pi,200);
xmin = min(real(q)-rad);
xmax = max(real(q)+rad);
ymin = min(imag(q)-rad);
ymax = max(imag(q)+rad);
pad = 0.1*max(rad);

figure();
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

vals = {v_body,Q};
titles = {'Given body voltages','Computed net charges'};
cbar_labels = {'Voltage','Net charge'};
cmaps = {parula(256),blueWhiteRedMap(256)};

for it = 1:2
    ax = nexttile;
    hold(ax,'on');

    for k = 1:P
        zk = q(k) + rad(k)*(cos(theta)+1i*sin(theta));
        fill(ax,real(zk),imag(zk),vals{it}(k), ...
            'EdgeColor',[0.2 0.2 0.2],'LineWidth',0.75);
    end

    axis(ax,'equal');
    xlim(ax,[xmin-pad xmax+pad]);
    ylim(ax,[ymin-pad ymax+pad]);
    xlabel(ax,'$x$','Interpreter','latex');
    ylabel(ax,'$y$','Interpreter','latex');
    title(ax,titles{it},'Interpreter','latex');
    box(ax,'on');
    colormap(ax,cmaps{it});
    ax.FontSize = font_size;
    ax.TitleFontSizeMultiplier = 1.0;
    ax.LabelFontSizeMultiplier = 1.0;
    ax.TickLabelInterpreter = 'latex';

    if min(vals{it}) < 0 && max(vals{it}) > 0
        vmax = max(abs(vals{it}));
        clim(ax,[-vmax vmax]);
    end

    c = colorbar(ax);
    c.Label.String = cbar_labels{it};
    c.FontSize = font_size;
    c.Label.FontSize = font_size;
    c.TickLabelInterpreter = 'latex';
    c.Label.Interpreter = 'latex';
end

sgtitle('Bodywise voltages and net charges', ...
    'FontSize',font_size,'Interpreter','latex');
end

function cmap = blueWhiteRedMap(n)
if nargin < 1 || isempty(n)
    n = 256;
end

x = linspace(0,1,n)';
cmap = zeros(n,3);

mid = 0.5;
left = x <= mid;
right = x > mid;

tleft = x(left)/mid;
tright = (x(right)-mid)/(1-mid);

cmap(left,1) = tleft;
cmap(left,2) = tleft;
cmap(left,3) = 0.3 + 0.7*tleft;

cmap(right,1) = 1.0;
cmap(right,2) = 1.0 - 0.8*tright;
cmap(right,3) = 1.0 - tright;
end
