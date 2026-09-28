clear;
close all;
clc;
set(0,'DefaultFigureVisible','off');

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Laplace Pair-Rotation Comparison (Mar 20, 2026) ===\n');

rng(20);

R = 2;
gap = 0.001;
theta = pi/5;
q = [0; (2*R + gap)*exp(1i*theta)];
P = numel(q);
N_c = 60;
N_f = 150;

opt = getLaplace2Dparams(P,R,N_c,N_f);
opt.use_fmm = false;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.debug = 0;
opt.delta_pair = 0.2*R;
opt.N_peanut = 200;
opt.compress_cmap = 0;
opt.reuse_pair_basis_by_sep = true;
opt.check_rotations = true;
opt.project_charge = true;

[rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] = build_circle_grids( ...
    opt.rad,opt.N_c,opt.N_f,opt.a_c,opt.a_f,opt.Rp_c,opt.Rp_f);

rvec_out = zeros(P*nout,1);
for k = 1:P
    rvec_out((k-1)*nout+1:k*nout) = q(k) + rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
assert(size(pairs,1)==1,'Expected exactly one close pair for this test.');

[~,~,~,~,~,~,pair_cache] = getPairBasisLaplace( ...
    q,rbase_in_c,rbase_in_f,rout_base_f,rbase_out_c, ...
    rimage_vec,refine,pairs,opt);
[Uii,Yii] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout],true);

tau = randn(P*nout,1);
lam_c = zeros(P*opt.N_c,1);
lam_c_nonp = zeros(P*opt.N_c,1);
for k = 1:P
    tau_k = tau((k-1)*nout+1:k*nout);
    lam_k_nonp = Yii{1}*(Uii{1}*tau_k);
    lam_k = project_charge_mode(lam_k_nonp,opt.project_charge);

    idx = (k-1)*opt.N_c+1:k*opt.N_c;
    lam_c(idx) = lam_k;
    lam_c_nonp(idx) = lam_k_nonp;
end

lam_i = lam_c(1:opt.N_c);
lam_p2 = lam_c(opt.N_c+1:2*opt.N_c);

meta = pair_cache.meta(1);
canon = pair_cache.groups(1);
exact = pair_cache.check_pairs(1);

fprintf('Geometry: R=%.3f, gap=%.3f, theta=%.3f rad\n',R,gap,theta);
fprintf('Pair cache groups: %d\n',pair_cache.n_groups);

for cmap_mode = [1 0]
    opt.cmap = cmap_mode;

    [beta_canon,tau_canon,rot_stats] = compare_pair_outputs( ...
        canon,meta,lam_i,lam_p2,opt,true,cmap_mode);
    [beta_exact,tau_exact] = compare_pair_outputs(exact,meta,lam_i,lam_p2,opt,false,cmap_mode);

    fprintf('\n  cmap=%d\n',cmap_mode);
    fprintf('    rhs sum relerr (pre/post FFT rotation)  = %.3e\n', ...
        relerr(rot_stats.rhs_pre_sum,rot_stats.rhs_post_sum));
    fprintf('    beta sum relerr (pre/post FFT rotation) = %.3e\n', ...
        relerr(rot_stats.beta_pre_sum,rot_stats.beta_post_sum));
    fprintf('    tau sum relerr (pre/post FFT rotation)  = %.3e\n', ...
        relerr(rot_stats.tau_pre_sum,rot_stats.tau_post_sum));
    fprintf('    beta projected relerr     = %.3e\n',relerr(beta_canon(:,1:2),beta_exact(:,1:2)));
    fprintf('    beta non-projected relerr = %.3e\n',relerr(beta_canon(:,3:4),beta_exact(:,3:4)));
    fprintf('    tau projected relerr      = %.3e\n',relerr(tau_canon(:,1:2),tau_exact(:,1:2)));
    fprintf('    tau non-projected relerr  = %.3e\n',relerr(tau_canon(:,3:4),tau_exact(:,3:4)));
end

fprintf('\nDone.\n');

function [rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] = build_circle_grids( ...
    rad,N_c,N_f,a_c,a_f,Rp_c,Rp_f)
nout = ceil(a_c*N_c);

t = linspace(0,2*pi,N_c+1)';
t = t(1:end-1);
rbase_in_c = Rp_c*(cos(t)+1i*sin(t));

t = linspace(0,2*pi,nout+1)';
t = t(1:end-1);
rbase_out_c = rad*(cos(t)+1i*sin(t));

t = linspace(0,2*pi,N_f+1)';
t = t(1:end-1);
rbase_in_f = Rp_f*(cos(t)+1i*sin(t));

% Fine collocation points are not needed explicitly in this comparison, but
% the rout_base_f grid is required when we build the pair basis.
t = linspace(0,2*pi,ceil(a_f*N_f)+1)';
t = t(1:end-1);
rout_base_f = rad*(cos(t)+1i*sin(t));
end

function [beta_bundle,tau_bundle,rot_stats] = compare_pair_outputs(pair,meta,lam_i,lam_p2,opt,use_rotation,use_cmap)
N_f = opt.N_f;
N_c = opt.N_c;
project_charge = logical(opt.project_charge);
rot_stats = struct('rhs_pre_sum',[],'rhs_post_sum',[], ...
    'beta_pre_sum',[],'beta_post_sum',[], ...
    'tau_pre_sum',[],'tau_post_sum',[]);

if use_rotation
    rhs_pre = [lam_i lam_p2];
    rot_stats.rhs_pre_sum = sum(rhs_pre,1);
    rhs_pair = rotateUniformCircleData(rhs_pre,[],meta.phase_c);
    rot_stats.rhs_post_sum = sum(rhs_pair,1);
    rhs_pair = rhs_pair(:);
else
    rhs_pair = [lam_i; lam_p2];
end

pair_mapped = pair.Upf*rhs_pair;
beta_tot_nonp_local = pair.Ypf*pair_mapped;

im_i = numel(pair.rimage_canon{1});
im_p2 = numel(pair.rimage_canon{2});

s_i = 1:N_f;
e_i = N_f+1:N_f+im_i;
s_p2 = N_f+im_i+1:2*N_f+im_i;
e_p2 = 2*N_f+im_i+1:2*N_f+im_i+im_p2;

beta_i_nonp_local = [beta_tot_nonp_local(s_i); beta_tot_nonp_local(e_i)];
beta_p2_nonp_local = [beta_tot_nonp_local(s_p2); beta_tot_nonp_local(e_p2)];
beta_i_local = project_charge_mode(beta_i_nonp_local,project_charge);
beta_p2_local = project_charge_mode(beta_p2_nonp_local,project_charge);

beta_bundle = [beta_i_local(1:N_f) beta_p2_local(1:N_f) ...
               beta_i_nonp_local(1:N_f) beta_p2_nonp_local(1:N_f)];
rot_stats.beta_pre_sum = sum(beta_bundle,1);
if use_rotation
    beta_bundle = rotateUniformCircleData(beta_bundle,[],meta.phase_f_inv);
end
rot_stats.beta_post_sum = sum(beta_bundle,1);

if use_cmap
    tau_peanut_nonp_local = pair.Cmap*rhs_pair;
else
    tau_peanut_nonp_local = pair.YC*(pair.DC*beta_tot_nonp_local);
end

tau_peanut_nonp_pair = reshape(tau_peanut_nonp_local,N_c,2);
tau_peanut_pair = [project_charge_mode(tau_peanut_nonp_pair(:,1),project_charge) ...
    project_charge_mode(tau_peanut_nonp_pair(:,2),project_charge)];

tau_bundle = [tau_peanut_pair tau_peanut_nonp_pair];
rot_stats.tau_pre_sum = sum(tau_bundle,1);
if use_rotation
    tau_bundle = rotateUniformCircleData(tau_bundle,[],meta.phase_c_inv);
end
rot_stats.tau_post_sum = sum(tau_bundle,1);
end

function lam_out = project_charge_mode(lam_in,do_project)
if do_project
    lam_out = lam_in - mean(lam_in);
else
    lam_out = lam_in;
end
end

function e = relerr(A,B)
e = norm(A-B,'fro') / max(1,norm(B,'fro'));
end
