function beta_rot = rotateStokesPairSourceVector(beta,Nf,im_i,im_j,phase_f,vec_rot)
%ROTATESTOKESPAIRSOURCEVECTOR Rotate pair-ordered fine/image Stokes sources.
%
% Syntax:
%   beta_rot = rotateStokesPairSourceVector(beta,Nf,im_i,im_j,phase_f,vec_rot)
%   rotateStokesPairSourceVector()   % run built-in self-test
%
% Ordering:
%   [f_i_x; e_i_x; f_j_x; e_j_x; f_i_y; e_i_y; f_j_y; e_j_y]

if nargin == 0
    run_self_test();
    if nargout > 0
        beta_rot = [];
    end
    return
end

if nargin < 6 || isempty(vec_rot)
    vec_rot = 1;
end

beta = beta(:);
if isempty(beta)
    beta_rot = beta;
    return
end

f_ind1_x = 1:Nf;
e_ind1_x = Nf+1:Nf+im_i;
f_ind2_x = Nf+im_i+1:2*Nf+im_i;
e_ind2_x = 2*Nf+im_i+1:2*Nf+im_i+im_j;
f_ind1_y = 2*Nf+im_i+im_j+1:3*Nf+im_i+im_j;
e_ind1_y = 3*Nf+im_i+im_j+1:3*Nf+2*im_i+im_j;
f_ind2_y = 3*Nf+2*im_i+im_j+1:4*Nf+2*im_i+im_j;
e_ind2_y = 4*Nf+2*im_i+im_j+1:4*Nf+2*im_i+2*im_j;

fine_pair = [beta(f_ind1_x) beta(f_ind2_x) beta(f_ind1_y) beta(f_ind2_y)];
%fine_pair = beta;

% Rotate the fine data on the periodic grids
fine_pair = rotatePairOrderedStokesData(fine_pair,Nf,phase_f,vec_rot);

ei = vec_rot*(beta(e_ind1_x) + 1i*beta(e_ind1_y));
ej = vec_rot*(beta(e_ind2_x) + 1i*beta(e_ind2_y));

beta_rot = zeros(size(beta));
beta_rot(f_ind1_x) = fine_pair(:,1);
beta_rot(f_ind2_x) = fine_pair(:,2);
beta_rot(f_ind1_y) = fine_pair(:,3);
beta_rot(f_ind2_y) = fine_pair(:,4);
beta_rot(e_ind1_x) = real(ei);
beta_rot(e_ind1_y) = imag(ei);
beta_rot(e_ind2_x) = real(ej);
beta_rot(e_ind2_y) = imag(ej);
end

function run_self_test()
repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

rng(1);

opt = get2Dparams(3);
opt.rad = 1;
opt.tol_c = 1e-12;
opt.use_fmm = false;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.get_bndry_field = 0;
opt.precomp = 1;
opt.cmap = 1;
opt.N_c = 16;
opt.N_f = 24;
opt.N_peanut = 0;
opt.Nclust = 60;
opt.delta_pair = 0.2;
opt.a_c = 1.2;
opt.a_f = 1.2;
opt.project_force = true;
opt.reuse_pair_basis_by_sep = true;

geoms = { ...
    struct('name','collinear', ...
           'q',[0; 2.001; 4.002]), ...
    struct('name','equilateral', ...
           'q',[0; 2.001; 0.5*2.001 + 1i*(sqrt(3)/2*2.001)])};

tol = 1e-12;
max_err = 0;

fprintf('rotateStokesPairSourceVector self-test:\n');
for g = 1:numel(geoms)
    data = build_self_test_case(geoms{g}.q,opt);
    fprintf('  Geometry: %s\n',geoms{g}.name);
    for row = 1:size(data.pairs,1)
        [fine_err,e_err] = check_pair_rotation(data,row);
        max_err = max([max_err, fine_err, e_err]);
        i = data.pairs(row,1);
        j = data.pairs(row,2);
        theta = atan2(imag(data.pair_cache.meta(row).rot), ...
            real(data.pair_cache.meta(row).rot));
        fprintf('    row %d pair (%d,%d), angle %.6f rad: fine %.3e, e %.3e\n', ...
            row,i,j,theta,fine_err,e_err);
    end
end

fprintf('  max field-equivariance error: %.3e\n',max_err);
if max_err > tol
    error('rotateStokesPairSourceVector:selftest', ...
        'Self-test failed with max error %.3e.',max_err);
end
end

function data = build_self_test_case(q,opt_in)
q = q(:);
P = numel(q);
opt = opt_in;
opt.P = P;

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

N_c = opt.N_c;
N_f = opt.N_f;
tol_c = opt.tol_c;
rad = opt.rad;

tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
sep_c = (1/N_c)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',rad*max([1-sep_c,0.01]));
rbase_in_c = Rp_c*cos(tin_c) + 1i*Rp_c*sin(tin_c);

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
sep_f = (1/N_f)*log(1/tol_c);
Rp_f = getOptField(opt,'Rp_f',rad*max([1-sep_f,0.01]));
rbase_in_f = Rp_f*cos(tin_f) + 1i*Rp_f*sin(tin_f);

[~,~,Lc] = getSelfPseudoMobilityStokes(rad,q,rbase_in_c, ...
    rad*(cos(linspace(0,2*pi,ceil(opt.a_c*N_c)+1)).' + ...
    1i*sin(linspace(0,2*pi,ceil(opt.a_c*N_c)+1)).'),[],[0,ceil(opt.a_c*N_c)]);
[~,~,~,~,~,~,pair_cache] = getPairBasisStokes(q,rbase_in_c,rbase_in_f, ...
    rimage_vec,refine,pairs,opt,Lc{1});

data = struct();
data.q = q;
data.opt = opt;
data.pairs = pairs;
data.rbase_in_f = rbase_in_f;
data.pair_cache = pair_cache;
end

function [fine_err,e_err] = check_pair_rotation(data,row)
pair = getStokesPairInstance(data.pair_cache,row);
idx = get_pair_source_indices(data.opt.N_f,numel(pair.rimage_i),numel(pair.rimage_j));

targets_canon = [pair.group.refine_canon{1}; pair.group.refine_canon{2}];
targets_actual = [pair.refine_i; pair.refine_j];

beta_f_loc = zeros(idx.ntot,1);
beta_f_loc(idx.f1x) = randn(numel(idx.f1x),1);
beta_f_loc(idx.f2x) = randn(numel(idx.f2x),1);
beta_f_loc(idx.f1y) = randn(numel(idx.f1y),1);
beta_f_loc(idx.f2y) = randn(numel(idx.f2y),1);
    beta_f_rot = rotateStokesPairSourceVector(beta_f_loc,data.opt.N_f, ...
        numel(pair.rimage_i),numel(pair.rimage_j),pair.meta.phase_f_inv,pair.meta.rot);

src_f_canon = [pair.group.q_pair(1)+data.rbase_in_f; ...
               pair.group.q_pair(2)+data.rbase_in_f];
src_f_actual = [data.q(pair.meta.i)+data.rbase_in_f; ...
                data.q(pair.meta.j)+data.rbase_in_f];
field_f_canon = eval_stokes_field(src_f_canon, ...
    [beta_f_loc(idx.f1x); beta_f_loc(idx.f2x)], ...
    [beta_f_loc(idx.f1y); beta_f_loc(idx.f2y)],targets_canon);
field_f_actual = eval_stokes_field(src_f_actual, ...
    [beta_f_rot(idx.f1x); beta_f_rot(idx.f2x)], ...
    [beta_f_rot(idx.f1y); beta_f_rot(idx.f2y)],targets_actual);
fine_err = norm(field_f_actual - pair.meta.rot*field_f_canon) / ...
    max(norm(field_f_actual),eps);

beta_e_loc = zeros(idx.ntot,1);
beta_e_loc(idx.e1x) = randn(numel(idx.e1x),1);
beta_e_loc(idx.e2x) = randn(numel(idx.e2x),1);
beta_e_loc(idx.e1y) = randn(numel(idx.e1y),1);
beta_e_loc(idx.e2y) = randn(numel(idx.e2y),1);
    beta_e_rot = rotateStokesPairSourceVector(beta_e_loc,data.opt.N_f, ...
        numel(pair.rimage_i),numel(pair.rimage_j),pair.meta.phase_f_inv,pair.meta.rot);

src_e_canon = [pair.group.rimage_canon{1}; pair.group.rimage_canon{2}];
src_e_actual = [pair.rimage_i; pair.rimage_j];
field_e_canon = eval_stokes_field(src_e_canon, ...
    [beta_e_loc(idx.e1x); beta_e_loc(idx.e2x)], ...
    [beta_e_loc(idx.e1y); beta_e_loc(idx.e2y)],targets_canon);
field_e_actual = eval_stokes_field(src_e_actual, ...
    [beta_e_rot(idx.e1x); beta_e_rot(idx.e2x)], ...
    [beta_e_rot(idx.e1y); beta_e_rot(idx.e2y)],targets_actual);
e_err = norm(field_e_actual - pair.meta.rot*field_e_canon) / ...
    max(norm(field_e_actual),eps);
end

function idx = get_pair_source_indices(Nf,im_i,im_j)
idx = struct();
idx.f1x = 1:Nf;
idx.e1x = Nf+1:Nf+im_i;
idx.f2x = Nf+im_i+1:2*Nf+im_i;
idx.e2x = 2*Nf+im_i+1:2*Nf+im_i+im_j;
idx.f1y = 2*Nf+im_i+im_j+1:3*Nf+im_i+im_j;
idx.e1y = 3*Nf+im_i+im_j+1:3*Nf+2*im_i+im_j;
idx.f2y = 3*Nf+2*im_i+im_j+1:4*Nf+2*im_i+im_j;
idx.e2y = 4*Nf+2*im_i+im_j+1:4*Nf+2*im_i+2*im_j;
idx.ntot = 4*Nf + 2*(im_i+im_j);
end

function field = eval_stokes_field(zsrc,fx,fy,ztar)
if isempty(zsrc)
    field = zeros(numel(ztar),1);
    return
end
[u,v] = stokSLPdirect(real(zsrc),imag(zsrc),real(ztar),imag(ztar),fx,fy,numel(zsrc));
field = u + 1i*v;
end
