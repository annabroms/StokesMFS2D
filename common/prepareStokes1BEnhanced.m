function [geom,basis,opt,rad] = prepareStokes1BEnhanced(q,opt)
%PREPARESTOKES1BENHANCED Build 1-body Stokeslet-only enhanced geometry and basis.
%
% Syntax:
%   [geom,basis,opt,rad] = prepareStokes1BEnhanced(q,opt)
%
% Inputs:
%   q   - Complex particle centers (P x 1).
%   opt - Options struct.
%
% Outputs:
%   geom  - Geometry struct with bodywise source/target blocks.
%   basis - Per-body pseudoinverse factors for the enhanced 1B blocks.
%   opt   - Options struct with derived defaults filled in.
%   rad   - Physical body radius used for all particles.
%
% Notes:
%   This helper uses the ellipse-segment enhancement nodes produced by
%   getEnhancedGrid and keeps the solver Stokeslet-only. It currently
%   assumes equal body radii.

if nargin < 2 || ~isstruct(opt)
    error('prepareStokes1BEnhanced requires q and an options struct opt.');
end

q = q(:);
P = numel(q);

rad = 1;
rad0 = rad(1);

N_c = getOptField(opt,'N_c',100);
a_c = getOptField(opt,'a_c',1.2);
tol_c = getOptField(opt,'tol_c',1e-12);
self_tol = getOptField(opt,'self_tol',1e-11);
delta_pair = getOptField(opt,'delta_pair',0.2);
Nclust = getOptField(opt,'Nclust',100);
beta = getOptField(opt,'beta',0.3);
use_fmm = logical(getOptField(opt,'use_fmm',true));
visualise_grid = logical(getOptField(opt,'visualise_grid',false));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
if isfield(opt,'solve_resistance')
    solve_resistance = logical(opt.solve_resistance);
elseif isfield(opt,'one_body_mode')
    solve_resistance = strcmpi(char(opt.one_body_mode),'resistance');
else
    solve_resistance = false;
end

sep_c = (1/N_c)*log(1/tol_c);

Rp_c = getOptField(opt,'Rp_c',rad0*max([1-sep_c,0.01]));
if isscalar(Rp_c)
    Rp_c = repmat(Rp_c,P,1);
else
    Rp_c = Rp_c(:);
    if numel(Rp_c) ~= P
        error('prepareStokes1BEnhanced:badRpC', ...
            'opt.Rp_c must be scalar or have one entry per particle.');
    end
    if any(abs(Rp_c-Rp_c(1)) > 1e-12)
        error('prepareStokes1BEnhanced:nonuniformRpC', ...
            ['prepareStokes1BEnhanced currently assumes equal proxy ', ...
             'radii for all particles.']);
    end
    Rp_c = repmat(Rp_c(1),P,1);
end

% Fill the commonly used fields so downstream helpers can rely on them.
opt_clean = opt;
for fld = {'N_f','a_f','Rp_f','one_body_mode'}
    if isfield(opt_clean,fld{1})
        opt_clean = rmfield(opt_clean,fld{1});
    end
end

opt_clean.rad = rad0;
opt_clean.N_c = N_c;
opt_clean.a_c = a_c;
opt_clean.Rp_c = Rp_c(1);
opt_clean.delta_pair = delta_pair;
opt_clean.Nclust = Nclust;
opt_clean.beta = beta;
opt_clean.use_fmm = use_fmm;
opt_clean.visualise_grid = visualise_grid;
opt_clean.gmres_verbose = gmres_verbose;
opt_clean.solve_resistance = solve_resistance;
opt_clean.image = 0;
opt_clean.pc = 0;

opt = opt_clean;
opt_grid = opt_clean;
opt_grid.r_proxy = Rp_c(1);
[cent_clust_cells,~,coll_clust_cells,~,~,pairs] = getEnhancedGrid(q,opt_grid);

nout = ceil(a_c*N_c);
tout = linspace(0,2*pi,nout+1)';
tout = tout(1:end-1);
rbase_out_c = rad0*(cos(tout)+1i*sin(tout));

tin = linspace(0,2*pi,N_c+1)';
tin = tin(1:end-1);
rbase_in_c = Rp_c(1)*(cos(tin)+1i*sin(tin));

rin_body = cell(P,1);
rout_body = cell(P,1);
source_count = zeros(P,1);
target_count = zeros(P,1);

for k = 1:P
    rin_body{k} = [q(k) + rbase_in_c; cent_clust_cells{k}];
    rout_body{k} = [q(k) + rbase_out_c; coll_clust_cells{k}];
    source_count(k) = numel(rin_body{k});
    target_count(k) = numel(rout_body{k});
end

source_offsets = cumsum([1; source_count(1:end-1)]);
target_offsets = cumsum([1; target_count(1:end-1)]);
total_source_count = sum(source_count);
total_target_count = sum(target_count);

source_x_ind = cell(P,1);
source_y_ind = cell(P,1);
target_x_ind = cell(P,1);
target_y_ind = cell(P,1);
Ksrc = cell(P,1);
Aii = cell(P,1);
U = cell(P,1);
Y = cell(P,1);

for k = 1:P
    ns = source_count(k);
    nt = target_count(k);
    sx = source_offsets(k):source_offsets(k)+ns-1;
    tx = target_offsets(k):target_offsets(k)+nt-1;

    source_x_ind{k} = sx;
    source_y_ind{k} = sx;
    target_x_ind{k} = tx;
    target_y_ind{k} = total_target_count + tx;

    Ksrc{k} = getKmat2D(rin_body{k},q(k));
    Nio = stokSLPmat(rin_body{k},rout_body{k},1);

    if solve_resistance
        Aii{k} = Nio;
        [Yk,Uk] = getPseudoFactors(Nio,self_tol,1);
    else
        Ktar = getKmat2D(rout_body{1},q(1));
        L = Ksrc{k}*((Ksrc{k}'*Ksrc{k})\Ksrc{k}');
        A = Nio*(eye(size(L)) - L);
        [Yk,Uk] = getPseudoFactors(A+Ktar*Ksrc{k}',self_tol,1);
        Aii{k} = A;
    end

    
    U{k} = Uk';
    Y{k} = Yk;
end

geom = struct();
geom.q = q;
geom.rad = rad0;
geom.rbase_in_c = rbase_in_c;
geom.rbase_out_c = rbase_out_c;
geom.rin_body = rin_body;
geom.rout_body = rout_body;
geom.source_count = source_count;
geom.target_count = target_count;
geom.source_offsets = source_offsets;
geom.target_offsets = target_offsets;
geom.source_x_ind = source_x_ind;
geom.source_y_ind = source_y_ind;
geom.target_x_ind = target_x_ind;
geom.target_y_ind = target_y_ind;
geom.total_source_count = total_source_count;
geom.total_target_count = total_target_count;
geom.N_c = N_c;
geom.rvec_in = vertcat(rin_body{:});
geom.rout = vertcat(rout_body{:});
geom.pairs = pairs;
geom.Ksrc = Ksrc;
geom.cent_clust_cells = cent_clust_cells;
geom.coll_clust_cells = coll_clust_cells;
geom.use_fmm = use_fmm;
geom.solve_resistance = solve_resistance;
geom.opt = opt_clean;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Aii = Aii;
basis.Ksrc = Ksrc;
basis.self_tol = self_tol;
basis.use_fmm = use_fmm;
basis.solve_resistance = solve_resistance;

end
