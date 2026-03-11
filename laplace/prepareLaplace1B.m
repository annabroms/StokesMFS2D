function [geom,basis,opt,rad] = prepareLaplace1B(q,use_fmm,visualise,project_charge)
%PREPARELAPLACE1B Build shared one-body Laplace geometry and pseudo-inverse factors.
%
% Syntax:
%   [geom,basis,opt,rad] = prepareLaplace1B(q,use_fmm)
%   [geom,basis,opt,rad] = prepareLaplace1B(q,use_fmm,visualise,project_charge)
%
% Inputs:
%   q              - Complex particle centers (P x 1).
%   use_fmm        - Logical flag propagated to field evaluations.
%   visualise      - Optional logical flag for enhanced-grid visualisation.
%   project_charge - Optional logical flag. If true, build elastance-style
%                    one-body factors from S(I-Lq)+Lr; otherwise use S.
%
% Outputs:
%   geom, basis - One-body geometry and pseudo-inverse data.
%   opt         - Laplace option struct used for grid construction.
%   rad         - Physical particle radius.
%
% See also: solve_cap_1B, solve_elast_1B, getEnhancedGrid.
%
% Anna Broms, Mar 2026

if nargin < 3 || isempty(visualise)
    visualise = false;
end
if nargin < 4 || isempty(project_charge)
    project_charge = false;
end

q = q(:);
P = numel(q);

opt = getLaplace2Dparams();
rad = opt.rad;

N_c = 80;
N_f = 150;
a_c = 1.2;
a_f = 1.2;

tol_c = 1e-10;
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = rad*max([1-sep_c,0.01]);
Rp_f = rad*max([1-sep_f,0.01]);

accstop = (rad-Rp_c)^2/Rp_c;

opt.Rp_c = Rp_c;
opt.Rp_f = Rp_f;
opt.a_c = a_c;
opt.a_f = a_f;
opt.N_c = N_c;
opt.N_f = N_f;
opt.N_peanut = 0;
opt.precomp = 1;
opt.pc = 0;
opt.delta_pair = accstop;
opt.P = P;
opt.Nclust = 100;
opt.cmap = 0;
opt.use_fmm = use_fmm;
opt.show_counter = false;
opt.visualise_grid = logical(visualise);
opt.project_charge = logical(project_charge);
opt.rads = rad*ones(P,1);

nout = ceil(a_c*N_c);
tout = linspace(0,2*pi,nout+1)';
tout = tout(1:end-1);
rbase_out_c = rad*(cos(tout)+1i*sin(tout));

tin = linspace(0,2*pi,N_c+1)';
tin = tin(1:end-1);
rbase_in_c = Rp_c*(cos(tin)+1i*sin(tin));

[cent_clust_cells,~,coll_clust_cells,~,~,pairs] = getEnhancedGrid(q,opt);

rin_body = cell(P,1);
rout_body = cell(P,1);
source_ind = cell(P,1);
target_ind = cell(P,1);
Aii = cell(P,1);
U = cell(P,1);
Y = cell(P,1);

r_chunks = cell(P,1);
t_chunks = cell(P,1);
source_start = 1;
target_start = 1;

for k = 1:P
    rin_k = [q(k)+rbase_in_c; cent_clust_cells{k}];
    rout_k = [q(k)+rbase_out_c; coll_clust_cells{k}];

    rin_body{k} = rin_k;
    rout_body{k} = rout_k;

    ns = length(rin_k);
    nt = length(rout_k);

    source_ind{k} = source_start:source_start+ns-1;
    target_ind{k} = target_start:target_start+nt-1;
    source_start = source_start+ns;
    target_start = target_start+nt;

    SLP = lapSLPmat(rin_k,rout_k);
    if project_charge
        Aii{k} = apply_charge_projection_with_closure(SLP,ns);
        [Yk,Uk] = getPseudoFactors(Aii{k}+1,1e-14,0);
    else
        Aii{k} = SLP;
        [Yk,Uk] = getPseudoFactors(Aii{k},1e-14,0);
    end

    U{k} = Uk';
    Y{k} = Yk;

    r_chunks{k} = rin_k;
    t_chunks{k} = rout_k;
end

geom = struct();
geom.q = q;
geom.rvec_in = vertcat(r_chunks{:});
geom.rout = vertcat(t_chunks{:});
geom.rin_body = rin_body;
geom.rout_body = rout_body;
geom.source_ind = source_ind;
geom.target_ind = target_ind;
geom.pairs = pairs;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Aii = Aii;
basis.use_fmm = use_fmm;

end

function Aproj = apply_charge_projection_with_closure(A,n)
% Apply A*(I-Lq) for Kq=ones(n,1), then add Lr.
Aproj = A;
row_mean = sum(Aproj,2)/n;
Aproj = Aproj - row_mean*ones(1,n);
end
