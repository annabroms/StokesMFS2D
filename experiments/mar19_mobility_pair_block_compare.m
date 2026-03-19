clear;
close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Mobility Pair-Block Comparison (Mar 19, 2026) ===\n');

rad = 1;
q_pair = [0; 2.2];

N_c = 150;
N_f = 150;
a_c = 1.2;
a_f = 1.2;
tol_c = 1e-12;
N_peanut = 200;
s = [0 0 1 1 0 0 0];

[nout_c,Rp_c,rbase_in_c,rbase_out_c] = build_circle_grids(rad,N_c,a_c,tol_c);
[nout_f,Rp_f,rbase_in_f,rout_base_f] = build_circle_grids(rad,N_f,a_f,tol_c);

fprintf('Single pair geometry: q = [%.3f, %.3f], gap = %.3f\n', ...
    q_pair(1),q_pair(2),abs(q_pair(2)-q_pair(1))-2*rad);
fprintf('  N_c=%d, N_f=%d, N_peanut=%d, Rp_c=%.4f, Rp_f=%.4f\n', ...
    N_c,N_f,N_peanut,Rp_c,Rp_f);

Kc = getKmat2D(rbase_in_c,0);
Lc_body = Kc*((Kc'*Kc)\Kc');
Kf = getKmat2D(rbase_in_f,0);
Lf_body = Kf*((Kf'*Kf)\Kf');

Nio_self = stokSLPmat(rbase_in_c,rbase_out_c,1);
Lr_self = getKmat2D(rbase_out_c,0) * Kc';
Nii_self = Nio_self * (eye(size(Lc_body)) - Lc_body) + Lr_self;

fprintf('\nSelf map:\n');
fprintf('  raw operator rank @1e-10 = %d, @1e-14 = %d\n', ...
    svd_rank(Nii_self,1e-10),svd_rank(Nii_self,1e-14));
fprintf('  for this geometry the raw self operator is rank-stable across those cutoffs.\n');

% Fine pair projection conventions.
Lf_pair_old = getILpair(Lf_body);
Kf1 = getKmat2D(q_pair(1)+rbase_in_f,q_pair(1));
Kf2 = getKmat2D(q_pair(2)+rbase_in_f,q_pair(2));
Lf_pair_new = getLfPair(Kf1,Kf2);

pair_proj_err = relerr(Lf_pair_old,eye(size(Lf_pair_new)) - Lf_pair_new);

rin_pair_f = [q_pair(1)+rbase_in_f; q_pair(2)+rbase_in_f];
rout_f = [q_pair(1)+rout_base_f; q_pair(2)+rout_base_f];
B1 = getKmat2D(q_pair(1)+rout_base_f,q_pair(1));
B2 = getKmat2D(q_pair(2)+rout_base_f,q_pair(2));
Lr_pair = getLrPair(B1,B2,Kf1,Kf2);

[~,~,N_pair_old] = getPairBlock( ...
    q_pair,rin_pair_f,rout_f,[],[],s,Lf_pair_old,Lr_pair,0);
[~,~,N_pair_new] = getPairBlockStokes( ...
    rin_pair_f,rout_f,Lf_pair_new,Lr_pair);

fprintf('\nFine pair block:\n');
fprintf('  projection complement relerr = %.3e\n', pair_proj_err);
fprintf('  raw operator relerr old/new   = %.3e\n', relerr(N_pair_old,N_pair_new));
fprintf('  raw operator rank @1e-14 = %d, @1e-11 = %d\n', ...
    svd_rank(N_pair_new,1e-14),svd_rank(N_pair_new,1e-11));
fprintf('  for this geometry the raw pair operator is also rank-stable here.\n');

% Peanut compression.
rout_peanut = createPeanut(q_pair(1),q_pair(2),N_peanut,0,rad);
rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
Lc_pair = getILpair(Lc_body);

Nf_peanut = stokSLPmat(rin_pair_f,rout_peanut,1);
Npeanut_old = Nf_peanut * Lf_pair_old;
Npeanut_new = Nf_peanut - Nf_peanut * Lf_pair_new;

[~,~] = getPeanutBlock( ...
    rin_pair_c,rin_pair_f,rout_peanut,[],[],[],Lc_pair,Lf_pair_old);
[~,~] = getPeanutBlockStokes( ...
    rin_pair_c,rin_pair_f,rout_peanut,Lc_pair,Lf_pair_new);

fprintf('\nPeanut block:\n');
fprintf('  fine projection complement relerr = %.3e\n', relerr(Npeanut_old,Npeanut_new));
fprintf('  raw operator rank @1e-14 = %d, @1e-11 = %d\n', ...
    svd_rank(Npeanut_new,1e-14),svd_rank(Npeanut_new,1e-11));

fprintf('\nDone.\n');

function [nout,Rp,rbase_in,rbase_out] = build_circle_grids(rad,N,a_c,tol_c)
nout = ceil(a_c*N);
sep = (1/N) * log(1/tol_c);
Rp = rad * max([1-sep,0.01]);

tin = linspace(0,2*pi,N+1)';
tin = tin(1:end-1);
rbase_in = Rp * (cos(tin) + 1i*sin(tin));

tout = linspace(0,2*pi,nout+1)';
tout = tout(1:end-1);
rbase_out = rad * (cos(tout) + 1i*sin(tout));
end

function e = relerr(A,B)
e = norm(A-B,'fro') / max(1,norm(B,'fro'));
end

function r = svd_rank(A,tol)
s = svd(A,'econ');
r = sum(s > max(s)*tol);
end
