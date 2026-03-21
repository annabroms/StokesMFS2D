function pair_proj = projectOutRigidPair2D(pair_vec, rsrc_i, q_i, rsrc_j, q_j)
%PROJECTOUTRIGIDPAIR2D Matrix-free pair analogue of getILpair(L).
%
% Syntax:
%   pair_proj = projectOutRigidPair2D(pair_vec, rsrc_i, q_i, rsrc_j, q_j)
%   projectOutRigidPair2D()   % run built-in self-test
%
% Input ordering:
%   pair_vec = [x_i; x_j; y_i; y_j]
%
% The two bodies are projected independently with projectOutRigid2D and
% then reassembled into pair ordering. This is the matrix-free equivalent
% of multiplying by getILpair(L) for the single-body projector L.
%
% Anna Broms, Mar 2026

if nargin == 0
    run_self_test();
    return;
end

if nargin < 5
    error('projectOutRigidPair2D requires pair_vec, rsrc_i, q_i, rsrc_j, and q_j.');
end

pair_vec = pair_vec(:);
rsrc_i = rsrc_i(:);
rsrc_j = rsrc_j(:);

n = numel(rsrc_i);
if numel(rsrc_j) ~= n
    error('projectOutRigidPair2D:SizeMismatch', ...
        'The two bodies must have the same number of source points.');
end
if numel(pair_vec) ~= 4*n
    error('projectOutRigidPair2D:SizeMismatch', ...
        'pair_vec must have length 4*numel(rsrc_i).');
end

proj_i = projectOutRigid2D([pair_vec(1:n); pair_vec(2*n+1:3*n)], rsrc_i, q_i);
proj_j = projectOutRigid2D([pair_vec(n+1:2*n); pair_vec(3*n+1:4*n)], rsrc_j, q_j);

pair_proj = [proj_i(1:n); proj_j(1:n); proj_i(n+1:end); proj_j(n+1:end)];
end

function run_self_test()
rng(2);
tol_eq = 5e-12;
tol_mom = 5e-12;

q_i = -0.35 + 0.15i;
q_j =  0.40 - 0.10i;
N = 48;
t = linspace(0,2*pi,N+1)';
t(end) = [];
r_i = q_i + 0.9*(cos(t) + 1i*sin(t));
r_j = q_j + 0.9*(cos(t) + 1i*sin(t));

K = getKmat2D(r_i,q_i);
L = K*((K'*K)\K');
ILpair = getILpair(L);

fprintf('projectOutRigidPair2D self-test:\n');
for k = 1:5
    pair_vec = randn(4*N,1);
    pair_ref = ILpair*pair_vec;
    pair_proj = projectOutRigidPair2D(pair_vec,r_i,q_i,r_j,q_j);

    rel_err = norm(pair_proj - pair_ref,inf)/max(1,norm(pair_ref,inf));
    mom_i = norm(K'*[pair_proj(1:N); pair_proj(2*N+1:3*N)],inf);
    mom_j = norm(K'*[pair_proj(N+1:2*N); pair_proj(3*N+1:4*N)],inf);

    fprintf('  sample %d: rel err = %.3e, mom_i = %.3e, mom_j = %.3e\n', ...
        k,rel_err,mom_i,mom_j);
    assert(rel_err < tol_eq,'Pair rigid projection mismatch against dense reference.');
    assert(mom_i < tol_mom && mom_j < tol_mom, ...
        'Projected pair vector still has rigid moments.');
end

fprintf('  PASSED\n');
end
