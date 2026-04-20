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
% The two bodies are projected independently and reassembled into pair
% ordering. This is the matrix-free equivalent of multiplying by
% getILpair(L) for the single-body projector L.
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

if n == 0
    pair_proj = pair_vec;
    return
end

rel_i = rsrc_i - q_i;
rel_j = rsrc_j - q_j;

x_i = pair_vec(1:n);
x_j = pair_vec(n+1:2*n);
y_i = pair_vec(2*n+1:3*n);
y_j = pair_vec(3*n+1:4*n);

if isequal(rel_i,rel_j)
    dx = real(rel_i);
    dy = imag(rel_i);
    G = rigid_projection_gram(dx,dy,n);
    alpha = G \ [sum(x_i), sum(x_j); ...
                 sum(y_i), sum(y_j); ...
                 -dy.'*x_i + dx.'*y_i, -dy.'*x_j + dx.'*y_j];

    pair_proj = pair_vec;
    pair_proj(1:n) = x_i - alpha(1,1) + alpha(3,1)*dy;
    pair_proj(n+1:2*n) = x_j - alpha(1,2) + alpha(3,2)*dy;
    pair_proj(2*n+1:3*n) = y_i - alpha(2,1) - alpha(3,1)*dx;
    pair_proj(3*n+1:4*n) = y_j - alpha(2,2) - alpha(3,2)*dx;
else
    [proj_i_x,proj_i_y] = project_body_components(x_i,y_i,rel_i);
    [proj_j_x,proj_j_y] = project_body_components(x_j,y_j,rel_j);
    pair_proj = [proj_i_x; proj_j_x; proj_i_y; proj_j_y];
end
end

function G = rigid_projection_gram(dx,dy,n)
g13 = -sum(dy);
g23 = sum(dx);
g33 = sum(dx.^2 + dy.^2);
G = [n, 0, g13; ...
     0, n, g23; ...
     g13, g23, g33];
end

function [proj_x,proj_y] = project_body_components(tau_x,tau_y,rel)
dx = real(rel);
dy = imag(rel);
n = numel(rel);
G = rigid_projection_gram(dx,dy,n);
alpha = G \ [sum(tau_x); ...
             sum(tau_y); ...
             -dy.'*tau_x + dx.'*tau_y];

proj_x = tau_x - alpha(1) + alpha(3)*dy;
proj_y = tau_y - alpha(2) - alpha(3)*dx;
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
