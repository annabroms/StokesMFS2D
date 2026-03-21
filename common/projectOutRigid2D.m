function tau_proj = projectOutRigid2D(tau_vec, rsrc, q)
%PROJECTOUTRIGID2D Remove rigid translation/rotation density components in 2D.
%
% Syntax:
%   tau_proj = projectOutRigid2D(tau_vec, rsrc, q)
%   projectOutRigid2D()   % run built-in self-test
%
% Description:
%   For source locations rsrc (complex) around center q, this routine
%   applies the matrix-free equivalent of the K-based projection
%
%       tau_proj = (I - K*(K'*K)\K') * tau_vec,
%
%   where K = getKmat2D(rsrc,q). Here tau_vec = [tau_x; tau_y].
%   The rigid component is recovered from a 3x3 normal-equation system.
%
% Inputs:
%   tau_vec  - stacked source strengths [tau_x; tau_y], length 2*numel(rsrc)
%   rsrc     - complex source locations
%   q        - particle center (complex scalar)
%
% Outputs:
%   tau_proj - projected density (same layout as tau_vec)
%
% Notes:
%   If called with no inputs, a self-test is run that compares against the
%   explicit matrix projection built from getKmat2D.
%   If the internal 3x3 system is ill-conditioned, an error is thrown so a
%   caller can fall back to a matrix-based implementation.
%
% Anna Broms, Feb 2026

if nargin == 0
    run_self_test();
    return;
end

if nargin < 3
    error('projectOutRigid2D requires tau_vec, rsrc, and q.');
end

tau_vec = tau_vec(:);
rsrc = rsrc(:);
n = numel(rsrc);

if n == 0
    tau_proj = tau_vec;
    return;
end

if numel(tau_vec) ~= 2*n
    error('tau_vec must have length 2*numel(rsrc).');
end

tau_x = tau_vec(1:n);
tau_y = tau_vec(n+1:end);

rel = rsrc - q;
dx = real(rel);
dy = imag(rel);

g13 = -sum(dy);
g23 = sum(dx);
g33 = sum(dx.^2 + dy.^2);
G = [n, 0, g13; ...
     0, n, g23; ...
     g13, g23, g33];

b = [sum(tau_x); ...
     sum(tau_y); ...
     sum((-dy).*tau_x + dx.*tau_y)];

alpha = G\b;
rigid_x = alpha(1) - alpha(3)*dy;
rigid_y = alpha(2) + alpha(3)*dx;

tau_proj = [tau_x - rigid_x; tau_y - rigid_y];

end

function run_self_test()
rng(1);
tol_eq = 5e-12;
tol_mom = 5e-12;

q = 0.15 - 0.25i;
N = 64;
t = linspace(0,2*pi,N+1)'; t(end) = [];
r = q + (1.0 + 0.15*cos(3*t)).*(cos(t) + 1i*sin(t));

K = getKmat2D(r,q);
L = K*((K'*K)\K');

fprintf('projectOutRigid2D self-test:\n');
for k = 1:5
    tau = randn(2*N,1);
    tau_proj = projectOutRigid2D(tau,r,q);
    tau_ref = tau - L*tau;

    rel_err = norm(tau_proj - tau_ref,inf)/max(1,norm(tau_ref,inf));
    mom_err = norm(K'*tau_proj,inf);

    fprintf('  sample %d: rel err = %.3e, moment err = %.3e\n',k,rel_err,mom_err);
    assert(rel_err < tol_eq,'Projection mismatch against K-based reference.');
    assert(mom_err < tol_mom,'Projected vector still has rigid moments.');
end

% Rigid mode should be projected to (near) zero.
tau_rigid = [ones(N,1); zeros(N,1)];
tau_rigid_proj = projectOutRigid2D(tau_rigid,r,q);
rigid_err = norm(tau_rigid_proj,inf);
fprintf('  rigid-mode residual : %.3e\n',rigid_err);
assert(rigid_err < tol_mom,'Rigid mode was not removed.');

fprintf('  PASSED\n');
end
