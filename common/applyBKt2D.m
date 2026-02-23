function u = applyBKt2D(r_tar, q_tar, r_src, q_src, f_x, f_y)
%APPLYBKT2D Apply B*K'*[f_x; f_y] without building the RBM K and B matrices, 
% 
% Remember K^T maps from sources to net force/torque and B maps from RBM data to surface velocities. 
%
% Syntax:
%   u = applyBKt2D(r_tar, q_tar, r_src, q_src, f_x, f_y)
%   applyBKt2D()   % run self-test
%
% Description:
%   Computes the action
%
%       u = K_tar * (K_src' * [f_x; f_y]),
%
%   where K_src = getKmat2D(r_src, q_src) and K_tar = getKmat2D(r_tar, q_tar),
%   for one particle in 2D. This is exactly the same map as
%   B*K'*[f_x; f_y], but evaluated matrix-free via total force/torque moments.
%
% Inputs:
%   r_tar - complex target points (column vector)
%   q_tar - complex target center (scalar)
%   r_src - complex source points (column vector)
%   q_src - complex source center (scalar)
%   f_x   - x-component source strengths at r_src
%   f_y   - y-component source strengths at r_src
%
% Output:
%   u     - stacked velocity-like vector [u_x; u_y] at r_tar
%
% Anna Broms, Feb 2026

if nargin == 0
    run_self_test();
    return;
end

r_tar = r_tar(:);
r_src = r_src(:);
f_x = f_x(:);
f_y = f_y(:);

rel_src = r_src - q_src;
Fx = sum(f_x);
Fy = sum(f_y);
Tz = sum((-imag(rel_src)).*f_x + real(rel_src).*f_y);

rel_tar = r_tar - q_tar;
u_x = Fx + Tz*(-imag(rel_tar));
u_y = Fy + Tz*( real(rel_tar));
u = [u_x; u_y];

end

function run_self_test()
rng(2);
tol = 5e-12;

fprintf('applyBKt2D self-test:\n');

n_src = 200;
n_tar = 300; 

t_src = linspace(0,2*pi,n_src+1)'; t_src(end) = [];
t_tar = linspace(0,2*pi,n_tar+1)'; t_tar(end) = [];

r_src = (cos(t_src) + 1i*sin(t_src));
r_tar = (cos(t_tar) + 1i*sin(t_tar));

f_x = randn(n_src,1);
f_y = randn(n_src,1);

u_fast = applyBKt2D(r_tar,0,r_src,0,f_x,f_y);

K_src = getKmat2D(r_src,0);
K_tar = getKmat2D(r_tar,0);
u_ref = K_tar*(K_src'*[f_x; f_y]);

rel_err = norm(u_fast-u_ref,inf)/max(1,norm(u_ref,inf));
fprintf('  rel.err = %.3e\n', rel_err);
assert(rel_err < tol, 'applyBKt2D mismatch against matrix-based reference.');


fprintf('  PASSED\n');
end
