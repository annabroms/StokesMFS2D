function A = lapSLPmat(rsrc,rtar)
%LAPSLPMAT Dense target-from-source matrix for 2D Laplace SLP.
%
% Syntax:
%   A = lapSLPmat(rsrc,rtar)
%
% Inputs:
%   rsrc - Complex source locations.
%   rtar - Complex target locations.
%
% Output:
%   A    - Dense matrix such that u = A*sigma corresponds to
%          u(x) = -(1/(2*pi)) sum_j sigma_j log|x-x_j|.
%
% Self-test:
%   lapSLPmat()
%
% See also: lapSLPfield, lapSLPdirect.
%
% Anna Broms, Mar 2026

if nargin==0
    test_lapSLPmat();
    A = [];
    return
end

rsrc = rsrc(:);
rtar = rtar(:);

R = abs(bsxfun(@minus,rtar,rsrc.'));
R(R==0) = 1; % avoid log(0) in pathological self-eval calls
A = -(1/(2*pi))*log(R);

end

function test_lapSLPmat
fprintf('lapSLPmat self-test\n');

ns = 120;
nt = 170;

t1 = linspace(0,2*pi,ns+1)';
t1 = t1(1:end-1);
rsrc = 0.7*(cos(t1)+1i*sin(t1));

t2 = linspace(0,2*pi,nt+1)';
t2 = t2(1:end-1);
rtar = 2 + cos(t2)+1i*sin(t2);

sigma = randn(ns,1);

A = lapSLPmat(rsrc,rtar);
u_mat = A*sigma;
u_dir = lapSLPfield(rsrc,rtar,sigma,false);

rel = norm(u_mat-u_dir,inf)/max(1,norm(u_dir,inf));
fprintf('  rel. err (A*sigma vs direct): %.3e\n',rel);

end
