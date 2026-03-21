function [udirect,vdirect] = stokSLPdirect(xsrc,ysrc,xtar,ytar,f1,f2,N,use_compiled)
%STOKSLPDIRECT  Evaluate 2D Stokeslet flow from discrete sources.
%
%   [udirect, vdirect] = STOKSLPDIRECT(xsrc, ysrc, xtar, ytar, f1, f2, N)
%   [udirect, vdirect] = STOKSLPDIRECT(..., use_compiled)
%
%   Computes the velocity field (u,v) at a set of target points due to N
%   Stokeslet sources in 2D, located at (xsrc, ysrc) with force strengths (f1, f2).
%   The velocity is evaluated using the free-space 2D Stokeslet kernel:
%
%       u(x) = 1/(4π) ∑_k [ -½ log(r_k²) I + (r_k ⊗ r_k)/r_k² ] * f_k,
%
%   where r_k = x - x_k is the vector from source to target, and f_k is the force.
%
%   INPUTS:
%     xsrc    - N×1 vector of x-coordinates of sources.
%     ysrc    - N×1 vector of y-coordinates of sources.
%     xtar    - M×1 vector of x-coordinates of targets.
%     ytar    - M×1 vector of y-coordinates of targets.
%     f1      - N×nd array of x-components of source strengths.
%     f2      - N×nd array of y-components of source strengths.
%     N       - Number of source points (should match size of xsrc, ysrc, f1, f2).
%     use_compiled - Logical flag. If true and fmm2d's compiled direct
%                    evaluator is available, use it; otherwise use the
%                    MATLAB loop below.
%
%   OUTPUTS:
%     udirect - M×nd array of x-components of velocity at each target point.
%     vdirect - M×nd array of y-components of velocity at each target point.
%
%   NOTES:
%     - Sources and targets must not coincide (no self-interaction).
%     - The kernel is the fundamental solution to the 2D Stokes equations
%       in free space with viscosity μ = 1.
%     - Multiple force vectors (nd > 1) are supported (batched columns).
%     - The compiled branch calls fmm2d's `st2ddir` wrapper.
%
% Anna Broms, Oct 19, 2025

if nargin==0
    test_stokSLPdirect;
    return
elseif nargin<8
    use_compiled = 0;
end


xsrc = xsrc(:);
ysrc = ysrc(:);
xtar = xtar(:);
ytar = ytar(:);

assert(numel(xsrc)==N,'Wrong number of x sources');
assert(numel(ysrc)==N,'Wrong number of y sources');
assert(size(f1,1)==N,'Wrong number of x source strengths');
assert(size(f2,1)==N,'Wrong number of y source strengths');

if use_compiled
    [udirect,vdirect] = stokSLPdirect_compiled(xsrc,ysrc,xtar,ytar,f1,f2,N);
else
    [udirect,vdirect] = stokSLPdirect_matlab(xsrc,ysrc,xtar,ytar,f1,f2,N);
end

end

function [udirect,vdirect] = stokSLPdirect_matlab(xsrc,ysrc,xtar,ytar,f1,f2,N)
nd = size(f1,2);

udirect = zeros(length(ytar),nd);
vdirect = zeros(length(ytar),nd);

for k = 1:N % loop over sources
    rx = xtar - xsrc(k);
    ry = ytar - ysrc(k);
    rho2 = rx.^2 + ry.^2;
    rdotf = rx.*f1(k,:) + ry.*f2(k,:);
    logrho2 = log(rho2);
    udirect = udirect-0.5*logrho2.*f1(k,:) + rdotf./rho2.*rx;
    vdirect = vdirect-0.5*logrho2.*f2(k,:) + rdotf./rho2.*ry;
end

vdirect = vdirect/4/pi;
udirect = udirect/4/pi;
end

function [udirect,vdirect] = stokSLPdirect_compiled(xsrc,ysrc,xtar,ytar,f1,f2,N)
nd = size(f1,2);

srcinfo = struct();
srcinfo.sources = [xsrc.'; ysrc.'];
srcinfo.nd = nd;
srcinfo.stoklet = permute(cat(3,f1.',f2.'),[1 3 2]);

U = st2ddir(srcinfo,[xtar.'; ytar.'],1);
pottarg = U.pottarg./(2*pi);

if nd == 1
    udirect = pottarg(1,:).';
    vdirect = pottarg(2,:).';
else
    udirect = squeeze(pottarg(:,1,:)).';
    vdirect = squeeze(pottarg(:,2,:)).';
end
end

function test_stokSLPdirect
fprintf('stokSLPdirect self-test and benchmark\n');

if exist('st2ddir','file') == 0 && exist('startup','file') == 2
    startup;
end
if exist('st2ddir','file') == 0
    error('stokSLPdirect:st2ddirMissing', ...
        'st2ddir is not on the MATLAB path. Run startup.m and make sure ../fmm2d is available.');
end

rng(21);
ns = 100;
nt = 120;
nrepeat = 20000;

xsrc = rand(ns,1);
ysrc = rand(ns,1);
xtar = 1.5 + rand(nt,1);
ytar = 1.5 + rand(nt,1);
f1 = randn(ns,1);
f2 = randn(ns,1);

[u_matlab,v_matlab] = stokSLPdirect(xsrc,ysrc,xtar,ytar,f1,f2,ns,false);
[u_compiled,v_compiled] = stokSLPdirect(xsrc,ysrc,xtar,ytar,f1,f2,ns,true);

field_matlab = [u_matlab; v_matlab];
field_compiled = [u_compiled; v_compiled];
abs_err = norm(field_compiled-field_matlab,inf);
rel_err = abs_err/max(1,norm(field_matlab,inf));

fprintf('sources=%d, targets=%d, repeats=%d\n',ns,nt,nrepeat);
fprintf('warm-up check: abs err %.3e, rel err %.3e\n',abs_err,rel_err);

tic;
for k = 1:nrepeat
    [u_matlab,v_matlab] = stokSLPdirect(xsrc,ysrc,xtar,ytar,f1,f2,ns,false);
end
t_matlab = toc;

tic;
for k = 1:nrepeat
    [u_compiled,v_compiled] = stokSLPdirect(xsrc,ysrc,xtar,ytar,f1,f2,ns,true);
end
t_compiled = toc;

fprintf('stokSLPdirect (MATLAB) total time: %.6f s\n',t_matlab);
fprintf('stokSLPdirect (compiled) total time: %.6f s\n',t_compiled);
fprintf('stokSLPdirect (MATLAB) time/call: %.3e s\n',t_matlab/nrepeat);
fprintf('stokSLPdirect (compiled) time/call: %.3e s\n',t_compiled/nrepeat);
fprintf('speedup (MATLAB / compiled): %.2fx\n',t_matlab/max(t_compiled,eps));
end
