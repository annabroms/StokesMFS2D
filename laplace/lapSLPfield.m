function u = lapSLPfield(rsrc,rtar,sigma,use_fmm)
%LAPSLPFIELD Evaluate scalar Laplace SLP field.
%
% Syntax:
%   u = lapSLPfield(rsrc,rtar,sigma,use_fmm)
%
% Inputs:
%   rsrc    - Complex valued array of source locations.
%   rtar    - Complex valued array of target locations.
%   sigma   - Source strengths, size numel(rsrc) x nd.
%   use_fmm - Logical flag. If true, use fmm2d (of flatiron) when available,
%             otherwise use the compiled direct evaluator r2ddir.
%
% Output:
%   u       - Potential values at targets, size numel(rtar) x nd.
%
% Self-test:
%   lapSLPfield()
%   Compares fmm2d output to compiled direct summation.
%
% See also: lapSLPdirect, lapSLPmat.
%
% Anna Broms, Mar 2026

if nargin==0
    test_lapSLPfield;
    return
end

if nargin < 4 || isempty(use_fmm)
    use_fmm = true;
end

rsrc = rsrc(:);
rtar = rtar(:);
ns = numel(rsrc);

if ns==0 || isempty(rtar)
    u = zeros(numel(rtar),size(sigma,2));
    return
end

assert(size(sigma,1)==ns,'sigma must have one row per source location.');

if ~use_fmm
    srcinfo = struct();
    srcinfo.sources = [real(rsrc)'; imag(rsrc)'];
    srcinfo.nd = size(sigma,2);
    srcinfo.charges = (-sigma.')/(2*pi);

    targ = [real(rtar)'; imag(rtar)'];
    U = r2ddir(srcinfo,targ,1);
    u = U.pottarg.';
    % u = lapSLPdirect(rsrc,rtar,sigma); % direct summation, for testing
    return
end


nd = size(sigma,2);
targ = [real(rtar)'; imag(rtar)'];

srcinfo = struct();
srcinfo.sources = [real(rsrc)'; imag(rsrc)'];
srcinfo.nd = nd;
srcinfo.charges = (-sigma.')/(2*pi); % fmm2d (of flatiron) uses log(r)

eps_fmm = 1e-10;
U = rfmm2d(eps_fmm,srcinfo,0,targ,1);

if nd==1
    u = U.pottarg(:);
else
    u = reshape(U.pottarg,nd,numel(rtar)).';
end

end

function test_lapSLPfield
fprintf('lapSLPfield self-test\n');

ns = 200;
nt = 220;

rsrc = rand(ns,1)+1i*rand(ns,1);
rtar = rand(nt,1)+1i*rand(nt,1);

sigma = randn(ns,1);

u_dir = lapSLPfield(rsrc,rtar,sigma,false);
u_fmm = lapSLPfield(rsrc,rtar,sigma,true);

rel = norm(u_fmm-u_dir,inf)/max(1,norm(u_dir,inf));
fprintf('  rel. err (fmm vs direct): %.3e\n',rel);

end
