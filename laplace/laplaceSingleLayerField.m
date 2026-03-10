function u = laplaceSingleLayerField(rsrc,rtar,sigma,use_fmm)
%LAPLACESINGLELAYERFIELD Evaluate scalar Laplace SLP field.
%
% Syntax:
%   u = laplaceSingleLayerField(rsrc,rtar,sigma,use_fmm)
%
% Inputs:
%   rsrc    - Complex source locations.
%   rtar    - Complex target locations.
%   sigma   - Source strengths, size numel(rsrc) x nd.
%   use_fmm - Logical flag. If true, use rfmm2d when available.
%
% Output:
%   u       - Potential values at targets, size numel(rtar) x nd.
%
% Self-test:
%   laplaceSingleLayerField()
%   Compares rfmm2d output to loop-based direct summation.
%
% Anna Broms style adaptation, Mar 2026

if nargin==0
    test_laplaceSingleLayerField;
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
    u = laplaceSLPDirect(rsrc,rtar,sigma);
    return
end

if exist('rfmm2d','file')~=2
    persistent warned_missing_rfmm2d
    if isempty(warned_missing_rfmm2d) || ~warned_missing_rfmm2d
        warning('laplaceSingleLayerField:MissingRfmm2d', ...
            'rfmm2d not found. Falling back to direct Laplace SLP evaluation.');
        warned_missing_rfmm2d = true;
    end
    u = laplaceSLPDirect(rsrc,rtar,sigma);
    return
end


nd = size(sigma,2);
targ = [real(rtar)'; imag(rtar)'];

srcinfo = struct();
srcinfo.sources = [real(rsrc)'; imag(rsrc)'];
srcinfo.nd = nd;
srcinfo.charges = (-sigma.')/(2*pi); % rfmm2d uses log(r)

eps_fmm = 1e-10;
U = rfmm2d(eps_fmm,srcinfo,0,targ,1);

if nd==1
    u = U.pottarg(:);
else
    u = reshape(U.pottarg,nd,numel(rtar)).';
end

end

function test_laplaceSingleLayerField
fprintf('laplaceSingleLayerField self-test\n');

ns = 200;
nt = 220;

rsrc = rand(ns,1)+1i*rand(ns,1);
rtar = rand(nt,1)+1i*rand(nt,1);

sigma = randn(ns,1);

u_dir = laplaceSingleLayerField(rsrc,rtar,sigma,false);
u_fmm = laplaceSingleLayerField(rsrc,rtar,sigma,true);

rel = norm(u_fmm-u_dir,inf)/max(1,norm(u_dir,inf));
fprintf('  rel. err (fmm vs direct): %.3e\n',rel);

end
