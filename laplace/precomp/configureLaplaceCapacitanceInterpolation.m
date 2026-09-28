function [opt,mode] = configureLaplaceCapacitanceInterpolation(opt)
%CONFIGURELAPLACECAPACITANCEINTERPOLATION Validate the smooth model profile.
%
% This is deliberately limited to the 160-by-160 canonical map. Setting
% use_interpolation to 'full', 'reduced', or the preferred
% 'reduced_noconst' selects the smooth construction used in the documented
% experiments: Tikhonov regularisation, a frozen 150-node ellipse
% discretisation, and no separate TSVD compression of Cmap. Panels, orders,
% and reduced ranks are trained for the requested tolerances.

mode = getLaplaceInterpolationMode(opt);
opt.use_interpolation = mode;
if strcmp(mode,'none')
    return
end

R = getOptField(opt,'rad',[]);
if isempty(R) || ~isscalar(R) || ~isreal(R) || ~isfinite(R) || R <= 0
    error('configureLaplaceCapacitanceInterpolation:BadRadius', ...
        'Interpolation currently requires one positive scalar radius.');
end

if logical(getOptField(opt,'project_charge',false))
    error('configureLaplaceCapacitanceInterpolation:CapacitanceOnly', ...
        ['opt.use_interpolation is currently implemented only for the ', ...
         'Laplace capacitance problem (opt.project_charge=false).']);
end
% The interpolation model is canonical, so every physical pair must be
% rotated/refitted to and from the aligned reference frame.
opt.cmap = true;
opt.reuse_pair_basis_by_sep = true;
opt.compress_cmap = false;

% Smooth profile required by the interpolation trainer.  The ellipse
% discretisation remains fixed so the operator family has constant size.
opt.ellipse_constant = true;
opt.Nclust = 150;
opt.use_tikhonov = true;
opt.tikhonov_tol = 1e-11;
if ~isfield(opt,'smallest_delta') || isempty(opt.smallest_delta)
    opt.smallest_delta = 1e-3*R;
end

interpolation_tol = getOptField(opt,'interpolation_tol',1e-6);
charge_interpolation_tol = getOptField(opt, ...
    'charge_interpolation_tol',1e-8);
validateattributes(interpolation_tol,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'opt.interpolation_tol');
validateattributes(charge_interpolation_tol,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'opt.charge_interpolation_tol');
opt.interpolation_tol = interpolation_tol;
opt.charge_interpolation_tol = charge_interpolation_tol;

N_c = getOptField(opt,'N_c',80);
N_cmap = getOptField(opt,'N_cmap',N_c);
refit = logical(getOptField(opt,'refit',false));
validateattributes(N_cmap,{'numeric'},{'scalar','integer','positive'}, ...
    mfilename,'opt.N_cmap');
if ~refit && N_c ~= N_cmap
    error('configureLaplaceCapacitanceInterpolation:FourierGridMismatch', ...
        'refit=false requires N_c=N_cmap for Fourier rotation.');
end

delta_pair = getOptField(opt,'delta_pair',nan);
if ~isscalar(delta_pair) || ~isfinite(delta_pair) || ...
        delta_pair <= opt.smallest_delta
    error('configureLaplaceCapacitanceInterpolation:BadGapRange', ...
        'Require 0 < opt.smallest_delta < opt.delta_pair.');
end
end
