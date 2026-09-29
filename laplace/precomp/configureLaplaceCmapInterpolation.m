function [opt,mode] = configureLaplaceCmapInterpolation(opt,problem)
%CONFIGURELAPLACECMAPINTERPOLATION Validate a smooth Cmap model profile.

if nargin < 2
    problem = 'capacitance';
end
problem = resolveLaplaceInterpolationProblem(problem,mfilename);
mode = getLaplaceInterpolationMode(opt);
opt.use_interpolation = mode;
opt.interpolation_problem = problem;
opt.project_charge = strcmp(problem,'elastance');
if strcmp(mode,'none')
    return
end
if strcmp(problem,'elastance') && ~strcmp(mode,'full')
    error('configureLaplaceCmapInterpolation:ElastanceFullOnly', ...
        ['Laplace elastance supports only use_interpolation=''full''; ', ...
         'reduced interpolation modes are capacitance-only.']);
end

R = getOptField(opt,'rad',[]);
if isempty(R) || ~isscalar(R) || ~isreal(R) || ~isfinite(R) || R <= 0
    error('configureLaplaceCmapInterpolation:BadRadius', ...
        'Interpolation currently requires one positive scalar radius.');
end

opt.cmap = true;
opt.reuse_pair_basis_by_sep = true;
opt.compress_cmap = false;
opt.ellipse_constant = true;
Nclust = getOptField(opt,'Nclust',150);
if isempty(Nclust)
    Nclust = 150;
end
validateattributes(Nclust,{'numeric'}, ...
    {'scalar','integer','positive','finite'},mfilename, ...
    'opt.Nclust');
opt.Nclust = Nclust;
opt.use_tikhonov = true;
tikhonov_tol = getOptField(opt,'tikhonov_tol',[]);
if isempty(tikhonov_tol)
    tikhonov_tol = 1e-11;
end
validateattributes(tikhonov_tol,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'opt.tikhonov_tol');
opt.tikhonov_tol = tikhonov_tol;
if ~isfield(opt,'smallest_delta') || isempty(opt.smallest_delta)
    opt.smallest_delta = 1e-3*R;
end

interpolation_tol = getOptField(opt,'interpolation_tol',1e-6);
has_new = isfield(opt,'volt_charge_interp_tol') && ...
    ~isempty(opt.volt_charge_interp_tol);
has_legacy = isfield(opt,'charge_interpolation_tol') && ...
    ~isempty(opt.charge_interpolation_tol);
if has_new && has_legacy && ...
        ~isequaln(opt.volt_charge_interp_tol,opt.charge_interpolation_tol)
    error('configureLaplaceCmapInterpolation:ConflictingTolerance', ...
        ['opt.volt_charge_interp_tol and the legacy ', ...
         'opt.charge_interpolation_tol must agree when both are supplied.']);
elseif has_new
    volt_charge_interp_tol = opt.volt_charge_interp_tol;
elseif has_legacy
    volt_charge_interp_tol = opt.charge_interpolation_tol;
else
    volt_charge_interp_tol = 1e-8;
end
validateattributes(interpolation_tol,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'opt.interpolation_tol');
validateattributes(volt_charge_interp_tol,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename, ...
    'opt.volt_charge_interp_tol');
opt.interpolation_tol = interpolation_tol;
opt.volt_charge_interp_tol = volt_charge_interp_tol;

N_c = getOptField(opt,'N_c',80);
N_cmap = getOptField(opt,'N_cmap',N_c);
refit = logical(getOptField(opt,'refit',false));
validateattributes(N_cmap,{'numeric'},{'scalar','integer','positive'}, ...
    mfilename,'opt.N_cmap');
if ~refit && N_c ~= N_cmap
    error('configureLaplaceCmapInterpolation:FourierGridMismatch', ...
        'refit=false requires N_c=N_cmap for Fourier rotation.');
end

delta_pair = getOptField(opt,'delta_pair',nan);
if ~isscalar(delta_pair) || ~isfinite(delta_pair) || ...
        delta_pair <= opt.smallest_delta
    error('configureLaplaceCmapInterpolation:BadGapRange', ...
        'Require 0 < opt.smallest_delta < opt.delta_pair.');
end
end
