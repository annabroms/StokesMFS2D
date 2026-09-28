% sep23_laplace_offgrid_pair_rotation.m
%
% Diagnose off-grid Fourier rotation of a reused two-body Laplace basis.
% The input follows the production path: smooth data are sampled on each
% particle boundary and mapped to coarse proxy coefficients with the
% one-body pseudoinverse before the pair correction is applied.
%
% The experiment separates three questions:
%   1. Does rotating the one-body coefficients preserve their field?
%   2. Does rotating the compressed pair output preserve its field?
%   3. Does the complete reused pair correction agree with a pair basis
%      built directly at the physical orientation?
%
% The direct FFT shift and the zero-padded "oversampled_fft" shift are
% compared explicitly. Two least-squares controls are included: a reusable
% per-particle refit based on the one-body pseudoinverse, and a joint fit on
% the peanut surface that tests the full two-particle representation.
%
% Anna Broms, Sep 23, 2026

script_name = mfilename;
script_date = 'Sep 23, 2026';
repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

fprintf('=== %s (%s) ===\n',script_name,script_date);

% -------------------------------------------------------------------------
% Configuration. Values may be assigned before running this script.
% -------------------------------------------------------------------------
if ~exist('N_values','var') || isempty(N_values)
    N_values = [40 80 120];
end
if ~exist('shift_fractions','var') || isempty(shift_fractions)
    shift_fractions = 0:0.05:1;
end
if ~exist('R','var') || isempty(R)
    R = 2;
end
if ~exist('relative_gap','var') || isempty(relative_gap)
    relative_gap = 1e-3;
end
if ~exist('proxy_ratio','var') || isempty(proxy_ratio)
    % Hold the proxy geometry fixed at the N=80 default while N changes.
    proxy_ratio = max(1-log(1/1e-12)/80,0.01);
end
if ~exist('N_peanut','var') || isempty(N_peanut)
    N_peanut = 200;
end
if ~exist('rotation_oversample','var') || isempty(rotation_oversample)
    rotation_oversample = 8;
end
if ~exist('plotfig','var') || isempty(plotfig)
    plotfig = true;
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end
if ~exist('field_tolerance','var') || isempty(field_tolerance)
    field_tolerance = 1e-6;
end

N_values = N_values(:).';
shift_fractions = shift_fractions(:).';
data_names = {'fourier','harmonic'};
d = 2*R*(1+relative_gap);
assert(all(mod(N_values,2)==0),'N_values must be even for the Nyquist diagnostic.');
assert(mod(N_peanut,8)==0,'N_peanut must be divisible by 8.');

config = struct();
config.N_values = N_values;
config.shift_fractions = shift_fractions;
config.R = R;
config.relative_gap = relative_gap;
config.surface_gap = d-2*R;
config.proxy_ratio = proxy_ratio;
config.N_peanut = N_peanut;
config.rotation_oversample = rotation_oversample;
config.data_names = data_names;
config.field_tolerance = field_tolerance;

fprintf(['Pair: R=%.3g, surface gap/R=%.3e, fixed proxy radius/R=%.6f\n', ...
    'Offsets: %d fractions of one grid step, plus 60 degrees\n\n'], ...
    R,config.surface_gap/R,proxy_ratio,numel(shift_fractions));

% Verify the interpolation primitive independently of any pair basis.
interpolator = run_interpolator_checks(N_values,rotation_oversample);

records = repmat(empty_record(),0,1);

for N = N_values
    nout = ceil(1.2*N);
    nout_f = ceil(1.2*N);

    t_c = (0:N-1)'*(2*pi/N);
    t_out = (0:nout-1)'*(2*pi/nout);
    t_f = t_c;
    t_out_f = (0:nout_f-1)'*(2*pi/nout_f);

    rbase_in_c = R*proxy_ratio*exp(1i*t_c);
    rbase_out_c = R*exp(1i*t_out);
    rbase_in_f = R*proxy_ratio*exp(1i*t_f);
    rout_base_f = R*exp(1i*t_out_f);

    [Uself,Yself] = getSelfPseudoLaplace(1,rbase_in_c, ...
        rbase_out_c,[0 nout],false);

    angle_values = [shift_fractions*(2*pi/N) pi/3];
    angle_kinds = [repmat({'fractional'},size(shift_fractions)) ...
        {'sixty_degree'}];

    if verbose
        fprintf('N=%d: %d orientations x %d smooth inputs\n', ...
            N,numel(angle_values),numel(data_names));
    end

    for ia = 1:numel(angle_values)
        theta = angle_values(ia);
        rot = exp(1i*theta);
        q = [-0.5*d*rot; 0.5*d*rot];
        mid = mean(q);

        opt = make_options(N,R,proxy_ratio,N_peanut,rotation_oversample);
        [~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
        assert(size(pairs,1)==1,'Expected exactly one close pair.');

        % check_rotations builds the reused canonical group and a directly
        % discretized pair at this physical orientation in the same call.
        [~,~,~,~,~,~,pair_cache] = getPairBasisLaplace( ...
            q,rbase_in_c,rbase_in_f,rout_base_f,rbase_out_c, ...
            rimage_vec,refine,pairs,opt);
        canonical = pair_cache.groups(1);
        direct = pair_cache.check_pairs(1);

        opt_fft = opt;
        opt_fft.rotation_mode = 'fft';
        opt_over = opt;
        opt_over.rotation_mode = 'oversampled_fft';
        phase_fft = make_phases(N,rot,opt_fft);
        phase_over = make_phases(N,rot,opt_over);

        grids = make_physical_grids(q,mid,rot,R,N,rbase_in_c, ...
            rbase_in_f,canonical,d,N_peanut);
        joint_fit_factors = make_refit_factors(grids.coarse_actual, ...
            grids.peanut_train);
        particle_refit_map = make_particle_refit_map(rbase_in_c, ...
            rbase_out_c,rot,Yself{1},Uself{1});

        for idata = 1:numel(data_names)
            data_name = data_names{idata};
            tau_blocks = smooth_boundary_data(data_name,q,rbase_out_c, ...
                t_out,R);

            % This is exactly the one-body application in
            % transform_lap_peanut.
            lambda_blocks = Yself{1}*(Uself{1}*tau_blocks);
            lambda = lambda_blocks(:);

            direct_data = apply_direct_pair(direct,lambda,N,grids);
            fft_data = apply_reused_pair(canonical,lambda_blocks, ...
                phase_fft,N,grids,joint_fit_factors,particle_refit_map);
            over_data = apply_reused_pair(canonical,lambda_blocks, ...
                phase_over,N,grids,joint_fit_factors,particle_refit_map);

            rec = empty_record();
            rec.N = N;
            rec.data_name = data_name;
            rec.angle_kind = angle_kinds{ia};
            rec.theta = theta;
            rec.theta_degrees = theta*180/pi;
            rec.shift_fraction = mod(theta/(2*pi/N),1);
            if strcmp(angle_kinds{ia},'fractional') && ...
                    abs(shift_fractions(ia)-1) < 100*eps
                rec.shift_fraction = 1;
            end

            [rec.input_tail,rec.input_nyquist] = ...
                spectrum_metrics(lambda_blocks);
            [rec.output_tail_fft,rec.output_nyquist_fft] = ...
                spectrum_metrics(fft_data.tau_canonical);

            rec.input_transfer_boundary_fft = fft_data.input_boundary;
            rec.input_transfer_boundary_over = over_data.input_boundary;
            rec.coarse_coeff_error_fft = relative_error( ...
                fft_data.tau_actual,direct_data.tau);
            rec.coarse_coeff_error_over = relative_error( ...
                over_data.tau_actual,direct_data.tau);
            rec.fine_coeff_error_fft = relative_error( ...
                fft_data.beta_actual,direct_data.beta);
            rec.fine_coeff_error_over = relative_error( ...
                over_data.beta_actual,direct_data.beta);

            rec.output_transfer_peanut_fft = fft_data.output_peanut;
            rec.output_transfer_peanut_over = over_data.output_peanut;
            rec.output_refit_peanut_fft = fft_data.refit_peanut;
            rec.output_refit_peanut_over = over_data.refit_peanut;
            rec.output_particle_refit_peanut_fft = ...
                fft_data.particle_refit_peanut;
            rec.output_particle_refit_peanut_over = ...
                over_data.particle_refit_peanut;
            rec.output_transfer_far_fft = fft_data.output_far;
            rec.output_transfer_far_over = over_data.output_far;
            rec.output_refit_far_fft = fft_data.refit_far;
            rec.output_refit_far_over = over_data.refit_far;
            rec.output_particle_refit_far_fft = fft_data.particle_refit_far;
            rec.output_particle_refit_far_over = over_data.particle_refit_far;

            rec.u_corr_boundary_fft = scaled_error( ...
                fft_data.u_boundary,direct_data.u_boundary, ...
                direct_data.scale_boundary);
            rec.u_corr_boundary_over = scaled_error( ...
                over_data.u_boundary,direct_data.u_boundary, ...
                direct_data.scale_boundary);
            rec.u_corr_peanut_fft = scaled_error( ...
                fft_data.u_peanut,direct_data.u_peanut, ...
                direct_data.scale_peanut);
            rec.u_corr_peanut_over = scaled_error( ...
                over_data.u_peanut,direct_data.u_peanut, ...
                direct_data.scale_peanut);
            rec.u_corr_far_fft = scaled_error( ...
                fft_data.u_far,direct_data.u_far,direct_data.scale_far);
            rec.u_corr_far_over = scaled_error( ...
                over_data.u_far,direct_data.u_far,direct_data.scale_far);
            rec.u_corr_boundary_particle_refit_fft = scaled_error( ...
                fft_data.u_boundary_particle_refit,direct_data.u_boundary, ...
                direct_data.scale_boundary);
            rec.u_corr_boundary_particle_refit_over = scaled_error( ...
                over_data.u_boundary_particle_refit,direct_data.u_boundary, ...
                direct_data.scale_boundary);
            rec.u_corr_peanut_particle_refit_fft = scaled_error( ...
                fft_data.u_peanut_particle_refit,direct_data.u_peanut, ...
                direct_data.scale_peanut);
            rec.u_corr_peanut_particle_refit_over = scaled_error( ...
                over_data.u_peanut_particle_refit,direct_data.u_peanut, ...
                direct_data.scale_peanut);
            rec.u_corr_far_particle_refit_fft = scaled_error( ...
                fft_data.u_far_particle_refit,direct_data.u_far, ...
                direct_data.scale_far);
            rec.u_corr_far_particle_refit_over = scaled_error( ...
                over_data.u_far_particle_refit,direct_data.u_far, ...
                direct_data.scale_far);
            rec.u_corr_boundary_reference_norm = norm(direct_data.u_boundary);
            rec.u_corr_boundary_scale = direct_data.scale_boundary;
            rec.u_corr_peanut_reference_norm = norm(direct_data.u_peanut);
            rec.u_corr_peanut_scale = direct_data.scale_peanut;
            rec.u_corr_far_reference_norm = norm(direct_data.u_far);
            rec.u_corr_far_scale = direct_data.scale_far;

            rec.fft_over_coeff_difference = relative_error( ...
                fft_data.tau_actual,over_data.tau_actual);
            rec.fft_over_u_corr_difference = scaled_error( ...
                fft_data.u_peanut,over_data.u_peanut, ...
                direct_data.scale_peanut);

            records(end+1,1) = rec; %#ok<SAGROW>
        end
    end
end

cases = struct2table(records);
acceptance = assess_results(cases,interpolator,N_values,data_names);
diagnosis = classify_results(cases,N_values,data_names,field_tolerance);

result = struct();
result.config = config;
result.interpolator = interpolator;
result.cases = cases;
result.acceptance = acceptance;
result.diagnosis = diagnosis;
result.metric_notes = struct( ...
    'relative_error','norm(test-reference)/norm(reference)', ...
    'u_corr_error',[ ...
        'norm(reused-direct)/max(norm(direct fine field), ', ...
        'norm(direct compressed field), norm(direct correction))'], ...
    'particle_refit',[ ...
        'Each Cmap output block is mapped by A^+*B_theta; A^+ is the ', ...
        'reused one-body pseudoinverse and B_theta evaluates its ', ...
        'co-rotated canonical source circle on the fixed boundary grid'], ...
    'joint_refit','Both output blocks are fitted together on the peanut surface', ...
    'tail_energy','fraction of FFT energy in modes |k|>N/4');

print_summary(result);
if plotfig
    plot_results(result);
end

% =========================================================================
% Local helpers
% =========================================================================

function opt = make_options(N,R,proxy_ratio,N_peanut,rotation_oversample)
opt = getLaplace2Dparams(2,R,N,N);
opt.Rp_c = R*proxy_ratio;
opt.Rp_f = R*proxy_ratio;
opt.delta_pair = 0.2*R;
opt.N_peanut = N_peanut;
opt.cmap = true;
opt.compress_cmap = false;
opt.project_charge = false;
opt.reuse_pair_basis_by_sep = true;
opt.check_rotations = true;
opt.parallel_precomp = false;
opt.rotation_mode = 'fft';
opt.rotation_oversample = rotation_oversample;
opt.use_fmm = false;
opt.get_bndry_field = true;
opt.show_counter = false;
opt.visualise_grid = false;
end

function phases = make_phases(N,rot,opt)
phases.coarse = getUniformCircleRotationSpec(N,rot,opt);
phases.coarse_inv = invertUniformCircleRotationSpec(phases.coarse);
phases.fine = getUniformCircleRotationSpec(N,rot,opt);
phases.fine_inv = invertUniformCircleRotationSpec(phases.fine);
end

function tau = smooth_boundary_data(name,q,rbase_out,t,R)
nout = numel(t);
tau = zeros(nout,2);
switch name
    case 'fourier'
        tau(:,1) = cos(2*t) + 0.30*sin(3*t) + 0.10*cos(7*t);
        tau(:,2) = 0.80*sin(t-0.2) - 0.25*cos(4*t+0.1) + ...
            0.10*sin(6*t);
    case 'harmonic'
        z = [q(1)+rbase_out q(2)+rbase_out];
        pole = mean(q) + R*(4+1.3i);
        tau = -log(abs((z-pole)/R));
    otherwise
        error('Unknown data family "%s".',name);
end
end

function grids = make_physical_grids(q,mid,rot,R,N,rbase_c,rbase_f, ...
        canonical,d,N_peanut)
grids = struct();
grids.coarse_actual = [q(1)+rbase_c; q(2)+rbase_c];
grids.coarse_rotated = [q(1)+rot*rbase_c; q(2)+rot*rbase_c];

image_i = mid + rot*canonical.rimage_canon{1};
image_j = mid + rot*canonical.rimage_canon{2};
grids.nimage = [numel(image_i) numel(image_j)];
grids.fine_actual = [q(1)+rbase_f; image_i; ...
    q(2)+rbase_f; image_j];
grids.fine_rotated = [q(1)+rot*rbase_f; image_i; ...
    q(2)+rot*rbase_f; image_j];

n_boundary = 4*N;
t_boundary = (0:n_boundary-1)'*(2*pi/n_boundary) + ...
    0.37*(2*pi/n_boundary);
base_boundary = R*exp(1i*t_boundary);
grids.boundary = [q(1)+base_boundary; q(2)+base_boundary];

grids.peanut_train = createPeanut(q(1),q(2),2*N_peanut,0,R);
grids.peanut = createPeanut(q(1),q(2),3*N_peanut,0,R);

n_far = 4*N;
t_far = (0:n_far-1)'*(2*pi/n_far) + 0.29*(2*pi/n_far);
grids.far = mean(q) + 3*R*exp(1i*t_far);

% Keep the canonical separation in the result construction explicit.
grids.separation = d;
end

function factors = make_refit_factors(source_actual,targets_train)
A = lapSLPmat(source_actual,targets_train);
[Y,U] = getPseudoFactors(A,1e-14,0);
factors.Y = Y;
factors.U = U;
end

function transfer = make_particle_refit_map(rbase_in,rbase_out,rot,Yself,Uself)
% Re-express the field of a co-rotated source circle on the fixed source
% grid. Yself*Uself is the reusable one-body pseudoinverse A^+; only the
% inexpensive evaluation B_theta changes with the pair angle.
Btheta = lapSLPmat(rot*rbase_in,rbase_out);
transfer = Yself*(Uself*Btheta);
end

function data = apply_direct_pair(pair,lambda,N,grids)
data = struct();
data.tau = pair.Cmap*lambda;
data.beta = pair.Ypf*(pair.Upf*lambda);

q_pair = pair.q_pair(:);
rbase_c = grids.coarse_actual(1:N)-q_pair(1);
rbase_f = grids.fine_actual(1:N)-q_pair(1);
source_c = [q_pair(1)+rbase_c; q_pair(2)+rbase_c];
source_f = [q_pair(1)+rbase_f; pair.rimage_canon{1}; ...
    q_pair(2)+rbase_f; pair.rimage_canon{2}];

[data.u_boundary,data.scale_boundary] = correction_field( ...
    source_f,data.beta,source_c,data.tau,grids.boundary);
[data.u_peanut,data.scale_peanut] = correction_field( ...
    source_f,data.beta,source_c,data.tau,grids.peanut);
[data.u_far,data.scale_far] = correction_field( ...
    source_f,data.beta,source_c,data.tau,grids.far);
end

function data = apply_reused_pair(pair,lambda_blocks,phases,N,grids, ...
        joint_fit,particle_refit_map)
data = struct();

lambda_canonical = rotateUniformCircleData(lambda_blocks,[],phases.coarse);
rhs = lambda_canonical(:);

tau_canonical = reshape(pair.Cmap*rhs,N,2);
tau_actual = rotateUniformCircleData(tau_canonical,[],phases.coarse_inv);
tau_particle_refit = particle_refit_map*tau_canonical;
beta_canonical = pair.Ypf*(pair.Upf*rhs);
beta_actual = transfer_fine_coefficients(beta_canonical,N, ...
    grids.nimage,phases.fine_inv);

data.tau_canonical = tau_canonical;
data.tau_actual = tau_actual(:);
data.tau_particle_refit = tau_particle_refit(:);
data.beta_actual = beta_actual;

% Stage 1: does the one-body coefficient transfer preserve its field?
u_input_old = lapSLPmat(grids.coarse_rotated,grids.boundary)*rhs;
u_input_new = lapSLPmat(grids.coarse_actual,grids.boundary)* ...
    lambda_blocks(:);
data.input_boundary = relative_error(u_input_old,u_input_new);

% Stage 2: transfer the canonical compressed output to the fixed physical
% source grid, and compare with a fresh least-squares fit on that grid.
u_output_train = lapSLPmat(grids.coarse_rotated,grids.peanut_train)* ...
    tau_canonical(:);
tau_refit = joint_fit.Y*(joint_fit.U'*u_output_train);

u_output_old = lapSLPmat(grids.coarse_rotated,grids.peanut)* ...
    tau_canonical(:);
u_output_fft = lapSLPmat(grids.coarse_actual,grids.peanut)*tau_actual(:);
u_output_refit = lapSLPmat(grids.coarse_actual,grids.peanut)*tau_refit;
u_output_particle_refit = lapSLPmat(grids.coarse_actual,grids.peanut)* ...
    tau_particle_refit(:);
data.output_peanut = relative_error(u_output_fft,u_output_old);
data.refit_peanut = relative_error(u_output_refit,u_output_old);
data.particle_refit_peanut = relative_error( ...
    u_output_particle_refit,u_output_old);

u_output_old = lapSLPmat(grids.coarse_rotated,grids.far)* ...
    tau_canonical(:);
u_output_fft = lapSLPmat(grids.coarse_actual,grids.far)*tau_actual(:);
u_output_refit = lapSLPmat(grids.coarse_actual,grids.far)*tau_refit;
u_output_particle_refit = lapSLPmat(grids.coarse_actual,grids.far)* ...
    tau_particle_refit(:);
data.output_far = relative_error(u_output_fft,u_output_old);
data.refit_far = relative_error(u_output_refit,u_output_old);
data.particle_refit_far = relative_error(u_output_particle_refit,u_output_old);

% Stage 3: full fine-minus-compressed pair correction.
data.u_boundary = correction_field(grids.fine_actual,beta_actual, ...
    grids.coarse_actual,tau_actual(:),grids.boundary);
data.u_peanut = correction_field(grids.fine_actual,beta_actual, ...
    grids.coarse_actual,tau_actual(:),grids.peanut);
data.u_far = correction_field(grids.fine_actual,beta_actual, ...
    grids.coarse_actual,tau_actual(:),grids.far);
data.u_boundary_particle_refit = correction_field(grids.fine_actual, ...
    beta_actual,grids.coarse_actual,tau_particle_refit(:),grids.boundary);
data.u_peanut_particle_refit = correction_field(grids.fine_actual, ...
    beta_actual,grids.coarse_actual,tau_particle_refit(:),grids.peanut);
data.u_far_particle_refit = correction_field(grids.fine_actual, ...
    beta_actual,grids.coarse_actual,tau_particle_refit(:),grids.far);
end

function beta_actual = transfer_fine_coefficients(beta,N,nimage,phase_inv)
n1 = nimage(1);
n2 = nimage(2);
f1 = 1:N;
e1 = N+(1:n1);
f2 = N+n1+(1:N);
e2 = 2*N+n1+(1:n2);

fine = rotateUniformCircleData([beta(f1) beta(f2)],[],phase_inv);
beta_actual = beta;
beta_actual(f1) = fine(:,1);
beta_actual(f2) = fine(:,2);
% Image nodes rotate geometrically, so their scalar strengths are unchanged.
beta_actual(e1) = beta(e1);
beta_actual(e2) = beta(e2);
end

function [u,scale] = correction_field(source_f,beta,source_c,tau,targets)
u_fine = lapSLPmat(source_f,targets)*beta;
u_coarse = lapSLPmat(source_c,targets)*tau;
u = u_fine-u_coarse;
scale = max([norm(u_fine),norm(u_coarse),norm(u),realmin]);
end

function checks = run_interpolator_checks(N_values,oversample)
checks = repmat(struct('N',0,'bandlimited_fft',0, ...
    'bandlimited_oversampled',0,'fft_vs_oversampled',0, ...
    'integer_shift',0,'nyquist_fft',0,'nyquist_oversampled',0), ...
    numel(N_values),1);

for k = 1:numel(N_values)
    N = N_values(k);
    t = (0:N-1)'*(2*pi/N);
    theta = 0.37*(2*pi/N);
    values = cos(2*t) + 0.3*sin(3*t) + 0.1*cos(7*t);
    truth = cos(2*(t+theta)) + 0.3*sin(3*(t+theta)) + ...
        0.1*cos(7*(t+theta));

    opt_fft = struct('rotation_mode','fft', ...
        'rotation_oversample',oversample);
    opt_over = struct('rotation_mode','oversampled_fft', ...
        'rotation_oversample',oversample);
    fft_spec = getUniformCircleRotationSpec(N,exp(1i*theta),opt_fft);
    over_spec = getUniformCircleRotationSpec(N,exp(1i*theta),opt_over);
    vf = rotateUniformCircleData(values,[],fft_spec);
    vo = rotateUniformCircleData(values,[],over_spec);

    integer_spec = getUniformCircleRotationSpec(N,exp(1i*2*pi/N), ...
        opt_fft);
    vi = rotateUniformCircleData(values,[],integer_spec);

    nyquist = cos((N/2)*t);
    nyquist_truth = cos((N/2)*(t+theta));
    nyquist_fft = rotateUniformCircleData(nyquist,[],fft_spec);
    nyquist_over = rotateUniformCircleData(nyquist,[],over_spec);

    checks(k).N = N;
    checks(k).bandlimited_fft = relative_error(vf,truth);
    checks(k).bandlimited_oversampled = relative_error(vo,truth);
    checks(k).fft_vs_oversampled = relative_error(vf,vo);
    checks(k).integer_shift = relative_error(vi,circshift(values,-1));
    checks(k).nyquist_fft = relative_error(nyquist_fft,nyquist_truth);
    checks(k).nyquist_oversampled = relative_error( ...
        nyquist_over,nyquist_truth);
end

tol = 1e-11;
assert(max([checks.bandlimited_fft]) < tol, ...
    'Band-limited direct FFT rotation failed.');
assert(max([checks.bandlimited_oversampled]) < tol, ...
    'Band-limited oversampled FFT rotation failed.');
assert(max([checks.fft_vs_oversampled]) < tol, ...
    'Direct and oversampled FFT rotations differ.');
assert(max([checks.integer_shift]) < tol, ...
    'Integer grid rotation is not an exact permutation.');
end

function [tail,nyquist] = spectrum_metrics(values)
N = size(values,1);
F = fft(values,[],1);
energy = sum(abs(F).^2,1);
modes = [0:N/2 -N/2+1:-1]';
tail_energy = sum(abs(F(abs(modes)>N/4,:)).^2,1);
tail_each = tail_energy./max(energy,realmin);
tail = max(tail_each);

if mod(N,2)==0
    nyquist_each = abs(F(N/2+1,:)).^2./max(energy,realmin);
    nyquist = max(nyquist_each);
else
    nyquist = 0;
end
end

function err = relative_error(value,reference)
denom = norm(reference(:));
err = norm(value(:)-reference(:))/max(denom,realmin);
end

function err = scaled_error(value,reference,scale)
err = norm(value(:)-reference(:))/max(scale,realmin);
end

function rec = empty_record()
rec = struct( ...
    'N',0, ...
    'data_name','', ...
    'angle_kind','', ...
    'theta',0, ...
    'theta_degrees',0, ...
    'shift_fraction',0, ...
    'input_tail',0, ...
    'input_nyquist',0, ...
    'output_tail_fft',0, ...
    'output_nyquist_fft',0, ...
    'input_transfer_boundary_fft',0, ...
    'input_transfer_boundary_over',0, ...
    'coarse_coeff_error_fft',0, ...
    'coarse_coeff_error_over',0, ...
    'fine_coeff_error_fft',0, ...
    'fine_coeff_error_over',0, ...
    'output_transfer_peanut_fft',0, ...
    'output_transfer_peanut_over',0, ...
    'output_refit_peanut_fft',0, ...
    'output_refit_peanut_over',0, ...
    'output_particle_refit_peanut_fft',0, ...
    'output_particle_refit_peanut_over',0, ...
    'output_transfer_far_fft',0, ...
    'output_transfer_far_over',0, ...
    'output_refit_far_fft',0, ...
    'output_refit_far_over',0, ...
    'output_particle_refit_far_fft',0, ...
    'output_particle_refit_far_over',0, ...
    'u_corr_boundary_fft',0, ...
    'u_corr_boundary_over',0, ...
    'u_corr_peanut_fft',0, ...
    'u_corr_peanut_over',0, ...
    'u_corr_far_fft',0, ...
    'u_corr_far_over',0, ...
    'u_corr_boundary_particle_refit_fft',0, ...
    'u_corr_boundary_particle_refit_over',0, ...
    'u_corr_peanut_particle_refit_fft',0, ...
    'u_corr_peanut_particle_refit_over',0, ...
    'u_corr_far_particle_refit_fft',0, ...
    'u_corr_far_particle_refit_over',0, ...
    'u_corr_boundary_reference_norm',0, ...
    'u_corr_boundary_scale',0, ...
    'u_corr_peanut_reference_norm',0, ...
    'u_corr_peanut_scale',0, ...
    'u_corr_far_reference_norm',0, ...
    'u_corr_far_scale',0, ...
    'fft_over_coeff_difference',0, ...
    'fft_over_u_corr_difference',0);
end

function acceptance = assess_results(T,interpolator,N_values,data_names)
acceptance = struct();
acceptance.interpolator = max([interpolator.bandlimited_fft]) < 1e-11 && ...
    max([interpolator.bandlimited_oversampled]) < 1e-11 && ...
    max([interpolator.integer_shift]) < 1e-11;
acceptance.fft_matches_oversampled = ...
    max(T.fft_over_coeff_difference) < 1e-11 && ...
    max(T.fft_over_u_corr_difference) < 1e-11;

aligned = strcmp(T.angle_kind,'fractional') & ...
    (abs(T.shift_fraction) < 100*eps | abs(T.shift_fraction-1) < 100*eps);
acceptance.aligned_pair_rotation = ...
    max(T.u_corr_peanut_fft(aligned)) < 1e-8;

acceptance.spectral_convergence = true;
acceptance.half_step_errors = struct();
for idata = 1:numel(data_names)
    errors = nan(size(N_values));
    for k = 1:numel(N_values)
        rows = T.N==N_values(k) & strcmp(T.data_name,data_names{idata}) & ...
            strcmp(T.angle_kind,'fractional');
        fractions = T.shift_fraction(rows);
        values = T.u_corr_peanut_fft(rows);
        [~,idx] = min(abs(fractions-0.5));
        errors(k) = values(idx);
    end
    acceptance.half_step_errors.(data_names{idata}) = errors;
    if numel(errors)>1 && any(diff(errors) >= 0)
        acceptance.spectral_convergence = false;
    end
end
end

function diagnosis = classify_results(T,N_values,data_names,tol)
rows_out = repmat(struct('N',0,'data_name','', ...
    'max_input_transfer',0,'max_output_transfer',0, ...
    'max_particle_refit_error',0,'max_joint_refit_error',0, ...
    'max_u_corr_error',0,'max_u_corr_particle_refit_error',0, ...
    'classification',''),0,1);

for N = N_values
    for idata = 1:numel(data_names)
        name = data_names{idata};
        rows = T.N==N & strcmp(T.data_name,name) & ...
            strcmp(T.angle_kind,'fractional') & ...
            T.shift_fraction>0 & T.shift_fraction<1;

        item = struct();
        item.N = N;
        item.data_name = name;
        item.max_input_transfer = max(T.input_transfer_boundary_fft(rows));
        item.max_output_transfer = max(T.output_transfer_peanut_fft(rows));
        item.max_particle_refit_error = max( ...
            T.output_particle_refit_peanut_fft(rows));
        item.max_joint_refit_error = max(T.output_refit_peanut_fft(rows));
        item.max_u_corr_error = max(T.u_corr_peanut_fft(rows));
        item.max_u_corr_particle_refit_error = max( ...
            T.u_corr_peanut_particle_refit_fft(rows));

        if item.max_input_transfer > tol && item.max_u_corr_error > tol
            item.classification = 'input transfer unresolved';
        elseif item.max_output_transfer > tol && ...
                item.max_particle_refit_error <= tol
            item.classification = 'coefficient transfer limited; particle refit succeeds';
        elseif item.max_output_transfer > tol && ...
                item.max_joint_refit_error <= tol
            item.classification = 'coefficient transfer limited; joint refit succeeds';
        elseif item.max_output_transfer > tol && ...
                item.max_joint_refit_error > tol
            item.classification = 'shifted representation/resolution limited';
        elseif item.max_u_corr_error > tol && ...
                item.max_u_corr_particle_refit_error <= tol
            item.classification = 'Cmap output transfer limited; particle refit fixes correction';
        elseif item.max_u_corr_error > tol
            item.classification = 'pair map or fine-output transfer limited';
        else
            item.classification = 'resolved at requested tolerance';
        end
        rows_out(end+1,1) = item; %#ok<AGROW>
    end
end
diagnosis = struct2table(rows_out);
end

function print_summary(result)
T = result.cases;
fprintf('\nInterpolator checks\n');
fprintf('%6s %13s %13s %13s %13s\n', ...
    'N','FFT truth','Over truth','FFT vs over','integer');
for k = 1:numel(result.interpolator)
    c = result.interpolator(k);
    fprintf('%6d %13.3e %13.3e %13.3e %13.3e\n', ...
        c.N,c.bandlimited_fft,c.bandlimited_oversampled, ...
        c.fft_vs_oversampled,c.integer_shift);
end

fprintf('\nWorst off-grid errors over the one-grid-step sweep\n');
fprintf(['%6s %-10s %11s %11s %11s %11s %11s %11s ', ...
    '%11s %11s\n'], ...
    'N','data','input fld','coeff','output FFT','particle fit', ...
    'joint fit','u FFT','u pfit','FFT-vs-over');
for N = result.config.N_values
    for idata = 1:numel(result.config.data_names)
        name = result.config.data_names{idata};
        rows = T.N==N & strcmp(T.data_name,name) & ...
            strcmp(T.angle_kind,'fractional') & ...
            T.shift_fraction>0 & T.shift_fraction<1;
        fprintf(['%6d %-10s %11.3e %11.3e %11.3e %11.3e ', ...
            '%11.3e %11.3e %11.3e %11.3e\n'], ...
            N,name,max(T.input_transfer_boundary_fft(rows)), ...
            max(T.coarse_coeff_error_fft(rows)), ...
            max(T.output_transfer_peanut_fft(rows)), ...
            max(T.output_particle_refit_peanut_fft(rows)), ...
            max(T.output_refit_peanut_fft(rows)), ...
            max(T.u_corr_peanut_fft(rows)), ...
            max(T.u_corr_peanut_particle_refit_fft(rows)), ...
            max(T.fft_over_u_corr_difference(rows)));
    end
end

fprintf('\n60-degree comparison\n');
fprintf('%6s %-10s %10s %12s %12s %12s\n', ...
    'N','data','grid frac','input fld','output fld','u_corr');
rows60 = strcmp(T.angle_kind,'sixty_degree');
for row = find(rows60).'
    fprintf('%6d %-10s %10.3f %12.3e %12.3e %12.3e\n', ...
        T.N(row),T.data_name{row},T.shift_fraction(row), ...
        T.input_transfer_boundary_fft(row), ...
        T.output_transfer_peanut_fft(row), ...
        T.u_corr_peanut_fft(row));
end

fprintf('\nCoefficient spectral content (maximum over orientations)\n');
fprintf('%6s %-10s %12s %12s %12s %12s\n', ...
    'N','data','input tail','input Nyq','output tail','output Nyq');
for N = result.config.N_values
    for idata = 1:numel(result.config.data_names)
        name = result.config.data_names{idata};
        rows = T.N==N & strcmp(T.data_name,name);
        fprintf('%6d %-10s %12.3e %12.3e %12.3e %12.3e\n', ...
            N,name,max(T.input_tail(rows)),max(T.input_nyquist(rows)), ...
            max(T.output_tail_fft(rows)),max(T.output_nyquist_fft(rows)));
    end
end

fprintf('\nAcceptance summary\n');
fprintf('  interpolation primitives: %s\n', ...
    pass_fail(result.acceptance.interpolator));
fprintf('  fft equals oversampled_fft: %s\n', ...
    pass_fail(result.acceptance.fft_matches_oversampled));
fprintf('  aligned pair rotations: %s\n', ...
    pass_fail(result.acceptance.aligned_pair_rotation));
fprintf('  decreasing half-step error with N: %s\n', ...
    pass_fail(result.acceptance.spectral_convergence));
for idata = 1:numel(result.config.data_names)
    name = result.config.data_names{idata};
    values = result.acceptance.half_step_errors.(name);
    fprintf('    %-10s:',name);
    fprintf(' %.3e',values);
    fprintf('\n');
end


fprintf('\nField-based diagnosis at tolerance %.1e\n', ...
    result.config.field_tolerance);
for row = 1:height(result.diagnosis)
    fprintf('  N=%-4d %-10s %s\n',result.diagnosis.N(row), ...
        result.diagnosis.data_name{row}, ...
        result.diagnosis.classification{row});
end
end

function plot_results(result)
T = result.cases;
nN = numel(result.config.N_values);
nData = numel(result.config.data_names);

figure('Name','Sep 23 off-grid compressed-output transfer');
tiledlayout(nData,nN,'TileSpacing','compact','Padding','compact');
for idata = 1:nData
    name = result.config.data_names{idata};
    for k = 1:nN
        N = result.config.N_values(k);
        rows = T.N==N & strcmp(T.data_name,name) & ...
            strcmp(T.angle_kind,'fractional');
        [fraction,order] = sort(T.shift_fraction(rows));
        fft_err = T.output_transfer_peanut_fft(rows);
        over_err = T.output_transfer_peanut_over(rows);
        particle_refit_err = T.output_particle_refit_peanut_fft(rows);
        refit_err = T.output_refit_peanut_fft(rows);

        nexttile;
        semilogy(fraction,fft_err(order),'o-','DisplayName','FFT');
        hold on;
        semilogy(fraction,over_err(order),'x--', ...
            'DisplayName','oversampled FFT');
        semilogy(fraction,particle_refit_err(order),'d-', ...
            'DisplayName','particle refit');
        semilogy(fraction,refit_err(order),'s-', ...
            'DisplayName','joint peanut refit');
        grid on;
        xlabel('fraction of one grid step');
        ylabel('peanut field error');
        title(sprintf('%s, N=%d',name,N),'Interpreter','none');
        if idata==1 && k==1
            legend('Location','best');
        end
    end
end

figure('Name','Sep 23 off-grid full correction');
tiledlayout(nData,nN,'TileSpacing','compact','Padding','compact');
for idata = 1:nData
    name = result.config.data_names{idata};
    for k = 1:nN
        N = result.config.N_values(k);
        rows = T.N==N & strcmp(T.data_name,name) & ...
            strcmp(T.angle_kind,'fractional');
        [fraction,order] = sort(T.shift_fraction(rows));
        fft_err = T.u_corr_peanut_fft(rows);
        over_err = T.u_corr_peanut_over(rows);
        particle_refit_err = T.u_corr_peanut_particle_refit_fft(rows);

        nexttile;
        semilogy(fraction,fft_err(order),'o-','DisplayName','FFT');
        hold on;
        semilogy(fraction,over_err(order),'x--', ...
            'DisplayName','oversampled FFT');
        semilogy(fraction,particle_refit_err(order),'d-', ...
            'DisplayName','particle refit');
        grid on;
        xlabel('fraction of one grid step');
        ylabel('component-scaled correction error');
        title(sprintf('%s, N=%d',name,N),'Interpreter','none');
        if idata==1 && k==1
            legend('Location','best');
        end
    end
end
end

function value = pass_fail(tf)
if tf
    value = 'PASS';
else
    value = 'FAIL';
end
end
