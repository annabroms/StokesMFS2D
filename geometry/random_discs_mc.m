function [q, meta] = random_discs_mc(P, opt)
%RANDOM_DISCS_MC Generate equal-radius random disc centers by MC sweeps.
%
% Syntax:
%   [q, meta] = random_discs_mc(P, opt)
%
% Inputs:
%   P   - Number of discs.
%   opt - Options struct with fields:
%         rad               disc radius, default 1
%         domain            'periodic' or 'boxed', default 'periodic'
%         phi               target area fraction (optional if L is given)
%         L                 square side length (optional if phi is given)
%         min_gap           minimum surface gap, default 1e-3*rad
%         n_sweeps          number of Monte Carlo sweeps, default 10
%         target_acceptance target accepted-move fraction, default 0.5
%         move_scale0       initial proposal radius, optional
%         rng_seed          random seed, optional
%         visualise         plot final configuration, default false
%
% Outputs:
%   q    - Complex-valued disc centers (P x 1).
%   meta - Diagnostics and setup metadata.
%
% Notes:
%   The proposal distribution is uniform in a disk of radius move_scale.
%   When move_scale0 is omitted, it is derived from the initial ordered
%   layout via move_scale0 = max(0.05*rad, 0.5*(d_init_min-d_min_allowed)).
%
%   Calling RANDOM_DISCS_MC with no inputs runs a self-test for both
%   periodic and boxed domains.
%
%   Apr 2026

if nargin == 0
    random_discs_mc_selftest();
    q = [];
    meta = struct();
    return
end

if nargin < 2 || isempty(opt)
    opt = struct();
elseif ~isstruct(opt)
    error('random_discs_mc:badOpt','opt must be a struct.');
end

validateattributes(P,{'numeric'},{'scalar','integer','positive','finite'}, ...
    mfilename,'P',1);

rad = get_opt(opt,'rad',1);
validateattributes(rad,{'numeric'},{'scalar','positive','finite'}, ...
    mfilename,'opt.rad');

domain = lower(char(string(get_opt(opt,'domain','periodic'))));
if ~ismember(domain,{'periodic','boxed'})
    error('random_discs_mc:badDomain', ...
        'opt.domain must be ''periodic'' or ''boxed''.');
end

min_gap_default = max(1e-3*rad,1e-12);
min_gap = get_opt(opt,'min_gap',min_gap_default);
validateattributes(min_gap,{'numeric'},{'scalar','positive','finite'}, ...
    mfilename,'opt.min_gap');

n_sweeps = get_opt(opt,'n_sweeps',10);
validateattributes(n_sweeps,{'numeric'},{'scalar','integer','nonnegative','finite'}, ...
    mfilename,'opt.n_sweeps');

target_acceptance = get_opt(opt,'target_acceptance',0.5);
validateattributes(target_acceptance,{'numeric'},{'scalar','>',0,'<',1}, ...
    mfilename,'opt.target_acceptance');

visualise = as_logical(get_opt(opt,'visualise',false),'opt.visualise');

if isfield(opt,'rng_seed') && ~isempty(opt.rng_seed)
    rng_seed = opt.rng_seed;
    validateattributes(rng_seed,{'numeric'},{'scalar','integer','finite'}, ...
        mfilename,'opt.rng_seed');
    rng(rng_seed,'twister');
else
    rng_seed = [];
end

[L,phi_target,phi] = resolve_side_length(P,rad,opt);
d_min_allowed = 2*rad + min_gap;

rad_eff = rad + 0.5*min_gap;
phi_eff = P*pi*rad_eff^2 / L^2;
phi_max_hex = pi/(2*sqrt(3));
if phi_eff > phi_max_hex*(1+1e-12)
    error('random_discs_mc:packingTooDense', ...
        ['Requested P, rad, L, and min_gap exceed the hexagonal packing ', ...
         'upper bound for the effective radius rad + min_gap/2.']);
end

if strcmp(domain,'boxed') && L < 2*rad
    error('random_discs_mc:boxedTooSmall', ...
        'Boxed geometry needs L >= 2*rad to fit even one disc.');
end

[q, layout_meta] = build_initial_configuration(P,L,rad,domain,d_min_allowed);
d_init_min = layout_meta.d_init_min;

if isfield(opt,'move_scale0') && ~isempty(opt.move_scale0)
    move_scale = opt.move_scale0;
    validateattributes(move_scale,{'numeric'},{'scalar','positive','finite'}, ...
        mfilename,'opt.move_scale0');
else
    slack0 = max(d_init_min - d_min_allowed,0);
    move_scale = max(0.05*rad,0.5*slack0);
end
move_scale0 = move_scale;

[bins, bin_state] = initialise_cell_list(q,L,d_min_allowed,domain);

acceptance_by_sweep = zeros(n_sweeps,1);
move_scale_by_sweep = zeros(n_sweeps,1);

for sweep = 1:n_sweeps
    accepted = 0;
    move_scale_by_sweep(sweep) = move_scale;

    for k = 1:P
        z_trial = q(k) + sample_uniform_disk(move_scale);
        if strcmp(domain,'periodic')
            z_trial = wrap_complex(z_trial,L);
        else
            if ~inside_box(z_trial,L,rad)
                continue
            end
        end

        if can_accept_move(k,z_trial,q,bins,bin_state,d_min_allowed,L,domain)
            old_bin = get_bin_index(q(k),bin_state);
            new_bin = get_bin_index(z_trial,bin_state);
            q(k) = z_trial;
            if any(old_bin ~= new_bin)
                bins = remove_from_bin(bins,old_bin,k);
                bins = add_to_bin(bins,new_bin,k);
            end
            accepted = accepted + 1;
        end
    end

    acceptance_by_sweep(sweep) = accepted / P;

    if acceptance_by_sweep(sweep) > target_acceptance + 0.05
        move_scale = 1.1*move_scale;
    elseif acceptance_by_sweep(sweep) < target_acceptance - 0.05
        move_scale = 0.9*move_scale;
    end
end

min_center_distance = compute_min_center_distance(q,domain,L);
min_surface_gap = min_center_distance - 2*rad;

meta = struct();
meta.L = L;
meta.phi = phi;
meta.phi_target = phi_target;
meta.phi_error = phi - phi_target;
if abs(phi_target) > 0
    meta.phi_rel_error = (phi - phi_target) / phi_target;
else
    meta.phi_rel_error = 0;
end
meta.rad = rad;
meta.domain = domain;
meta.min_gap = min_gap;
meta.d_min_allowed = d_min_allowed;
meta.n_sweeps = n_sweeps;
meta.acceptance_by_sweep = acceptance_by_sweep;
meta.move_scale_by_sweep = move_scale_by_sweep;
meta.min_center_distance = min_center_distance;
meta.min_surface_gap = min_surface_gap;
meta.rng_seed = rng_seed;
meta.target_acceptance = target_acceptance;
meta.move_scale0 = move_scale0;
meta.final_move_scale = move_scale;
if n_sweeps > 0
    meta.final_acceptance = acceptance_by_sweep(end);
else
    meta.final_acceptance = nan;
end
meta.initial_layout = layout_meta.layout;
meta.initial_layout_details = layout_meta.details;
meta.d_init_min = d_init_min;

if visualise
    plot_configuration(q,meta,sprintf('%s (%s)',mfilename,domain));
end

end

function [L,phi_target,phi] = resolve_side_length(P,rad,opt)
has_phi = isfield(opt,'phi') && ~isempty(opt.phi);
has_L = isfield(opt,'L') && ~isempty(opt.L);

if ~has_phi && ~has_L
    error('random_discs_mc:missingSize', ...
        'Provide either opt.phi or opt.L.');
end

if has_phi
    phi_target = opt.phi;
    validateattributes(phi_target,{'numeric'},{'scalar','positive','finite','<',1}, ...
        mfilename,'opt.phi');
else
    phi_target = [];
end

if has_L
    L = opt.L;
    validateattributes(L,{'numeric'},{'scalar','positive','finite'}, ...
        mfilename,'opt.L');
else
    L = [];
end

if has_phi && ~has_L
    L = sqrt(P*pi*rad^2/phi_target);
    phi = phi_target;
elseif has_L && ~has_phi
    phi = P*pi*rad^2 / L^2;
    phi_target = phi;
else
    phi = P*pi*rad^2 / L^2;
    rel_mismatch = abs(phi - phi_target) / max(abs(phi_target),eps);
    if rel_mismatch > 1e-10
        error('random_discs_mc:inconsistentPhiL', ...
            'opt.phi and opt.L are inconsistent for the given P and rad.');
    end
end
end

function [q_best, meta] = build_initial_configuration(P,L,rad,domain,d_min_allowed)
square = best_square_layout(P,L,rad,domain);
staggered = best_staggered_layout(P,L,rad,domain);

tol = 1e-12*max(d_min_allowed,1);
valid_square = square.d_init_min + tol >= d_min_allowed;
valid_staggered = staggered.d_init_min + tol >= d_min_allowed;

if ~valid_square && ~valid_staggered
    error('random_discs_mc:noValidInitialLayout', ...
        ['Could not place %d discs with the requested L, rad, and ', ...
         'min_gap using either the square grid or staggered-row layout.'],P);
end

if valid_staggered && (~valid_square || staggered.d_init_min > square.d_init_min)
    q_best = staggered.q;
    meta = rmfield(staggered,'q');
else
    q_best = square.q;
    meta = rmfield(square,'q');
end
end

function layout = best_square_layout(P,L,rad,domain)
best_d = -inf;
best_rows = 1;
best_cols = P;

if strcmp(domain,'boxed')
    width = L - 2*rad;
else
    width = L;
end

for nrows = 1:P
    ncols = ceil(P/nrows);
    dmin = square_layout_min_distance(width,nrows,ncols,domain);
    if dmin > best_d
        best_d = dmin;
        best_rows = nrows;
        best_cols = ncols;
    end
end

q = build_square_layout(P,L,rad,domain,best_rows,best_cols);
layout = struct();
layout.q = q;
layout.layout = 'square_grid';
layout.details = struct('nrows',best_rows,'ncols',best_cols, ...
    'search_score',best_d);
layout.d_init_min = compute_min_center_distance(q,domain,L);
end

function dmin = square_layout_min_distance(width,nrows,ncols,domain)
if nrows == 1 && ncols == 1
    dmin = inf;
    return
end

vals = [];
if strcmp(domain,'periodic')
    if ncols > 1
        vals(end+1) = width / ncols; %#ok<AGROW>
    end
    if nrows > 1
        vals(end+1) = width / nrows; %#ok<AGROW>
    end
else
    if ncols > 1
        vals(end+1) = width / (ncols-1); %#ok<AGROW>
    end
    if nrows > 1
        vals(end+1) = width / (nrows-1); %#ok<AGROW>
    end
end

if isempty(vals)
    dmin = inf;
else
    dmin = min(vals);
end
end

function q = build_square_layout(P,L,rad,domain,nrows,ncols)
if strcmp(domain,'periodic')
    sx = L / ncols;
    sy = L / nrows;
    x = -L/2 + (0.5 + (0:ncols-1))*sx;
    y = -L/2 + (0.5 + (0:nrows-1))*sy;
else
    width = L - 2*rad;
    x = centered_line(ncols,width);
    y = centered_line(nrows,width);
end

q = zeros(P,1);
idx = 1;
for row = 1:nrows
    for col = 1:ncols
        q(idx) = x(col) + 1i*y(row);
        idx = idx + 1;
        if idx > P
            return
        end
    end
end
end

function layout = best_staggered_layout(P,L,rad,domain)
best_d = -inf;
best_rows = 1;
best_cols = P;

if strcmp(domain,'boxed')
    width = L - 2*rad;
    for nrows = 1:P
        [ncols_long, dmin] = staggered_boxed_signature(P,width,nrows);
        if dmin > best_d
            best_d = dmin;
            best_rows = nrows;
            best_cols = ncols_long;
        end
    end
    q = build_staggered_boxed(P,L,rad,best_rows,best_cols);
    details = struct('nrows',best_rows,'ncols_long',best_cols, ...
        'ncols_short',max(best_cols-1,0), ...
        'search_score',best_d);
else
    width = L;
    for nrows = 1:P
        ncols = ceil(P/nrows);
        dmin = staggered_periodic_min_distance(width,nrows,ncols);
        if dmin > best_d
            best_d = dmin;
            best_rows = nrows;
            best_cols = ncols;
        end
    end
    q = build_staggered_periodic(P,L,best_rows,best_cols);
    details = struct('nrows',best_rows,'ncols',best_cols, ...
        'search_score',best_d);
end

layout = struct();
layout.q = q;
layout.layout = 'staggered_rows';
layout.details = details;
layout.d_init_min = compute_min_center_distance(q,domain,L);
end

function [ncols_long, dmin] = staggered_boxed_signature(P,width,nrows)
nlong = ceil(nrows/2);
nshort = floor(nrows/2);
ncols_long = 1;
while nlong*ncols_long + nshort*max(ncols_long-1,0) < P
    ncols_long = ncols_long + 1;
end

if nrows == 1 && ncols_long == 1
    dmin = inf;
elseif nrows == 1
    dmin = width / (ncols_long-1);
elseif ncols_long == 1
    if nlong <= 1
        dmin = inf;
    else
        dmin = 2*width / (nrows-1);
    end
else
    sx = width / (ncols_long-1);
    sy = width / (nrows-1);
    vals = [sx, hypot(0.5*sx,sy)];
    if nrows > 2
        vals(end+1) = 2*sy; %#ok<AGROW>
    end
    dmin = min(vals);
end
end

function dmin = staggered_periodic_min_distance(L,nrows,ncols)
if nrows == 1 && ncols == 1
    dmin = inf;
elseif nrows == 1
    dmin = L / ncols;
elseif ncols == 1
    sy = L / nrows;
    vals = [hypot(0.5*L,sy)];
    if nrows > 2
        vals(end+1) = 2*sy; %#ok<AGROW>
    end
    dmin = min(vals);
else
    sx = L / ncols;
    sy = L / nrows;
    vals = [sx, hypot(0.5*sx,sy)];
    if nrows > 2
        vals(end+1) = 2*sy; %#ok<AGROW>
    end
    dmin = min(vals);
end
end

function q = build_staggered_periodic(P,L,nrows,ncols)
sx = L / ncols;
sy = L / nrows;
x_base = -L/2 + (0.5 + (0:ncols-1))*sx;
y = -L/2 + (0.5 + (0:nrows-1))*sy;

q = zeros(P,1);
idx = 1;
for row = 1:nrows
    shift = 0.5*sx*mod(row-1,2);
    x = wrap_real(x_base + shift,L);
    x = sort(x);
    for col = 1:numel(x)
        q(idx) = x(col) + 1i*y(row);
        idx = idx + 1;
        if idx > P
            return
        end
    end
end
end

function q = build_staggered_boxed(P,L,rad,nrows,ncols_long)
width = L - 2*rad;
y = centered_line(nrows,width);

if ncols_long == 1
    x_long = 0;
    x_short = zeros(1,0);
else
    x_long = centered_line(ncols_long,width);
    sx = x_long(2) - x_long(1);
    x_short = x_long(1:end-1) + 0.5*sx;
end

q = zeros(P,1);
idx = 1;
for row = 1:nrows
    if mod(row,2) == 1
        x_row = x_long;
    else
        x_row = x_short;
    end
    for col = 1:numel(x_row)
        q(idx) = x_row(col) + 1i*y(row);
        idx = idx + 1;
        if idx > P
            return
        end
    end
end
end

function line = centered_line(npts,width)
if npts <= 0
    line = zeros(1,0);
elseif npts == 1
    line = 0;
else
    line = linspace(-width/2,width/2,npts);
end
end

function [bins, state] = initialise_cell_list(q,L,d_min_allowed,domain)
P = numel(q);
max_bins = max(1,ceil(sqrt(P)));
n_bins = max(1,min(max_bins,floor(L/d_min_allowed)));
if n_bins < 1
    n_bins = 1;
end

state = struct();
state.n_bins = n_bins;
state.bin_width = L / n_bins;
state.L = L;
state.domain = domain;

bins = cell(n_bins,n_bins);
for k = 1:P
    ij = get_bin_index(q(k),state);
    bins = add_to_bin(bins,ij,k);
end
end

function bins = add_to_bin(bins,ij,k)
bins{ij(1),ij(2)} = [bins{ij(1),ij(2)}, k];
end

function bins = remove_from_bin(bins,ij,k)
entries = bins{ij(1),ij(2)};
entries(entries == k) = [];
bins{ij(1),ij(2)} = entries;
end

function ij = get_bin_index(z,state)
x = real(z);
y = imag(z);
if strcmp(state.domain,'periodic')
    x = wrap_real(x,state.L);
    y = wrap_real(y,state.L);
end

ix = floor((x + state.L/2) / state.bin_width) + 1;
iy = floor((y + state.L/2) / state.bin_width) + 1;
ix = min(max(ix,1),state.n_bins);
iy = min(max(iy,1),state.n_bins);
ij = [ix, iy];
end

function tf = can_accept_move(k,z_trial,q,bins,state,d_min_allowed,L,domain)
ij = get_bin_index(z_trial,state);
nb_x = neighbor_bins(ij(1),state.n_bins,strcmp(domain,'periodic'));
nb_y = neighbor_bins(ij(2),state.n_bins,strcmp(domain,'periodic'));
d2_min = d_min_allowed^2;
tol = 1e-14*max(d2_min,1);

for ax = 1:numel(nb_x)
    for ay = 1:numel(nb_y)
        occupants = bins{nb_x(ax),nb_y(ay)};
        for m = 1:numel(occupants)
            j = occupants(m);
            if j == k
                continue
            end
            if pair_distance_sq(z_trial,q(j),domain,L) < d2_min - tol
                tf = false;
                return
            end
        end
    end
end

tf = true;
end

function bins_out = neighbor_bins(i,n_bins,is_periodic)
if n_bins == 1
    bins_out = 1;
    return
end

candidates = i-1:i+1;
if is_periodic
    bins_out = mod(candidates-1,n_bins) + 1;
else
    bins_out = candidates(candidates >= 1 & candidates <= n_bins);
end
bins_out = unique(bins_out,'stable');
end

function d2 = pair_distance_sq(z1,z2,domain,L)
dx = real(z1) - real(z2);
dy = imag(z1) - imag(z2);
if strcmp(domain,'periodic')
    dx = wrap_delta(dx,L);
    dy = wrap_delta(dy,L);
end
d2 = dx.^2 + dy.^2;
end

function delta = sample_uniform_disk(radius)
rho = radius * sqrt(rand(1));
theta = 2*pi*rand(1);
delta = rho * exp(1i*theta);
end

function tf = inside_box(z,L,rad)
half_width = L/2 - rad;
tf = abs(real(z)) <= half_width && abs(imag(z)) <= half_width;
end

function z = wrap_complex(z,L)
z = wrap_real(real(z),L) + 1i*wrap_real(imag(z),L);
end

function x = wrap_real(x,L)
x = mod(x + L/2, L) - L/2;
end

function dx = wrap_delta(dx,L)
dx = dx - L*round(dx / L);
end

function dmin = compute_min_center_distance(q,domain,L)
P = numel(q);
if P < 2
    dmin = inf;
    return
end

pts = [real(q(:)), imag(q(:))];
dmin = inf;
for i = 1:P-1
    dx = pts(i+1:end,1) - pts(i,1);
    dy = pts(i+1:end,2) - pts(i,2);
    if strcmp(domain,'periodic')
        dx = wrap_delta(dx,L);
        dy = wrap_delta(dy,L);
    end
    dmin = min(dmin,min(hypot(dx,dy)));
end
end

function plot_configuration(q,meta,fig_name)
t = linspace(0,2*pi,240).';
figure('Name',fig_name,'Color','w');
hold on
for k = 1:numel(q)
    zk = q(k) + meta.rad*(cos(t) + 1i*sin(t));
    plot(real(zk),imag(zk),'k-','LineWidth',0.9);
    plot(real(q(k)),imag(q(k)),'k.','MarkerSize',8);
end

half_L = meta.L/2;
plot([-half_L, half_L, half_L, -half_L, -half_L], ...
     [-half_L, -half_L, half_L, half_L, -half_L], ...
     'r--','LineWidth',1.2);
axis equal
grid on
xlabel('x')
ylabel('y')
title(sprintf('%s: P=%d, %s, min gap %.3g', ...
    fig_name,numel(q),meta.domain,meta.min_gap),'Interpreter','none')
end

function value = get_opt(opt,name,default_value)
if isfield(opt,name) && ~isempty(opt.(name))
    value = opt.(name);
else
    value = default_value;
end
end

function value = as_logical(value,name)
if ~(isscalar(value) && (islogical(value) || isnumeric(value)))
    error('random_discs_mc:badLogical','%s must be a logical scalar.',name);
end
value = logical(value);
end

function random_discs_mc_selftest()
delta = 0.001;
phi = 0.45;
P = 200;
cases = {
    struct('label','periodic', ...
           'P',P, ...
           'opt',struct('domain','periodic','phi',phi,'rad',1, ...
                        'min_gap',delta,'n_sweeps',12,'rng_seed',1, ...
                        'visualise',true)), ...
    struct('label','boxed', ...
           'P',P, ...
           'opt',struct('domain','boxed','phi',phi,'rad',1, ...
                        'min_gap',delta,'n_sweeps',12,'rng_seed',2, ...
                        'visualise',true)) ...
    };

tol = 1e-10;
fprintf('== %s self-test ==\n',mfilename);
for idx = 1:numel(cases)
    test_case = cases{idx};
    [q,meta] = random_discs_mc(test_case.P,test_case.opt); %#ok<NASGU>
    pass = meta.min_center_distance + tol >= meta.d_min_allowed && ...
        meta.min_surface_gap + tol >= meta.min_gap;
    if pass
        verdict = 'PASS';
    else
        verdict = 'FAIL';
    end

    fprintf('\nCase: %s\n',test_case.label);
    fprintf('  P = %d, rad = %.6g, domain = %s\n', ...
        test_case.P,meta.rad,meta.domain);
    fprintf('  target area fraction = %.12g\n',meta.phi_target);
    fprintf('  achieved area fraction = %.12g\n',meta.phi);
    fprintf('  area fraction error = %.12g\n',meta.phi_error);
    fprintf('  L = %.6f\n',meta.L);
    fprintf('  min allowed center distance = %.12g\n',meta.d_min_allowed);
    fprintf('  measured min center distance = %.12g\n',meta.min_center_distance);
    fprintf('  minimum allowed surface gap = %.12g\n',meta.min_gap);
    fprintf('  measured min surface gap = %.12g\n',meta.min_surface_gap);
    fprintf('  difference to threshold = %.12g\n', ...
        meta.min_center_distance - meta.d_min_allowed);
    fprintf('  final acceptance rate = %.6f\n',meta.final_acceptance);
    fprintf('  final move scale = %.12g\n',meta.final_move_scale);
    fprintf('  initial layout = %s\n',meta.initial_layout);
    fprintf('  verdict = %s\n',verdict);

    if isfield(test_case.opt,'visualise') && test_case.opt.visualise
        plot_configuration(q,meta,sprintf('%s self-test (%s)',mfilename,test_case.label));
    end
end
end
