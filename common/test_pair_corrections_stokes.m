function report = test_pair_corrections_stokes(opt, solve_mob_fun, solve_res_fun)
%TEST_PAIR_CORRECTIONS_STOKES Preflight check for Stokes pair-correction solvers.
%
% Usage:
%   report = test_pair_corrections_stokes(opt, solve_mob_fun, solve_res_fun)
%
% Inputs:
%   opt             - Options struct passed to the pair-correction solvers.
%   solve_mob_fun   - Mobility solver under test, e.g. @solve_mob_peanut_enhanced
%                     or @solve_mob_2B_enhanced.
%   solve_res_fun   - Resistance solver under test, e.g. @solve_res_peanut_enhanced
%                     or @solve_res_2B_enhanced.
%
% The 1B companion solvers are inferred from the supplied pair-correction
% solver names by replacing the suffix '_peanut', '_peanut_enhanced',
% '_2B', or '_2B_enhanced' with '_1B'.
%
% The test runs three checks:
%   1. Two bodies at delta = delta_pair with one-body preconditioning.
%   2. A representative two-body separation delta_star where getEnhancedGrid
%      finds no ellipse-segment enhancement nodes.
%   3. A single-body solve using N_c and Rp_c from opt.
%
% For each case, both directions are run:
%   mobility -> resistance and resistance -> mobility.
%
% Anna Broms, Mar 2026

if nargin < 3
    error('test_pair_corrections_stokes requires opt, solve_mob_fun, and solve_res_fun.');
end

mob_name = solverLabel(solve_mob_fun);
res_name = solverLabel(solve_res_fun);
mob_1b_fun = inferOneBodyName(mob_name);
res_1b_fun = inferOneBodyName(res_name);

R = getScalarRadius(opt, 1);
delta_pair = getOptField(opt,'delta_pair',0.2*R);
if ~isfinite(delta_pair) || delta_pair <= 0
    error('opt.delta_pair must be a positive finite scalar.');
end

opt_eval = prepareStokesOpt(opt, R);
Rp_f = getOptField(opt_eval,'Rp_f',[]);

% For equal radii, the accumulation point measured from the left center is
% zacc1 = (D - sqrt(D^2 - 4R^2))/2 with D = 2R + delta.
% Setting zacc1 = Rp_f gives delta_star = (R - Rp_f)^2 / Rp_f.
delta_star = (R - Rp_f)^2 / Rp_f;

q_pair = [0; 2*R + delta_pair];
q_star = [0; 2*R + delta_star];
q_single = 0;

F_pair = [1 0; -1 0];
T_pair = [1; -1];
U_pair = [1 0; -1 0];
W_pair = [1; -1];

F_single = [1 0];
T_single = 1;
U_single = [1 0];
W_single = 1;

fprintf('=== test_pair_corrections_stokes ===\n');
fprintf('mob solver  : %s\n', mob_name);
fprintf('res solver  : %s\n', res_name);
fprintf('1B mob     : %s\n', mob_1b_fun);
fprintf('1B res     : %s\n', res_1b_fun);
fprintf('R = %.3g, delta_pair = %.3g, delta_star = %.3g\n\n', R, delta_pair, delta_star);

report = struct();
report.radius = R;
report.delta_pair = delta_pair;
report.delta_star = delta_star;
report.solve_mob = mob_name;
report.solve_res = res_name;
report.solve_mob_1B = mob_1b_fun;
report.solve_res_1B = res_1b_fun;

printSection('1. Two bodies at delta = delta_pair, one-body preconditioning');
report.delta_pair_case = runStokesCase(q_pair, mob_1b_fun, res_1b_fun, ...
    F_pair, T_pair, U_pair, W_pair, opt_eval);
printStokesCase(report.delta_pair_case, '1B');

if isfinite(delta_star) && delta_star < delta_pair
    printSection('2. Two bodies at delta = delta_star, pair corrections');
    report.delta_star_case = runStokesCase(q_star, mob_name, res_name, ...
        F_pair, T_pair, U_pair, W_pair, opt_eval);
    printStokesCase(report.delta_star_case, 'pair');
else
    report.delta_star_case = [];
    printSection('2. Two bodies at delta = delta_star, pair corrections');
    fprintf('No representative no-ellipse-node separation was found below delta_pair.\n');
    fprintf('Skipping the pair-correction check.\n');
end

printSection('3. Single body, one-body preconditioning');
report.single_case = runStokesCase(q_single, mob_1b_fun, res_1b_fun, ...
    F_single, T_single, U_single, W_single, opt_eval);
printStokesCase(report.single_case, '1B');

end

function result = runStokesCase(q, mob_fun, res_fun, F_ref, T_ref, U_ref, W_ref, opt)
[UW_mob, sol_mob, coeff_mob] = callMobilitySolver(mob_fun, q, F_ref, T_ref, opt);
[FT_res, sol_res, coeff_res] = callResistanceSolver(res_fun, q, U_ref, W_ref, opt);

[F_back, T_back] = unpackFT(FT_res);
[U_back, W_back] = unpackUW(UW_mob);

[FT_back, ~, ~] = callResistanceSolver(res_fun, q, U_back, W_back, opt);
[UW_back, ~, ~] = callMobilitySolver(mob_fun, q, F_back, T_back, opt);

result = struct();
result.mobility = summarizeStokesSol(sol_mob, coeff_mob);
result.resistance = summarizeStokesSol(sol_res, coeff_res);
result.two_way.mob_to_res = relerr(FT_back, packFT(F_ref, T_ref));
result.two_way.res_to_mob = relerr(UW_back, packUW(U_ref, W_ref));
end

function [UW, sol, coeff] = callMobilitySolver(fun, q, F, T, opt)
name = solverLabel(fun);
if contains(name, '_1B')
    rad = getOptField(opt,'rad',1);
    if isempty(rad)
        rad = 1;
    end
    if isscalar(rad)
        rad = rad(1) * ones(numel(q),1);
    else
        rad = rad(:);
    end
    image = getOptField(opt,'image',1);
    lr = getOptField(opt,'lr',0);
    visualise = false;
    gmres_tol = getOptField(opt,'gmres_tol',1e-10);
    debug = false;
    surface_error_mode = char(getOptField(opt,'surface_error_mode','rel'));
    gmres_verbose = getOptField(opt,'gmres_verbose',0);
    [UW, coeff, it, ~, rel_res, abs_res] = feval(fun, q, F, T, rad, image, lr, ...
        visualise, gmres_tol, debug, surface_error_mode, gmres_verbose);
    sol = struct('it',it,'rel_res',rel_res,'abs_res',abs_res);
else
    [UW, sol] = feval(fun, q, F, T, opt);
    coeff = extractCoefficients(sol);
end
end

function [FT, sol, coeff] = callResistanceSolver(fun, q, U, W, opt)
name = solverLabel(fun);
if contains(name, '_1B')
    rad = getOptField(opt,'rad',1);
    if isempty(rad)
        rad = 1;
    end
    if isscalar(rad)
        rad = rad(1) * ones(numel(q),1);
    else
        rad = rad(:);
    end
    image = getOptField(opt,'image',1);
    lr = getOptField(opt,'lr',0);
    visualise = false;
    gmres_tol = getOptField(opt,'gmres_tol',1e-10);
    debug = false;
    gmres_verbose = getOptField(opt,'gmres_verbose',0);
    [FT, coeff, it, ~, maxres] = feval(fun, q, U, W, rad, image, lr, ...
        visualise, gmres_tol, debug, gmres_verbose);
    sol = struct('it',it,'maxres',maxres);
else
    [FT, sol] = feval(fun, q, U, W, opt);
    coeff = extractCoefficients(sol);
end
end

function s = summarizeStokesSol(sol, coeffs)
s = struct();
s.it = getFieldOr(sol,'it',NaN);
if isfield(sol,'rel_res') && ~isempty(sol.rel_res)
    s.residual = sol.rel_res;
elseif isfield(sol,'maxres') && ~isempty(sol.maxres)
    s.residual = sol.maxres;
else
    s.residual = NaN;
end
s.coeff_inf = max(abs(coeffs));
s.coeff_l2 = norm(coeffs);
end

function printStokesCase(result, solver_tag)
fprintf('  %-10s %12s %14s %14s\n','solver','it','residual','|coeff|_inf');
fprintf('  %-10s %12d %14.3e %14.3e\n',sprintf('%s mob',solver_tag), ...
    result.mobility.it, result.mobility.residual, result.mobility.coeff_inf);
fprintf('  %-10s %12d %14.3e %14.3e\n',sprintf('%s res',solver_tag), ...
    result.resistance.it, result.resistance.residual, result.resistance.coeff_inf);
fprintf('  %-10s %16s %16s\n','two-way','mob->res','res->mob');
fprintf('  %-10s %16.3e %16.3e\n',' ',result.two_way.mob_to_res,result.two_way.res_to_mob);
fprintf('\n');
end

function opt_out = prepareStokesOpt(opt, R)
opt_out = opt;
opt_out.rad = R;
opt_out.visualise_sol = false;
opt_out.visualise = 0;
opt_out.visualise_grid = false;
opt_out.debug = false;
opt_out.gmres_verbose = 0;
opt_out.get_bndry_field = true;
end

function R = getScalarRadius(opt, defaultR)
R = getOptField(opt,'rad',defaultR);
if isempty(R)
    R = defaultR;
end
R = R(1);
end

function name = solverLabel(fun)
if isa(fun,'function_handle')
    name = func2str(fun);
else
    name = char(fun);
end
end

function name_1B = inferOneBodyName(name)
name_1B = regexprep(name,'(_peanut(_enhanced)?|_2B(_enhanced)?)$','_1B');
if strcmp(name_1B, name)
    error('Could not infer a 1B solver name from "%s".', name);
end
end

function coeffs = extractCoefficients(sol)
if isfield(sol,'lambda_proxy') && ~isempty(sol.lambda_proxy)
    coeffs = sol.lambda_proxy;
elseif isfield(sol,'lambda_c') && ~isempty(sol.lambda_c)
    coeffs = sol.lambda_c;
elseif isfield(sol,'lambda_all') && ~isempty(sol.lambda_all)
    coeffs = sol.lambda_all;
elseif isfield(sol,'lambda') && ~isempty(sol.lambda)
    coeffs = sol.lambda;
else
    coeffs = [];
end
end

function uw = packUW(U,W)
P = size(U,1);
uw = zeros(3*P,1);
uw(1:3:end) = U(:,1);
uw(2:3:end) = U(:,2);
uw(3:3:end) = W;
end

function [U,W] = unpackUW(UW)
P = numel(UW)/3;
U = [UW(1:3:end), UW(2:3:end)];
W = UW(3:3:end);
assert(size(U,1)==P,'Unexpected size in unpackUW.');
end

function ft = packFT(F,T)
P = size(F,1);
ft = zeros(3*P,1);
ft(1:3:end) = F(:,1);
ft(2:3:end) = F(:,2);
ft(3:3:end) = T;
end

function [F,T] = unpackFT(FT)
P = numel(FT)/3;
F = [FT(1:3:end), FT(2:3:end)];
T = FT(3:3:end);
assert(size(F,1)==P,'Unexpected size in unpackFT.');
end

function v = relerr(a,b)
v = norm(a-b,inf)/max(1,norm(b,inf));
end

function val = getFieldOr(s, field, default_val)
if isfield(s, field) && ~isempty(s.(field))
    val = s.(field);
else
    val = default_val;
end
end

function printSection(title_str)
bar = repmat('=',1,78);
fprintf('%s\n',bar);
fprintf('%s\n',title_str);
fprintf('%s\n\n',bar);
end
