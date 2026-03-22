function report = test_pair_corrections_laplace(opt, solve_cap_fun, solve_elast_fun)
%TEST_PAIR_CORRECTIONS_LAPLACE Preflight check for Laplace pair-correction solvers.
%
% Usage:
%   report = test_pair_corrections_laplace(opt, solve_cap_fun, solve_elast_fun)
%
% Inputs:
%   opt              - Options struct passed to the pair-correction solvers.
%   solve_cap_fun    - Capacitance solver under test, e.g. @solve_cap_peanut
%                      or @solve_cap_2B.
%   solve_elast_fun  - Elastance solver under test, e.g. @solve_elast_peanut
%                      or @solve_elast_2B.
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
%   capacitance -> elastance and elastance -> capacitance.
%
% Anna Broms, Mar 2026

if nargin < 3
    error('test_pair_corrections_laplace requires opt, solve_cap_fun, and solve_elast_fun.');
end

cap_name = solverLabel(solve_cap_fun);
elast_name = solverLabel(solve_elast_fun);
cap_1b_fun = inferOneBodyName(cap_name);
elast_1b_fun = inferOneBodyName(elast_name);

R = opt.rad;
delta_pair = opt.delta_pair;

opt_eval = prepareLaplaceOpt(opt, R);
Rp_c = opt_eval.Rp_c;

% For equal radii, the accumulation point measured from the left center is
% zacc1 = (D - sqrt(D^2 - 4R^2))/2 with D = 2R + delta.
% Setting zacc1 = Rp_c gives delta_star = (R - Rp_c)^2 / Rp_c.
delta_star = (R - Rp_c)^2 / Rp_c;

q_pair = [0; 2*R + delta_pair];
q_star = [0; 2*R + delta_star];
q_single = 0;

v_pair = [1; -0.5];
Q_pair = [1; -1];
v_single = 1;
Q_single = 1;

fprintf('=== test_pair_corrections_laplace ===\n');
fprintf('cap solver  : %s\n', cap_name);
fprintf('elast solver : %s\n', elast_name);
fprintf('1B cap      : %s\n', cap_1b_fun);
fprintf('1B elast    : %s\n', elast_1b_fun);
fprintf('R = %.3g, delta_pair = %.3g, delta_star = %.3g\n\n', R, delta_pair, delta_star);

report = struct();
report.radius = R;
report.delta_pair = delta_pair;
report.delta_star = delta_star;
report.solve_cap = cap_name;
report.solve_elast = elast_name;
report.solve_cap_1B = cap_1b_fun;
report.solve_elast_1B = elast_1b_fun;

printSection('1. Two bodies at delta = delta_pair, one-body preconditioning');
report.delta_pair_case = runLaplaceCase(q_pair, cap_1b_fun, elast_1b_fun, v_pair, Q_pair, opt_eval);
printLaplaceCase(report.delta_pair_case, '1B');

if isfinite(delta_star) && delta_star < delta_pair
    printSection('2. Two bodies at delta = delta_star, pair corrections');
    report.delta_star_case = runLaplaceCase(q_star, cap_name, elast_name, v_pair, Q_pair, opt_eval);
    printLaplaceCase(report.delta_star_case, 'pair');
else
    report.delta_star_case = [];
    printSection('2. Two bodies at delta = delta_star, pair corrections');
    fprintf('No representative no-ellipse-node separation was found below delta_pair.\n');
    fprintf('Skipping the pair-correction check.\n');
end

printSection('3. Single body, one-body preconditioning');
report.single_case = runLaplaceCase(q_single, cap_1b_fun, elast_1b_fun, v_single, Q_single, opt_eval);
printLaplaceCase(report.single_case, '1B');

end

function result = runLaplaceCase(q, cap_fun, elast_fun, v_ref, Q_ref, opt)
[Q_cap, sol_cap] = feval(cap_fun, q, v_ref, opt);
[v_elast, sol_elast] = feval(elast_fun, q, Q_ref, opt);

[v_back, ~] = feval(elast_fun, q, Q_cap, opt);
[Q_back, ~] = feval(cap_fun, q, v_elast, opt);

result = struct();
result.cap = summarizeLaplaceSol(sol_cap, Q_cap);
result.elast = summarizeLaplaceSol(sol_elast, v_elast);
result.two_way.cap_to_elast = relerr(v_back, v_ref);
result.two_way.elast_to_cap = relerr(Q_back, Q_ref);
end

function s = summarizeLaplaceSol(sol, coeffs)
s = struct();
s.it = getFieldOr(sol,'it',NaN);
s.maxres = getFieldOr(sol,'maxres',NaN);
s.coeff_inf = max(abs(coeffs));
s.coeff_l2 = norm(coeffs);
end

function printLaplaceCase(result, solver_tag)
fprintf('  %-10s %12s %14s %14s\n','solver','it','maxres','|coeff|_inf');
fprintf('  %-10s %12d %14.3e %14.3e\n',sprintf('%s cap',solver_tag), ...
    result.cap.it, result.cap.maxres, result.cap.coeff_inf);
fprintf('  %-10s %12d %14.3e %14.3e\n',sprintf('%s elast',solver_tag), ...
    result.elast.it, result.elast.maxres, result.elast.coeff_inf);
fprintf('  %-10s %16s %16s\n','two-way','cap->elast','elast->cap');
fprintf('  %-10s %16.3e %16.3e\n',' ',result.two_way.cap_to_elast,result.two_way.elast_to_cap);
fprintf('\n');
end

function opt_out = prepareLaplaceOpt(opt, R)
opt_out = opt;
opt_out.rad = R;
opt_out.visualise_sol = false;
opt_out.visualise = 0;
opt_out.visualise_grid = false;
opt_out.debug = false;
opt_out.gmres_verbose = 0;
opt_out.get_bndry_field = true;
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
