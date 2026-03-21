%TEST_MOB_RES Compare Stokes mobility and resistance solvers.
%
% Enhanced solves in this script use Stokeslet sources only:
%   solve_mob_2B_enhanced, solve_mob_peanut_enhanced,
%   solve_res_2B_enhanced, solve_res_peanut_enhanced.
%
% Optional 1-body solves (enabled below) use the image representation with
% multiple source types:
%   solve_mob_1B, solve_res_1B  (from mobility/dev_image_based and
%   resistance/dev_image_based).
%
% If you want 2B/peanut image-based mixed-source solves instead of
% enhanced Stokeslet-only solves, switch calls to:
%   solve_mob_2B_images, solve_mob_peanut_images,
%   solve_res_2B_images, solve_res_peanut_images.
%
% Anna Broms, Mar 2026

clear;
close all;
clc;
rng(8);

%% Configuration
P = 20;
delta = 1e-3;                 % geometric spacing parameter
geometry_mode = 'cluster';       % 'line' | 'dumbbells' | 'cluster'
delta_pair = 0.2;             % near-pair threshold for 2B/peanut
N_peanut = 400;               % points per peanut separation surface


run_image_1B = false;          % include solve_mob_1B / solve_res_1B
image = 1;
lr = 0;

visualise = 0;
gmres_tol = 1e-8;
debug = 0;
gmres_verbose = 0;            % 0=silent, 1=summary, 2=per-iteration
surface_error_mode = 'rel';

q = buildGeometry(geometry_mode,P,delta);
P = numel(q);
rad = ones(P,1);

fprintf('=== test_mob_res ===\n');
fprintf('geometry=%s, P=%d, delta=%.2e, delta_pair=%.2e, N_peanut=%d\n', ...
    geometry_mode,P,delta,delta_pair,N_peanut);
fprintf('enhanced family: Stokeslet sources only\n');
if run_image_1B
    fprintf('1B family      : image representation with mixed source types\n');
end
fprintf('\n');

%% Geometry preview
printDivider('Geometry preview');
plotGeometryPreview(q,rad,sprintf('test_mob_res geometry (%s, P=%d)',geometry_mode,P));

printPauseDivider();
disp('Press key to continue...')
pause();

%% Initialize settings
opt_mob = get2Dparams(P);
opt_mob.delta_pair = delta_pair;
opt_mob.N_peanut = N_peanut;
opt_mob.visualise_sol = visualise;
opt_mob.gmres_tol = gmres_tol;
opt_mob.debug = debug;
opt_mob.surface_error_mode = surface_error_mode;
opt_mob.gmres_verbose = gmres_verbose;
opt_mob.rad = rad(1);
opt_mob.visualise_grid = 1; 
opt_mob.reuse_pair_basis_by_sep = 1;

%% Mobility comparison (input F,T, compare UW)
printDivider('Mobility comparison (input F,T -> output U,W)');
F_ref = randn(P,2);
T_ref = randn(P,1);

[UW_2B,sol_mob_2B] = solve_mob_2B_enhanced(q,F_ref,T_ref,opt_mob);
[UW_peanut,sol_mob_peanut] = solve_mob_peanut_enhanced(q,F_ref,T_ref,opt_mob);
gmres_iter_mob_2B = sol_mob_2B.it;
gmres_iter_mob_peanut = sol_mob_peanut.it;
boundary_rel_mob_2B = sol_mob_2B.rel_res;
boundary_rel_mob_peanut = sol_mob_peanut.rel_res;
boundary_abs_mob_2B = sol_mob_2B.abs_res;
boundary_abs_mob_peanut = sol_mob_peanut.abs_res;

if run_image_1B
    [UW_1B,~,gmres_iter_mob_1B,~,boundary_rel_mob_1B,boundary_abs_mob_1B] = ...
        solve_mob_1B(q,F_ref,T_ref,rad,image,lr,visualise,gmres_tol,debug,surface_error_mode,gmres_verbose);
end

printSummaryHeader('Mobility comparison');
fprintf('Mobility solves (input F,T)\n');
fprintf('  %-10s %12s %14s %14s\n','solver','gmres_iter','rel_res','abs_res');
if run_image_1B
    fprintf('  %-10s %12d %14.3e %14.3e\n','1B-image',gmres_iter_mob_1B,boundary_rel_mob_1B,boundary_abs_mob_1B);
end
fprintf('  %-10s %12d %14.3e %14.3e\n','2B-enh',gmres_iter_mob_2B,boundary_rel_mob_2B,boundary_abs_mob_2B);
fprintf('  %-10s %12d %14.3e %14.3e\n','peanut',gmres_iter_mob_peanut,boundary_rel_mob_peanut,boundary_abs_mob_peanut);

if run_image_1B
    mobility_diff_2B_vs_1B = relerr(UW_2B,UW_1B);
    mobility_diff_peanut_vs_1B = relerr(UW_peanut,UW_1B);
    fprintf('  diff vs 1B-image : ||UW_2B-UW_1B||/||UW_1B|| = %.3e\n',mobility_diff_2B_vs_1B);
    fprintf('                    ||UW_peanut-UW_1B||/||UW_1B|| = %.3e\n',mobility_diff_peanut_vs_1B);
end
mobility_diff_peanut_vs_2B = relerr(UW_peanut,UW_2B);
fprintf('  diff vs 2B-enh   : ||UW_peanut-UW_2B||/||UW_2B|| = %.3e\n\n',mobility_diff_peanut_vs_2B);

printPauseDivider();
disp('Press key to continue...')
pause();

%% Resistance comparison (input U,W, compare FT)
printDivider('Resistance comparison (input U,W -> output F,T)');
U_ref = randn(P,2);
W_ref = randn(P,1);

% opt_res = get2Dparams(P);
% opt_res.rad = rad;
% opt_res.delta_pair = delta_pair;
% opt_res.N_peanut = N_peanut;
% opt_res.lr = lr;
% opt_res.visualise = visualise;
% opt_res.gmres_tol = gmres_tol;
% opt_res.debug = debug;
% opt_res.gmres_verbose = gmres_verbose;
opt_res = opt_mob;

[FT_2B,sol_res_2B] = solve_res_2B_enhanced(q,U_ref,W_ref,opt_res);
[FT_peanut,sol_res_peanut] = solve_res_peanut_enhanced(q,U_ref,W_ref,opt_res);
gmres_iter_res_2B = sol_res_2B.it;
gmres_iter_res_peanut = sol_res_peanut.it;
boundary_rel_res_2B = sol_res_2B.rel_res;
boundary_rel_res_peanut = sol_res_peanut.rel_res;

if run_image_1B
    [FT_1B,~,gmres_iter_res_1B,~,boundary_rel_res_1B] = ...
        solve_res_1B(q,U_ref,W_ref,rad,image,lr,visualise,gmres_tol,debug,gmres_verbose);
end

printSummaryHeader('Resistance comparison');
fprintf('Resistance solves (input U,W)\n');
fprintf('  %-10s %12s %14s\n','solver','gmres_iter','rel_res');
if run_image_1B
    fprintf('  %-10s %12d %14.3e\n','1B-image',gmres_iter_res_1B,boundary_rel_res_1B);
end
fprintf('  %-10s %12d %14.3e\n','2B-enh',gmres_iter_res_2B,boundary_rel_res_2B);
fprintf('  %-10s %12d %14.3e\n','peanut',gmres_iter_res_peanut,boundary_rel_res_peanut);

if run_image_1B
    resistance_diff_2B_vs_1B = relerr(FT_2B,FT_1B);
    resistance_diff_peanut_vs_1B = relerr(FT_peanut,FT_1B);
    fprintf('  diff vs 1B-image : ||FT_2B-FT_1B||/||FT_1B|| = %.3e\n',resistance_diff_2B_vs_1B);
    fprintf('                    ||FT_peanut-FT_1B||/||FT_1B|| = %.3e\n',resistance_diff_peanut_vs_1B);
end
resistance_diff_peanut_vs_2B = relerr(FT_peanut,FT_2B);
fprintf('  diff vs 2B-enh   : ||FT_peanut-FT_2B||/||FT_2B|| = %.3e\n\n',resistance_diff_peanut_vs_2B);

printPauseDivider();
disp('Press key to continue...')
pause();

%% Two-way checks
printDivider('Two-way checks');
FT_ref = packFT(F_ref,T_ref);
UW_ref = packUW(U_ref,W_ref);

% Mobility -> Resistance
[U_m2,W_m2] = unpackUW(UW_2B);
[FT_back_m2,~] = solve_res_2B_enhanced(q,U_m2,W_m2,opt_res);

[U_mp,W_mp] = unpackUW(UW_peanut);
[FT_back_mp,~] = solve_res_peanut_enhanced(q,U_mp,W_mp,opt_res);

two_way_mob_res_2B = relerr(FT_back_m2,FT_ref);
two_way_mob_res_peanut = relerr(FT_back_mp,FT_ref);

% Resistance -> Mobility
[F_r2,T_r2] = unpackFT(FT_2B);
[UW_back_r2,~] = solve_mob_2B_enhanced(q,F_r2,T_r2,opt_mob);

[F_rp,T_rp] = unpackFT(FT_peanut);
[UW_back_rp,~] = solve_mob_peanut_enhanced(q,F_rp,T_rp,opt_mob);

two_way_res_mob_2B = relerr(UW_back_r2,UW_ref);
two_way_res_mob_peanut = relerr(UW_back_rp,UW_ref);

if run_image_1B
    [U_m1,W_m1] = unpackUW(UW_1B);
    [FT_back_m1,~,~,~,~] = solve_res_1B(q,U_m1,W_m1,rad,image,lr,visualise,gmres_tol,debug,gmres_verbose);
    [F_r1,T_r1] = unpackFT(FT_1B);
    [UW_back_r1,~,~,~,~,~] = solve_mob_1B(q,F_r1,T_r1,rad,image,lr,visualise,gmres_tol,debug,surface_error_mode,gmres_verbose);

    two_way_mob_res_1B = relerr(FT_back_m1,FT_ref);
    two_way_res_mob_1B = relerr(UW_back_r1,UW_ref);
end

printSummaryHeader('Two-way checks');
fprintf('  %-10s %16s %16s\n','solver','mob->res','res->mob');
if run_image_1B
    fprintf('  %-10s %16.3e %16.3e\n','1B-image',two_way_mob_res_1B,two_way_res_mob_1B);
end
fprintf('  %-10s %16.3e %16.3e\n','2B-enh',two_way_mob_res_2B,two_way_res_mob_2B);
fprintf('  %-10s %16.3e %16.3e\n','peanut',two_way_mob_res_peanut,two_way_res_mob_peanut);


function q = buildGeometry(mode,P,delta)
switch lower(mode)
    case 'line'
        spacing = 2 + delta;
        q = (0:P-1)'*spacing;
        q = complex(q,zeros(size(q)));
    case 'dumbbells'
        if mod(P,2)~=0
            error('For geometry_mode=''dumbbells'', P must be even.');
        end
        q = createDumbbells(P,delta);
        q = q(:);
    case 'cluster'
        q = grow_cluster(P,delta,2);
        q = q(:);
    otherwise
        error('Unknown geometry_mode: %s',mode);
end
end

function UW = packUW(U,W)
P = size(U,1);
UW = zeros(3*P,1);
UW(1:3:end) = U(:,1);
UW(2:3:end) = U(:,2);
UW(3:3:end) = W;
end

function [U,W] = unpackUW(UW)
P = numel(UW)/3;
U = [UW(1:3:end), UW(2:3:end)];
W = UW(3:3:end);
assert(size(U,1)==P,'Unexpected size in unpackUW.')
end

function FT = packFT(F,T)
P = size(F,1);
FT = zeros(3*P,1);
FT(1:3:end) = F(:,1);
FT(2:3:end) = F(:,2);
FT(3:3:end) = T;
end

function [F,T] = unpackFT(FT)
P = numel(FT)/3;
F = [FT(1:3:end), FT(2:3:end)];
T = FT(3:3:end);
assert(size(F,1)==P,'Unexpected size in unpackFT.')
end

function e = relerr(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end

function printDivider(title_str)
bar = repmat('=',1,78);
fprintf('%s\n',bar);
fprintf('%s\n',bar);
fprintf('%s\n\n',title_str);
end

function printPauseDivider()
bar = repmat('-',1,78);
fprintf('\n%s\n',bar);
fprintf('%s\n',bar);
end

function printSummaryHeader(test_name)
bar = repmat('-',1,78);
fprintf('\n%s\n',bar);
fprintf('Summary of test: %s\n',test_name);
fprintf('%s\n',bar);
end

function plotGeometryPreview(q,rad,title_str)
t = linspace(0,2*pi,240);
figure('Name','test_mob_res geometry','Color','w');
clf;
hold on;

P = numel(q);
show_labels = P <= 40;
for i = 1:P
    x = real(q(i)) + rad(i)*cos(t);
    y = imag(q(i)) + rad(i)*sin(t);
    plot(x,y,'k-','LineWidth',1.3);
    plot(real(q(i)),imag(q(i)),'k.','MarkerSize',12);
    if show_labels
        text(real(q(i)),imag(q(i)),sprintf(' %d',i), ...
            'FontSize',10,'Color',[0.10 0.10 0.10], ...
            'VerticalAlignment','bottom','HorizontalAlignment','left');
    end
end

axis equal;
grid on;
xlabel('x');
ylabel('y');
title(title_str,'Interpreter','none');
hold off;
drawnow;
end
