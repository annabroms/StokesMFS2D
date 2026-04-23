close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Random Discs Compare Solvers (Apr 14, 2026) ===\n');

% Random seeds
geom_seed = 1;
load_seed = 11;

% Geometry
P = 10;
rad = 1;
domain = 'boxed';
phi = 0.65;
min_gap = 1e-3; 
n_sweeps = 30;
visualise_geometry = false;

% Solver
N_c = 60;
N_f = 150;

N_c = 120; %120 for small 2-way error
N_f = 60;
N_peanut = 400;
delta_pair = 0.2;
gmres_tol = 1e-8;
maxit = 1000;

geom_opt = struct();
geom_opt.domain = domain;
geom_opt.phi = phi;
geom_opt.rad = rad;
geom_opt.min_gap = min_gap;
geom_opt.n_sweeps = n_sweeps;
geom_opt.rng_seed = geom_seed;
geom_opt.visualise = visualise_geometry;

[q,geom_meta] = random_discs_mc(P,geom_opt);
[n_close,pairs] = count_close_pairs(q,delta_pair,rad); 

opt = get2Dparams(P,N_c,N_f);
opt.rad = rad;
opt.delta_pair = delta_pair;
opt.N_peanut = N_peanut;
opt.gmres_tol = gmres_tol;
opt.maxit = maxit;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.surface_error_mode = 'rel';
opt.reuse_pair_basis_by_sep = false;
opt.show_counter = 1;
opt.cmap = 1;
opt.self_correct = 1;
opt.use_dense = 1;
opt.get_bndry_field = 1;
opt.single_threaded = 0;
opt.mob_big_sparse_build_mode = 'precomputed'; %streaming can't be used with get_bndry_field

rng(load_seed);
F = randn(P,2);
F = F - mean(F,1);
T = randn(P,1);
U = randn(P,2);
W = randn(P,1);

methods = struct( ...
    'name',{'one_body','peanut'}, ...
    'label',{'1-body preconditioning','Peanut compression'}, ...
    'mob_solver',{@solve_mob_1B_enhanced,@solve_mob_peanut_enhanced}, ...
    'res_solver',{@solve_res_1B_enhanced,@solve_res_peanut_enhanced});

results = repmat(struct( ...
    'name','', ...
    'label','', ...
    'UW',[], ...
    'FT',[], ...
    'mob_sol',struct(), ...
    'res_sol',struct(), ...
    'mob_time',NaN, ...
    'res_time',NaN, ...
    'mob_gmres_residual',NaN, ...
    'res_gmres_residual',NaN, ...
    'two_way_mob_to_res',NaN, ...
    'two_way_res_to_mob',NaN), numel(methods), 1);

fprintf('Geometry:\n');
fprintf('  domain=%s, P=%d, phi_target=%.6f, phi=%.6f, L=%.6f\n', ...
    geom_meta.domain,P,geom_meta.phi_target,geom_meta.phi,geom_meta.L);
fprintf('  min allowed gap=%.3e, measured min gap=%.3e\n', ...
    geom_meta.min_gap,geom_meta.min_surface_gap);
fprintf('  close pairs below delta_pair=%.3f: %d\n',delta_pair,n_close);
fprintf('  N_c=%d, N_f=%d, N_peanut=%d, gmres_tol=%.1e\n', ...
    opt.N_c,opt.N_f,opt.N_peanut,opt.gmres_tol);

FT_ref = packFT(F,T);
UW_ref = packUW(U,W);
it_start = 1; 
for k = it_start:numel(methods)
    results(k).name = methods(k).name;
    results(k).label = methods(k).label;

    fprintf('\nRunning %s mobility...\n',methods(k).label);
    [results(k).UW,results(k).mob_sol] = methods(k).mob_solver(q,F,T,opt);
    results(k).mob_time = results(k).mob_sol.solve_time.total;
    results(k).mob_residual = results(k).mob_sol.rel_res;

    fprintf('Running %s resistance...\n',methods(k).label);
    [results(k).FT,results(k).res_sol] = methods(k).res_solver(q,U,W,opt);
    results(k).res_time =  results(k).res_sol.solve_time.total;
    results(k).res_residual = results(k).res_sol.rel_res;


    [Uk,Wk] = unpackUW(results(k).UW);
    [FT_back,~] = methods(k).res_solver(q,Uk,Wk,opt);
    results(k).two_way_mob_to_res = relerr_inf(FT_back,FT_ref);

    [Fk,Tk] = unpackFT(results(k).FT);
    [UW_back,~] = methods(k).mob_solver(q,Fk,Tk,opt);
    results(k).two_way_res_to_mob = relerr_inf(UW_back,UW_ref);
end

fprintf('\nMobility results:\n');
fprintf('  %-22s %8s %10s %14s %14s %10s\n', ...
    'solver','it','unknowns','surf_rel','surf_abs','time(s)');
for k = it_start:numel(results)
    fprintf('  %-22s %8d %10d %14.3e %14.3e %10.2f\n', ...
        results(k).label, ...
        results(k).mob_sol.it, ...
        results(k).mob_sol.gmres_unknowns, ...
        results(k).mob_sol.rel_res, ...
        results(k).mob_sol.abs_res, ...
        results(k).mob_time);
end

fprintf('\nResistance results:\n');
fprintf('  %-22s %8s %10s %14s %14s %10s\n', ...
    'solver','it','unknowns','gmres_res','surf_rel','time(s)');
for k = it_start:numel(results)
    fprintf('  %-22s %8d %10d %14.3e %14.3e %10.2f\n', ...
        results(k).label, ...
        results(k).res_sol.it, ...
        results(k).res_sol.gmres_unknowns, ...
        results(k).res_gmres_residual, ...
        results(k).res_sol.rel_res, ...
        results(k).res_time);
end

fprintf('\nTwo-way checks:\n');
fprintf(['  two-way mob->res = ||R(M(F,T)) - [F;T]||_inf / max(1,||[F;T]||_inf)\n']);
fprintf(['  two-way res->mob = ||M(R(U,W)) - [U;W]||_inf / max(1,||[U;W]||_inf)\n']);
fprintf('  %-22s %16s %16s\n','solver','mob->res','res->mob');
for k = it_start:numel(results)
    fprintf('  %-22s %16.3e %16.3e\n', ...
        results(k).label, ...
        results(k).two_way_mob_to_res, ...
        results(k).two_way_res_to_mob);
end

fprintf('\nRatios (1-body / peanut):\n');
fprintf('  mobility iteration ratio   = %.3f\n', ...
    results(1).mob_sol.it / results(2).mob_sol.it);
fprintf('  mobility unknown ratio     = %.3f\n', ...
    results(1).mob_sol.gmres_unknowns / results(2).mob_sol.gmres_unknowns);
fprintf('  mobility solve time ratio   = %.3f\n', ...
    results(1).mob_time / results(2).mob_time);
fprintf('  resistance iteration ratio = %.3f\n', ...
    results(1).res_sol.it / results(2).res_sol.it);
fprintf('  resistance unknown ratio   = %.3f\n', ...
    results(1).res_sol.gmres_unknowns / results(2).res_sol.gmres_unknowns);
fprintf('  resistance solve time ratio   = %.3f\n', ...
    results(1).res_time / results(2).res_time);

labels = {results.label};
mob_iter = arrayfun(@(s) s.mob_sol.it, results);
mob_unknowns = arrayfun(@(s) s.mob_sol.gmres_unknowns, results);
res_iter = arrayfun(@(s) s.res_sol.it, results);
res_unknowns = arrayfun(@(s) s.res_sol.gmres_unknowns, results);
two_way_mr = arrayfun(@(s) s.two_way_mob_to_res, results);
two_way_rm = arrayfun(@(s) s.two_way_res_to_mob, results);

figure('Name','apr14 random discs compare solvers','Color','w');
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

nexttile;
theta = linspace(0,2*pi,200);
hold on;
for j = 1:numel(q)
    z = q(j) + rad*(cos(theta) + 1i*sin(theta));
    plot(real(z),imag(z),'k-','LineWidth',0.8);
end
half_L = geom_meta.L/2;
plot([-half_L half_L half_L -half_L -half_L], ...
     [-half_L -half_L half_L half_L -half_L], ...
     'r--','LineWidth',1.1);
axis equal;
grid on;
xlabel('x');
ylabel('y');
title(sprintf('Random discs: P=%d, close pairs=%d',P,n_close), ...
    'Interpreter','none');

nexttile;
bar(categorical(labels),[mob_iter(:), mob_unknowns(:)]);
ylabel('count');
title('Mobility: iterations and unknowns');
legend({'it','unknowns'},'Location','best');
grid on;

nexttile;
bar(categorical(labels),[res_iter(:), res_unknowns(:)]);
ylabel('count');
title('Resistance: iterations and unknowns');
legend({'it','unknowns'},'Location','best');
grid on;

nexttile;
for k = 1:numel(results)
    semilogy(results(k).mob_sol.resvec,'LineWidth',1.5,'DisplayName',results(k).label);
    hold on
end
grid on;
xlabel('iteration');
ylabel('estimated relative residual');
title('Mobility GMRES residual history');
legend('Location','best');

nexttile;
for k = 1:numel(results)
    semilogy(results(k).res_sol.resvec,'LineWidth',1.5,'DisplayName',results(k).label);
    hold on
end
grid on;
xlabel('iteration');
ylabel('estimated relative residual');
title('Resistance GMRES residual history');
legend('Location','best');

nexttile;
bar(categorical(labels),[two_way_mr(:), two_way_rm(:)]);
set(gca,'YScale','log');
ylabel('two-way error');
title('Round-trip consistency');
legend({'mob->res','res->mob'},'Location','best');
grid on;

title(t,'Random discs: 1-body vs peanut solvers','Interpreter','none');



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
end

function e = relerr_inf(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
