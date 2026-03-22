%LEFT_RIGHT_PRECOND Compare left- vs right-preconditioned enhanced 1-body
%resistance solves on randomized clusters of P disks with separation delta.
%
% For each geometry/BC case, both solvers use the same enhanced
% discretisation (same sources/targets and boundary data); only the side on
% which the 1-body preconditioner is applied differs.
%
% The script sweeps over:
%   delta = [0.5, 1, 5]
%   P     = 10:10:100
%
% For each (delta,P) pair, 5 randomized cluster/velocity samples are
% generated. The same geometry and boundary data are passed to both
% solvers:
%   solve_res_1B_enhanced_left
%   solve_res_1B_enhanced
%
% Three figures are produced:
%   1. GMRES iteration counts
%   2. Relative boundary residuals
%   3. Relative force/torque mismatch between the two solves
%
% In each plot, the shaded band spans the trialwise min/max values and the
% solid line shows the trial mean.
%
% This is a demo script only. It does not run automatically.

clear;
close all;
clc;

%% Sweep configuration
delta_vec = [0.5, 1, 5];
P_vec = 10:10:100;
%P_vec = [10 50];
%P_vec = 2:1:5;
n_repeat = 5;
base_seed = 8;

%% Solver options
gmres_tol = 1e-8;
maxit = 800;

%% Storage
n_delta = numel(delta_vec);
n_P = numel(P_vec);

iter_left = nan(n_delta,n_P,n_repeat);
iter_right = nan(n_delta,n_P,n_repeat);
res_left = nan(n_delta,n_P,n_repeat);
res_right = nan(n_delta,n_P,n_repeat);
ft_err = nan(n_delta,n_P,n_repeat);

%% Main sweep
for id = 1:n_delta
    delta = delta_vec(id);
    fprintf('=== delta = %.3g ===\n',delta);

    for ip = 1:n_P
        P = P_vec(ip);
        fprintf('  P = %3d\n',P);
        opt = buildDemoOptions(P,gmres_tol,maxit);

        for ir = 1:n_repeat
            seed = base_seed + 10000*id + 100*P + ir;
            [q,U,W] = buildRandomCase(P,delta,seed);

            [FT_left,sol_left] = runSolverQuietly(@solve_res_1B_enhanced_left,q,U,W,opt);
            [FT_right,sol_right] = runSolverQuietly(@solve_res_1B_enhanced,q,U,W,opt);

            iter_left(id,ip,ir) = sol_left.it;
            iter_right(id,ip,ir) = sol_right.it;
            res_left(id,ip,ir) = sol_left.rel_res;
            res_right(id,ip,ir) = sol_right.rel_res;
            ft_err(id,ip,ir) = relerr(FT_left,FT_right);
        end
    end
end

%% Visualisation
left_color = [0.12, 0.47, 0.71];
right_color = [0.85, 0.33, 0.10];
err_color = [0.20, 0.60, 0.25];

plotSolverMetricFigure(P_vec,delta_vec,iter_left,iter_right, ...
    left_color,right_color, ...
    'GMRES Iterations: Left vs Right 1B Enhanced Resistance', ...
    'GMRES iterations',false);

plotSolverMetricFigure(P_vec,delta_vec,res_left,res_right, ...
    left_color,right_color, ...
    'Relative Boundary Residual: Left vs Right 1B Enhanced Resistance', ...
    'relative residual',true);

plotErrorFigure(P_vec,delta_vec,ft_err,err_color, ...
    'Relative Force/Torque Error: Left vs Right 1B Enhanced Resistance');


function opt = buildDemoOptions(P,gmres_tol,maxit)
opt = get2Dparams(P);
opt.rad = 1;
opt.delta_pair = 0.2;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.get_bndry_field = 1;
opt.gmres_tol = gmres_tol;
opt.gmres_verbose = 0;
opt.debug = 0;
opt.maxit = maxit;
opt.use_fmm = true;
end

function [q,U,W] = buildRandomCase(P,delta,seed)
rng(seed);
evalc('q = grow_cluster(P,delta,2);');
q = q(:);
U = randn(P,2);
W = randn(P,1);
end

function [FT,sol] = runSolverQuietly(solver_fun,q,U,W,opt)
evalc('[FT,sol] = solver_fun(q,U,W,opt);');
end

function plotSolverMetricFigure(P_vec,delta_vec,left_data,right_data, ...
    left_color,right_color,fig_title,y_label,use_log)
figure('Color','w','Name',fig_title);
tiledlayout(1,numel(delta_vec),'TileSpacing','compact','Padding','compact');

for id = 1:numel(delta_vec)
    ax = nexttile;
    hold(ax,'on');

    plotBand(ax,P_vec,squeeze(left_data(id,:,:)),left_color,'left');
    plotBand(ax,P_vec,squeeze(right_data(id,:,:)),right_color,'right');

    grid(ax,'on');
    box(ax,'on');
    title(ax,sprintf('\\delta = %.3g',delta_vec(id)));
    ylabel(ax,y_label);
    if use_log
        set(ax,'YScale','log');
    end
    xlabel(ax,'number of particles P');
   
    if id == 1
        legend(ax,'Location','best');
    end
    axis tight
end

sgtitle(fig_title,'Interpreter','none');
end

function plotErrorFigure(P_vec,delta_vec,err_data,err_color,fig_title)
figure('Color','w','Name',fig_title);
tiledlayout(numel(delta_vec),1,'TileSpacing','compact','Padding','compact');

for id = 1:numel(delta_vec)
    ax = nexttile;
    hold(ax,'on');

    plotBand(ax,P_vec,squeeze(err_data(id,:,:)),err_color,'left/right FT error');

    grid(ax,'on');
    box(ax,'on');
    set(ax,'YScale','log');
    title(ax,sprintf('\\delta = %.3g',delta_vec(id)));
    ylabel(ax,'relative FT error');
    if id == numel(delta_vec)
        xlabel(ax,'number of particles P');
    end
    if id == 1
        legend(ax,'Location','best');
    end
end

sgtitle(fig_title,'Interpreter','none');
end

function plotBand(ax,x,data,color,label_text)
x = x(:).';
data = reshape(data,numel(x),[]);

ymin = min(data,[],2,'omitnan');
ymax = max(data,[],2,'omitnan');
ymean = mean(data,2,'omitnan');

valid = ~(isnan(ymin) | isnan(ymax) | isnan(ymean));
if ~any(valid)
    return;
end

xv = x(valid);
ymin = ymin(valid).';
ymax = ymax(valid).';
ymean = ymean(valid).';

fill(ax,[xv, fliplr(xv)],[ymin, fliplr(ymax)],color, ...
    'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
plot(ax,xv,ymean,'Color',color,'LineWidth',1.8,'DisplayName',label_text);
end

function e = relerr(a,b)
den = norm(b);
if den == 0
    e = norm(a-b);
else
    e = norm(a-b)/den;
end
end
