%TEST_CAP_ELAST Compare Laplace capacitance/elastance solver variants.
%
% Capacitance solve in this script is the Modified exterior BVP in d = 2
% described by Stein & Barnett (ACOM, 2022):
%   prescribed per-body voltage v_body -> unknown per-body net charge Q.
%
% Elastance solve:
%   prescribed per-body net charge Q_body -> unknown per-body voltage v_body.
%
% This script compares 1B, 2B and peanut outputs and runs two-way checks
% (cap -> elast and elast -> cap) for each solver family.
%
% Anna Broms, Mar 2026

clear;
close all;
clc;
rng(9);

%% Configuration
R = 2; %particle radius
delta = 1e-3;
%P = 20;
%geometry_mode = 'dumbbells';       % 'line' | 'dumbbells' | 'cluster'
%q = buildGeometry(geometry_mode,P,delta,R);
geometry_mode = 'close packing';
mode = 2;
q = hexagonal_lattice(delta,mode,R);
P = numel(q);
rad = R*ones(P,1);

delta_pair = 0.2;             % near-pair threshold for 2B/peanut
N_peanut = 400;               % points per peanut separation surface
visualise = 0;                % visualise grid and magnitude of sources etc
gmres_tol = 1e-10;
debug = 0; % determine system matrix by using matvec with columns of the identity matrix as input
gmres_verbose = 0;            % 0=silent, 1=summary, 2=per-iteration
use_fmm = true;

fprintf('=== test_cap_elast ===\n');
fprintf('geometry=%s, P=%d, R=%.3f, delta=%.2e, delta_pair=%.2e, N_peanut=%d\n\n', ...
    geometry_mode,P,R,delta,delta_pair,N_peanut);

%% Geometry preview
printDivider('Geometry preview');
plotGeometryPreview(q,rad,sprintf('test_cap_elast geometry (%s, P=%d, R=%.3f)',geometry_mode,P,R));

printPauseDivider();
disp('Press key to continue...')
pause();

%% Capacitance: v_body -> Q
printDivider('Capacitance comparison (input v_body -> output Q_body)');
v_ref = randn(P,1);
opt = getLaplace2Dparams(P,R);
opt_cap = opt;
opt_cap.visualise = visualise;
opt_cap.gmres_tol = gmres_tol;
opt_cap.debug = debug;
opt_cap.use_fmm = use_fmm;
opt_cap.gmres_verbose = gmres_verbose;
opt_cap.delta_pair = delta_pair;
opt_cap.N_peanut = N_peanut;

[Q_1B,sol_cap_1B] = solve_cap_1B(q,v_ref,opt_cap);
[Q_2B,sol_cap_2B] = solve_cap_2B(q,v_ref,opt_cap);
[Q_P,sol_cap_P] = solve_cap_peanut(q,v_ref,opt_cap);
gmres_iter_cap_1B = sol_cap_1B.it;
gmres_iter_cap_2B = sol_cap_2B.it;
gmres_iter_cap_peanut = sol_cap_P.it;
maxres_cap_1B = sol_cap_1B.maxres;
maxres_cap_2B = sol_cap_2B.maxres;
maxres_cap_peanut = sol_cap_P.maxres;

printSummaryHeader('Capacitance comparison');
fprintf('Capacitance (input v_body)\n');
fprintf('  %-8s %12s %14s\n','solver','gmres_iter','maxres');
fprintf('  %-8s %12d %14.3e\n','1B',gmres_iter_cap_1B,maxres_cap_1B);
fprintf('  %-8s %12d %14.3e\n','2B',gmres_iter_cap_2B,maxres_cap_2B);
fprintf('  %-8s %12d %14.3e\n','peanut',gmres_iter_cap_peanut,maxres_cap_peanut);
fprintf('  diff   : ||Q_2B-Q_1B||/||Q_1B|| = %.3e\n',relerr(Q_2B,Q_1B));
fprintf('           ||Q_peanut-Q_2B||/||Q_2B|| = %.3e\n\n',relerr(Q_P,Q_2B));

printPauseDivider();
disp('Press key to continue...')
pause();

%% Elastance: Q_body -> v_body
printDivider('Elastance comparison (input Q_body -> output v_body)');
Q_ref = randn(P,1);
opt_elast = opt_cap;
[v_1B,sol_elast_1B] = solve_elast_1B(q,Q_ref,opt_elast);
[v_2B,sol_elast_2B] = solve_elast_2B(q,Q_ref,opt_elast);
[v_P,sol_elast_P] = solve_elast_peanut(q,Q_ref,opt_elast);
gmres_iter_elast_1B = sol_elast_1B.it;
gmres_iter_elast_2B = sol_elast_2B.it;
gmres_iter_elast_peanut = sol_elast_P.it;
maxres_elast_1B = sol_elast_1B.maxres;
maxres_elast_2B = sol_elast_2B.maxres;
maxres_elast_peanut = sol_elast_P.maxres;

printSummaryHeader('Elastance comparison');
fprintf('Elastance (input Q_body)\n');
fprintf('  %-8s %12s %14s\n','solver','gmres_iter','maxres');
fprintf('  %-8s %12d %14.3e\n','1B',gmres_iter_elast_1B,maxres_elast_1B);
fprintf('  %-8s %12d %14.3e\n','2B',gmres_iter_elast_2B,maxres_elast_2B);
fprintf('  %-8s %12d %14.3e\n','peanut',gmres_iter_elast_peanut,maxres_elast_peanut);
fprintf('  diff   : ||v_2B-v_1B||/||v_1B|| = %.3e\n',relerr(v_2B,v_1B));
fprintf('           ||v_peanut-v_2B||/||v_2B|| = %.3e\n\n',relerr(v_P,v_2B));

printPauseDivider();
disp('Press key to continue...')
pause();

%% Two-way checks
printDivider('Two-way checks');
% Capacitance -> Elastance
opt_tw = opt_cap;
opt_tw.visualise = 0;
opt_tw.debug = 0;
[Q_cap_1B,~] = solve_cap_1B(q,v_ref,opt_tw);
[v_back_1B,~] = solve_elast_1B(q,Q_cap_1B,opt_tw);

[Q_cap_2B,~] = solve_cap_2B(q,v_ref,opt_tw);
[v_back_2B,~] = solve_elast_2B(q,Q_cap_2B,opt_tw);

[Q_cap_P,~] = solve_cap_peanut(q,v_ref,opt_tw);
[v_back_P,~] = solve_elast_peanut(q,Q_cap_P,opt_tw);

% Elastance -> Capacitance
[v_el_1B,~] = solve_elast_1B(q,Q_ref,opt_tw);
[Q_back_1B,~] = solve_cap_1B(q,v_el_1B,opt_tw);

[v_el_2B,~] = solve_elast_2B(q,Q_ref,opt_tw);
[Q_back_2B,~] = solve_cap_2B(q,v_el_2B,opt_tw);

[v_el_P,~] = solve_elast_peanut(q,Q_ref,opt_tw);
[Q_back_P,~] = solve_cap_peanut(q,v_el_P,opt_tw);

two_way_cap_elast_1B = relerr(v_back_1B,v_ref);
two_way_cap_elast_2B = relerr(v_back_2B,v_ref);
two_way_cap_elast_peanut = relerr(v_back_P,v_ref);
two_way_elast_cap_1B = relerr(Q_back_1B,Q_ref);
two_way_elast_cap_2B = relerr(Q_back_2B,Q_ref);
two_way_elast_cap_peanut = relerr(Q_back_P,Q_ref);

printSummaryHeader('Two-way checks');
fprintf('  %-8s %16s %16s\n','solver','cap->elast','elast->cap');
fprintf('  %-8s %16.3e %16.3e\n','1B',two_way_cap_elast_1B,two_way_elast_cap_1B);
fprintf('  %-8s %16.3e %16.3e\n','2B',two_way_cap_elast_2B,two_way_elast_cap_2B);
fprintf('  %-8s %16.3e %16.3e\n','peanut',two_way_cap_elast_peanut,two_way_elast_cap_peanut);



function q = buildGeometry(mode,P,delta,R)
switch lower(mode)
    case 'line'
        spacing = R*(2 + delta);
        q = (0:P-1)'*spacing;
        q = complex(q,zeros(size(q)));
    case 'dumbbells'
        if mod(P,2)~=0
            error('For geometry_mode=''dumbbells'', P must be even.');
        end
        q = createDumbbells(P,delta,R);
        q = q(:);
    case 'cluster'
        q = grow_cluster(P,delta,2,R);
        q = q(:);
    otherwise
        error('Unknown geometry_mode: %s',mode);
end
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
figure('Name','test_cap_elast geometry','Color','w');
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
