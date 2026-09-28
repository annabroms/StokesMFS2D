clear;
close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(repo_root,'startup.m'));

problem = 'mobility';

rad = 1;
gap = 0.01;
N_c = 60;
N_f = 150;
N_peanut_fit = 200;
N_peanut_test = 400;
qr_tol = 1e-3;
peanut_svd_tol = 1e-14;
%peanut_svd_tol = 1e-16;
left_weight = false;
column_weight = false;
cmap_check_tol = 1e-12;

problem = lower(char(problem));
if ~ismember(problem,{'mobility','resistance'})
    error('may7_pair_peanut_qr_nodes:BadProblem', ...
        'problem must be ''mobility'' or ''resistance''.');
end

fprintf('=== %s (May 7, 2026) ===\n',mfilename);
fprintf(['problem=%s, gap=%.3f, N_c=%d, N_f=%d, N_peanut_fit=%d, ', ...
    'N_peanut_test=%d, qr_tol=%.1e\n'], ...
    problem,gap,N_c,N_f,N_peanut_fit,N_peanut_test,qr_tol);

opt = get2Dparams(2,N_c,N_f);
opt.rad = rad;
opt.delta_pair = 0.2;
opt.N_peanut = N_peanut_fit;
opt.use_fmm = false;
opt.visualise_grid = 0;
opt.show_counter = 0;
opt.cmap = true;
opt.project_force = strcmp(problem,'mobility');
opt.reuse_pair_basis_by_sep = false;

svd_opts = struct('left_weight',left_weight,'column_weight',column_weight);
q = [0; 2*rad + gap];

rbase_in_c = build_circle_nodes(opt.Rp_c,N_c);
rbase_in_f = build_circle_nodes(opt.Rp_f,N_f);
rout_base_c = build_circle_nodes(rad,ceil(opt.a_c*N_c));
rout_base_f = build_circle_nodes(rad,ceil(opt.a_f*N_f));

[~,~,~,rimage_pairs,refine,pairs] = getEnhancedGrid(q,opt);
assert(size(pairs,1) == 1 && all(pairs(1,:) == [1 2]), ...
    'Expected exactly one close pair.');

refine_i = refine{1,2};
refine_j = refine{2,1};
rimage_i = rimage_pairs{1,2};
rimage_j = rimage_pairs{2,1};

rin_pair_f = [q(1)+rbase_in_f; rimage_i; q(2)+rbase_in_f; rimage_j];
rout_f = [q(1)+rout_base_f; refine_i; q(2)+rout_base_f; refine_j];
rin_pair_c = [q(1)+rbase_in_c; q(2)+rbase_in_c];

svd_pair = svd_opts;
if left_weight
    svd_pair.row_weights = [...
        getPeriodicCurveWeights([q(1)+rout_base_f; refine_i], q(1));...
        getPeriodicCurveWeights([q(2)+rout_base_f; refine_j], q(2))];
end

rout_peanut_fit = createPeanut(q(1),q(2),N_peanut_fit,false,rad);
rout_peanut_test = createPeanut(q(1),q(2),N_peanut_test,false,rad);

if strcmp(problem,'mobility')
    [~,~,Lc] = getSelfPseudoMobilityStokes( ...
        1,q,rbase_in_c,rout_base_c,[],[0 numel(rout_base_c)],svd_opts);
    Lc_pair = getILpair(Lc{1});

    Kf1 = getKmat2D(rin_pair_f(1:end/2),q(1));
    Kf2 = getKmat2D(rin_pair_f(end/2+1:end),q(2));
    B1 = getKmat2D(rout_f(1:end/2),q(1));
    B2 = getKmat2D(rout_f(end/2+1:end),q(2));

    pair_moment_map = getKftPair(Kf1,Kf2);
    pair_rbm_map = pair_moment_map';
    pair_moment_gram = pair_moment_map*pair_rbm_map;
    pair_target_rbm_map = getKftPair(B1,B2)';

    [Uf,Yf] = getPairBlockStokes(rin_pair_f,rout_f, ...
        pair_moment_map,pair_rbm_map,pair_moment_gram, ...
        pair_target_rbm_map,svd_pair);

    Nf_fit = stokSLPmat(rin_pair_f,rout_peanut_fit,1);
    Nf_test = stokSLPmat(rin_pair_f,rout_peanut_test,1);
    Ntot_fit = Nf_fit - (Nf_fit*pair_rbm_map)/pair_moment_gram * ...
        pair_moment_map;
    Ntot_test = Nf_test - (Nf_test*pair_rbm_map)/pair_moment_gram * ...
        pair_moment_map;
    A_fit = stokSLPmat(rin_pair_c,rout_peanut_fit,1) * Lc_pair;
    A_test = stokSLPmat(rin_pair_c,rout_peanut_test,1) * Lc_pair;
else
    [Uf,Yf] = getPairBlockStokes(rin_pair_f,rout_f,[],[],[],[],svd_pair);
    Ntot_fit = stokSLPmat(rin_pair_f,rout_peanut_fit,1);
    Ntot_test = stokSLPmat(rin_pair_f,rout_peanut_test,1);
    A_fit = stokSLPmat(rin_pair_c,rout_peanut_fit,1);
    A_test = stokSLPmat(rin_pair_c,rout_peanut_test,1);
end

Npair = evaluateCoarseOnPair(q,rbase_in_c,rout_f);
Upf = -Uf' * Npair;
rhs_fit = Ntot_fit * (Yf * Upf);
rhs_test = Ntot_test * (Yf * Upf);

svd_peanut = svd_opts;
if left_weight
    svd_peanut.row_weights = getPeriodicCurveWeights(rout_peanut_fit);
end

[Y_full,U_full] = getPseudoFactors(A_fit,peanut_svd_tol,0,svd_peanut);
C_full = Y_full * (U_full' * rhs_fit);

if strcmp(problem,'mobility')
    basis_Lc = Lc{1};
else
    basis_Lc = [];
end
[~,~,~,~,C_map,~,~] = getPairBasisStokes( ...
    q,rbase_in_c,rbase_in_f,rimage_pairs,refine,pairs,opt,basis_Lc, ...
    rout_base_c,svd_opts);
C_map = C_map{1,2};
cmap_err = relerr(C_full,C_map);

node_qr = build_node_qr_matrix(C_map); %A_fit);
[~,R,piv] = qr(node_qr,'vector');
rank_qr = qr_rank_from_R(R,qr_tol);
selected_nodes = balance_selected_nodes(piv(1:rank_qr),N_c);
% selected_nodes = 1:2:2*N_c; 
% a = 20;
%selected_nodes = [1:a N_c-a:N_c 3*N_c/2-a:3*N_c/2+a]; 
%selected_nodes = [N_c/2-a:N_c/2+a N_c:N_c+a 2*N_c-a:2*N_c]; 

if isempty(selected_nodes)
    error('may7_pair_peanut_qr_nodes:EmptySelection', ...
        'QR selection kept no balanced source nodes.');
end

selected_cols = sort([selected_nodes(:); selected_nodes(:)+2*N_c]);
[Y_sel,U_sel] = getPseudoFactors( ...
    A_fit(:,selected_cols),peanut_svd_tol,0,svd_peanut);
C_sel_small = Y_sel * (U_sel' * rhs_fit);
C_sel = zeros(size(C_full));
C_sel(selected_cols,:) = C_sel_small;

fit_full_test = A_test * C_full;
fit_sel_test = A_test * C_sel;
err_full = relerr(fit_full_test,rhs_test);
err_sel = relerr(fit_sel_test,rhs_test);
err_sel_vs_full = relerr(fit_sel_test,fit_full_test);

n_body_1 = sum(selected_nodes <= N_c);
n_body_2 = sum(selected_nodes > N_c);

fprintf('Selected nodes: %d total = %d on body 1 + %d on body 2\n', ...
    numel(selected_nodes),n_body_1,n_body_2);
fprintf('Off-grid peanut relative error, full columns:    %.3e\n',err_full);
fprintf('Off-grid peanut relative error, selected nodes:  %.3e\n',err_sel);
fprintf('Selected-vs-full off-grid relative difference:   %.3e\n', ...
    err_sel_vs_full);

[full_pw,sel_pw] = pointwise_relative_errors( ...
    fit_full_test-rhs_test,fit_sel_test-rhs_test,rhs_test);

figure('Name','May 7 pair peanut QR nodes');
subplot(1,2,1);
plot(real(rout_peanut_fit),imag(rout_peanut_fit),'k-','LineWidth',1.0);
hold on
plot(real(q(1)+rad*exp(1i*linspace(0,2*pi,400))), ...
    imag(q(1)+rad*exp(1i*linspace(0,2*pi,400))),'k--');
plot(real(q(2)+rad*exp(1i*linspace(0,2*pi,400))), ...
    imag(q(2)+rad*exp(1i*linspace(0,2*pi,400))),'k--');
plot(real(rin_pair_c),imag(rin_pair_c),'.','Color',[0.65 0.65 0.65], ...
    'MarkerSize',12);
plot(real(rin_pair_c(selected_nodes)),imag(rin_pair_c(selected_nodes)), ...
    'ro','MarkerFaceColor','r','MarkerSize',6);
axis equal
grid on
title(sprintf('%s: selected coarse nodes (%d = %d + %d)', ...
    problem,numel(selected_nodes),n_body_1,n_body_2), ...
    'Interpreter','none');
xlabel('x');
ylabel('y');
legend('peanut fit boundary','circle 1','circle 2', ...
    'all coarse nodes','selected nodes','Location','best');

subplot(1,2,2);
semilogy(full_pw+eps,'k-','LineWidth',1.2);
hold on
semilogy(sel_pw+eps,'r--','LineWidth',1.2);
grid on
axis tight
xlabel('test peanut node');
ylabel('pointwise relative error');
title(sprintf('off-grid peanut errors: full %.2e, selected %.2e', ...
    err_full,err_sel));
legend('full columns','selected nodes','Location','best');

function z = build_circle_nodes(R,N)
t = linspace(0,2*pi,N+1)';
t = t(1:end-1);
z = R*(cos(t)+1i*sin(t));
end

function B = build_node_qr_matrix(A)
n_nodes = size(A,2)/2;
B = zeros(2*size(A,1),n_nodes);
for k = 1:n_nodes
    B(:,k) = reshape(A(:,[k k+n_nodes]),[],1);
end
end

function rank_qr = qr_rank_from_R(R,qr_tol)
d = abs(diag(R));
if isempty(d) || d(1) == 0
    rank_qr = 0;
else
    rank_qr = sum(d >= qr_tol*d(1));
end
end

function selected = balance_selected_nodes(nodes,N_c)
selected = nodes(:).';
is_body_1 = selected <= N_c;

while sum(is_body_1) ~= sum(~is_body_1)
    if sum(is_body_1) > sum(~is_body_1)
        drop = find(is_body_1,1,'last');
    else
        drop = find(~is_body_1,1,'last');
    end
    selected(drop) = [];
    is_body_1(drop) = [];
end

selected = sort(selected(:));
end

function [full_pw,sel_pw] = pointwise_relative_errors(err_full,err_sel,ref)
n = size(ref,1)/2;
ref_pw = max(hypot(ref(1:n,:),ref(n+1:end,:)),[],2);
scale = max(ref_pw,1e-14*max(1,max(ref_pw)));
full_pw = max(hypot(err_full(1:n,:),err_full(n+1:end,:)),[],2) ./ scale;
sel_pw = max(hypot(err_sel(1:n,:),err_sel(n+1:end,:)),[],2) ./ scale;
end

function e = relerr(a,b)
e = norm(a(:)-b(:),inf) / max(1,norm(b(:),inf));
end
