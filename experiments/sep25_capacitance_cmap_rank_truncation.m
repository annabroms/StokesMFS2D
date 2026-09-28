%SEP25_CAPACITANCE_CMAP_RANK_TRUNCATION Spatial (low-rank) structure of the
% Cmap pair map at several representative gaps.
%
% Two identical discs of radius R are placed on the x axis with gap delta
% between their boundaries (centers at 0 and 2R+delta). For each of several
% representative gaps, the canonical coarse-to-coarse pair map A=Cmap for
% that pair is built directly with getPairBasisLaplace, using the default
% coarse/fine grid resolutions from getLaplace2Dparams. The SVD
%   A = U*Sigma*V'
% gives, for every rank r, the best rank-r approximation A_r (Eckart-Young),
% with relative truncation errors
%   ||A-A_r||_2 / ||A||_2 = sigma_{r+1} / sigma_1
%   ||A-A_r||_F / ||A||_F = sqrt(sum_{j>r} sigma_j^2) / sqrt(sum_j sigma_j^2)
% Plotting these against r for every gap shows how the numerical rank of
% Cmap (i.e. how compressible the pair interaction is) depends on the gap.
%
% Anna Broms, Sep 2026

close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;

fprintf('=== sep25_capacitance_cmap_rank_truncation ===\n\n');

%% Configuration
R = 2;
P = 2;

% Default coarse/fine grid resolutions (N_c=80, N_f=150).
opt = getLaplace2Dparams(P,R);
opt.cmap = 1;
opt.reuse_pair_basis_by_sep = false; % build Cmap directly for this one pair
opt.parallel_precomp = false;
opt.show_counter = 0;
opt.delta_pair = 0.5*R; % keep every representative gap inside the pair-correction regime
opt.use_tikhonov = true;
opt.tikhonov_tol = 1e-11;

% Representative gaps, from near-touching to the edge of the pair-
% correction regime.
delta_over_R = [1e-3 3e-3 1e-2 3e-2 1e-1 3e-1];
delta_vals = delta_over_R*R;
n_gap = numel(delta_vals);

fprintf('R=%.3g, N_c=%d, N_f=%d, N_peanut=%d\n', ...
    R,opt.N_c,opt.N_f,opt.N_peanut);
fprintf('Representative gaps delta/R: %s\n\n',mat2str(delta_over_R));

%% Build the discretisation grids that are fixed across delta
N_c = opt.N_c;
N_f = opt.N_f;
a_c = opt.a_c;
a_f = opt.a_f;
Rp_c = opt.Rp_c;
Rp_f = opt.Rp_f;

nout_c = ceil(a_c*N_c);
tout_c = linspace(0,2*pi,nout_c+1)'; tout_c = tout_c(1:end-1);
rbase_out_c = R*(cos(tout_c)+1i*sin(tout_c));

tin_c = linspace(0,2*pi,N_c+1)'; tin_c = tin_c(1:end-1);
rbase_in_c = Rp_c*(cos(tin_c)+1i*sin(tin_c));

tin_f = linspace(0,2*pi,N_f+1)'; tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*(cos(tin_f)+1i*sin(tin_f));

nout_f = ceil(a_f*N_f);
tout_f = linspace(0,2*pi,nout_f+1)'; tout_f = tout_f(1:end-1);
rout_base_f = R*(cos(tout_f)+1i*sin(tout_f));

N = 2*N_c; % Cmap size: N_c nodes per body, two bodies

%% Build Cmap and its SVD at every representative gap
sigma_all = nan(N,n_gap);
gap_labels = arrayfun(@(d) sprintf('\\delta/R=%.3g',d),delta_over_R,'UniformOutput',false);

fprintf('Building Cmap and computing its SVD for each gap...\n');
for g = 1:n_gap
    delta = delta_vals(g);
    q = [0; 2*R+delta];

    [~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
    assert(isequal(pairs,[1 2]),'Expected a single close pair (1,2).');

    [~,~,~,~,Cmap] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f, ...
        rout_base_f,rbase_out_c,rimage_vec,refine,pairs,opt,rbase_in_c);
    A = Cmap{1,2};

    sigma_all(:,g) = svd(A);
    fprintf('  delta/R=%.3g: sigma_1=%.3e, sigma_N=%.3e, cond=%.3e\n', ...
        delta_over_R(g),sigma_all(1,g),sigma_all(end,g), ...
        sigma_all(1,g)/max(sigma_all(end,g),eps));
end
fprintf('Done.\n\n');

%% Best rank-r truncation errors (Eckart-Young), for r=0,...,N
% r=0 means the zero matrix, so both relative errors are 1 there.
r_vals = (0:N)';
err_spectral = nan(N+1,n_gap);
err_frobenius = nan(N+1,n_gap);
for g = 1:n_gap
    sigma = sigma_all(:,g);
    tail_energy = [flipud(cumsum(flipud(sigma.^2))); 0]; % tail_energy(r+1)=sum_{j>r} sigma_j^2, r=0..N
    total_energy = tail_energy(1);
    sigma_tail = [sigma; 0]; % sigma_tail(r+1) = sigma_{r+1} (0 for r=N)

    err_spectral(:,g) = sigma_tail/sigma(1);
    err_frobenius(:,g) = sqrt(tail_energy/total_energy);
end

%% Report the numerical rank needed for a few representative tolerances
tol_report = [1e-3 1e-6 1e-9 1e-12];
fprintf('%-12s','delta/R');
for t = 1:numel(tol_report)
    fprintf('%12s',sprintf('tol=%.0e',tol_report(t)));
end
fprintf('\n');
for g = 1:n_gap
    fprintf('%-12.3g',delta_over_R(g));
    for t = 1:numel(tol_report)
        r_needed = find(err_spectral(:,g) <= tol_report(t),1,'first') - 1;
        fprintf('%12d',r_needed);
    end
    fprintf('\n');
end
fprintf('\n');

%% Visualise the truncation errors against rank, for every gap
figure('Name','Cmap rank truncation errors');
subplot(1,2,1);
semilogy(r_vals,err_spectral,'-');
xlabel('rank r');
ylabel('||A-A_r||_2 / ||A||_2 = \sigma_{r+1}/\sigma_1');
legend(gap_labels,'Location','best');
title('Best rank-r spectral-norm error');
grid on;

subplot(1,2,2);
semilogy(r_vals,err_frobenius,'-');
xlabel('rank r');
ylabel('||A-A_r||_F / ||A||_F');
legend(gap_labels,'Location','best');
title('Best rank-r Frobenius-norm error');
grid on;

fprintf(['A curve that drops steeply (many decades over a small increase\n' ...
    'in r) shows that Cmap is well approximated by a low-rank matrix at\n' ...
    'that gap; curves for smaller gaps that decay more slowly indicate\n' ...
    'that the pair interaction becomes harder to compress as the two\n' ...
    'discs get closer.\n']);
