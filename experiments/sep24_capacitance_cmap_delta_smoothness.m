%% 
%SEP24_CAPACITANCE_CMAP_DELTA_SMOOTHNESS Smoothness of the Cmap pair map in
% the gap parameter delta.
%
% Two identical discs of radius R are placed on the x axis with gap delta
% between their boundaries (centers at 0 and 2R+delta). For each delta the
% canonical coarse-to-coarse pair map Cmap for that pair is built directly
% with getPairBasisLaplace, using the default coarse/fine grid resolutions
% from getLaplace2Dparams. A handful of representative entries of Cmap, and
% its Frobenius norm, are then plotted against delta on a logarithmic axis
% to check whether they vary smoothly with log(delta/R) or show jumps.
% Two further views cover all N*N entries at once: a heatmap of every
% entry (each rescaled to [0,1]) against log(delta/R), and the power
% spectrum of each entry's delta-trace. The delta range is split at the
% ellipse-segment enhancement on/off switch (see delta_star below) before
% the spectrum is computed, so the large discretisation jump at that
% switch does not mask the intrinsic smoothness within each regime.
%
% Anna Broms, Sep 2026

close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
startup;

fprintf('=== sep24_capacitance_cmap_delta_smoothness ===\n\n');

%% Configuration
R = 2;
P = 2;

% Default coarse/fine grid resolutions (N_c=80, N_f=150).
opt = getLaplace2Dparams(P,R);
opt.cmap = 1;
opt.reuse_pair_basis_by_sep = false; % build Cmap directly for this one pair
opt.parallel_precomp = false;
opt.show_counter = 0;
opt.delta_pair = 0.5*R; % keep every tested delta inside the pair-correction regime
opt.N_f = 150;
 

n_delta = 300;
delta_min = 1e-2*R;
delta_max = 0.3*R;
delta_vals = logspace(log10(delta_min),log10(delta_max),n_delta);

% Set to true to freeze the ellipse-segment discretisation at the smallest
% delta swept here, instead of rebuilding it for every actual gap. This
% removes the node-count-driven kinks diagnosed below, at the cost of
% using an unnecessarily fine enhancement grid for the larger gaps.
opt.ellipse_constant = 1;


fprintf('R=%.3g, N_c=%d, N_f=%d, N_peanut=%d\n', ...
    R,opt.N_c,opt.N_f,opt.N_peanut);
fprintf('delta swept from %.3e to %.3e (%d log-spaced points)\n', ...
    delta_min,delta_max,n_delta);
fprintf('opt.ellipse_constant = %d\n\n',opt.ellipse_constant);

% getEnhancedGrid.m switches the whole ellipse-segment enhancement for a
% body on/off depending on whether its accumulation point has moved
% outside the fine proxy radius Rp_f (see pair_clusters_ellipse.m). For
% equal radii this happens analytically at delta_star below; expect the
% largest kink in the plots near this value (unless opt.ellipse_constant
% freezes the discretisation, in which case it should disappear).
delta_star = (R-opt.Rp_f)^2/opt.Rp_f;
fprintf(['Predicted enhancement on/off switch at delta_star/R = %.3e\n' ...
    '(accumulation point crosses the fine proxy radius Rp_f)\n\n'],delta_star/R);

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

%% Track a handful of representative Cmap entries as delta varies
% Cmap is the 2*N_c x 2*N_c pair map, rows/cols 1:N_c for body 1 and
% N_c+1:2*N_c for body 2. Node 1 (angle 0) on body 1 and node idx_facing
% (angle pi) on body 2 are the nodes closest to the gap, so their entries
% are the ones most sensitive to delta.
%
% Assign track_idx as an n x 2 array of [row col] indices (each in
% 1:2*N_c) before running this script to choose different entries, e.g.:
%track_idx = [1 1; 5 130];
track_idx = [randi(2*N_c,10,1) randi(2*N_c,10,1)];
idx_facing = round(N_c/2)+1;
if ~exist('track_idx','var') || isempty(track_idx)
track_idx = [1 1; ... % body-1 self, nearest-to-gap node
    1 N_c+idx_facing; ... % cross-body, both nodes nearest the gap
    N_c+idx_facing N_c+idx_facing; ... % body-2 self, nearest-to-gap node
    idx_facing N_c+1]; % cross-body, both nodes farthest from the gap
end
validateattributes(track_idx,{'numeric'},{'ncols',2,'integer','positive','<=',N_c*2}, ...
    mfilename,'track_idx');
n_track = size(track_idx,1);
track_labels = arrayfun(@(k) sprintf('C(%d,%d)',track_idx(k,1),track_idx(k,2)), ...
    1:n_track,'UniformOutput',false);

entries = nan(n_delta,n_track);
fro_norm = nan(n_delta,1);
max_abs = nan(n_delta,1);
N = 2*N_c; % Cmap size: N_c nodes per body, two bodies
C_all = nan(N,N,n_delta); % every entry, kept for the "all entries" views below
n_image = zeros(n_delta,2); % ellipse-segment source-node count feeding body 1/2
rank_fine = nan(n_delta,1); % numerical rank of the fine pair block (getPairBlockLaplace)
rank_peanut = nan(n_delta,1); % numerical rank of the coarse-to-peanut block (getPeanutBlockLaplace)
svd_tol = 1e-14; % matches the hardcoded tolerance in both blocks

fprintf('Building Cmap for each delta...\n');
for k = 1:n_delta
    delta = delta_vals(k);
    q = [0; 2*R+delta];

    [~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
    assert(isequal(pairs,[1 2]),'Expected a single close pair (1,2).');
    n_image(k,1) = numel(rimage_vec{1,2});
    n_image(k,2) = numel(rimage_vec{2,1});

    [Uf,~,~,~,Cmap] = getPairBasisLaplace(q,rbase_in_c,rbase_in_f, ...
        rout_base_f,rbase_out_c,rimage_vec,refine,pairs,opt,rbase_in_c);
    C = Cmap{1,2};
    rank_fine(k) = size(Uf{1,2},1); % Uf{i,j} = -U_A'*Npair, so rows = rank(A)

    % Not returned by getPairBasisLaplace when opt.cmap=1, so replicate the
    % same matrix getPeanutBlockLaplace.m factorises, purely to read off
    % its numerical rank (no effect on the production Cmap computation).
    rin_pair_c = [q(1)+rbase_in_c; q(2)+rbase_in_c];
    rout_peanut = createPeanut(q(1),q(2),opt.N_peanut,0,R);
    s_peanut = svd(lapSLPmat(rin_pair_c,rout_peanut));
    rank_peanut(k) = sum(s_peanut > max(s_peanut)*svd_tol);

    for m = 1:n_track
        entries(k,m) = C(track_idx(m,1),track_idx(m,2));
    end
    fro_norm(k) = norm(C,'fro');
    max_abs(k) = max(abs(C(:)));
    C_all(:,:,k) = C;
end
fprintf('Done.\n\n');

%% Where does the enhancement-grid node count change?
% Every such change is a candidate kink location: the first jump from 0 is
% the on/off switch of the whole ellipse-segment mechanism (matches
% delta_star above); later +/-1 changes are individual ellipse-segment
% nodes crossing the r_proxy filter in ellipse_cheb_segment.m as the
% ellipse geometry shifts with delta. opt.compress_cmap is 0 here, so the
% SVD-rank-truncation branch in buildLaplacePairGroup.m is inactive; if
% enabled, its rank cutoff (S > max(S)*opt.cmap_tol) is a further source of
% discrete jumps by the same mechanism (an integer count changing as a
% threshold is crossed).
switch_idx = find(any(diff(n_image) ~= 0,2));
delta_switch = delta_vals(switch_idx+1);
fprintf('Enhancement-grid node count changes at %d of %d delta samples.\n', ...
    numel(switch_idx),n_delta);
fprintf('First few switches, delta/R: %s\n\n', ...
    mat2str(round(delta_switch(1:min(5,end))/R,4)));

%% Where does the truncated-SVD rank in the pair blocks change?
% getPairBlockLaplace.m and getPeanutBlockLaplace.m both call
% getPseudoFactors(...,1e-14,...), which keeps singular values sigma with
% sigma > max(sigma)*1e-14. As delta varies, the pair/peanut matrices'
% singular values drift smoothly, but every time one of them crosses this
% fixed threshold the retained rank changes by +/-1 and the pseudoinverse
% factors (hence Cmap) jump discontinuously, even with opt.ellipse_constant
% freezing the node counts above. Unlike delta_star, these thresholds have
% no simple closed form and must be located empirically, as done here.
rank_switch_idx = find(diff(rank_fine) ~= 0 | diff(rank_peanut) ~= 0);
delta_rank_switch = delta_vals(rank_switch_idx+1);
fprintf('Pseudoinverse rank changes at %d of %d delta samples.\n', ...
    numel(rank_switch_idx),n_delta);
fprintf('First few rank switches, delta/R: %s\n\n', ...
    mat2str(round(delta_rank_switch(1:min(5,end))/R,4)));

%% Simple jump diagnostic: largest relative change between neighbouring
% log-spaced delta samples, for each tracked entry. The step is scaled by
% the entry's overall range (not its local value), because scaling by the
% local value blows up near a zero crossing of an otherwise smooth curve
% and falsely flags it as a jump. "explained" checks whether the largest
% jump sits within one delta sample of a known switch (ellipse node count
% or pseudoinverse rank); "false" means it is not accounted for by either
% mechanism tracked above (and is worth inspecting by hand).
known_switch_idx = unique([switch_idx(:); rank_switch_idx(:)]);
fprintf('%-10s %14s %14s %10s %10s\n', ...
    'entry','max rel. step','argmax delta/R','nearest sw.','explained');
for m = 1:n_track
    step = abs(diff(entries(:,m)));
    scale = max(max(abs(entries(:,m))),eps);
    rel_step = step/scale;
    [val,i_max] = max(rel_step);
    if isempty(known_switch_idx)
        nearest_gap = inf;
    else
        nearest_gap = min(abs(known_switch_idx-i_max));
    end
    fprintf('%-10s %14.3e %14.4g %10d %10s\n',track_labels{m}, ...
        val,delta_vals(i_max+1)/R,nearest_gap,mat2str(nearest_gap<=1));
end
fprintf('\n');

%% Visualise entries vs delta on a log axis (equivalent to vs log(delta/R))
% Grey dashed lines mark every enhancement-grid node-count change, so
% kinks that line up with them are explained by the discretisation switch
% rather than by the underlying (smooth) physics.
figure('Name','Cmap entries vs delta');
subplot(2,1,1);
semilogx(delta_vals/R,entries,'-o','MarkerSize',3);
hold on;
if ~isempty(delta_switch)
    xline(delta_switch/R,':','Color',[0.6 0.6 0.6],'HandleVisibility','off');
end
if ~isempty(delta_rank_switch)
    xline(delta_rank_switch/R,'--','Color',[0.85 0.33 0.1],'HandleVisibility','off');
end
xlabel('\delta / R');
ylabel('C_{map} entry value');
legend(track_labels,'Location','best');
title('Representative Cmap entries vs gap \delta');
grid on;

subplot(2,1,2);
semilogx(delta_vals/R,fro_norm,'-o','MarkerSize',3,'DisplayName','||Cmap||_F');
hold on;
semilogx(delta_vals/R,max_abs,'-s','MarkerSize',3,'DisplayName','max|Cmap|');
if ~isempty(delta_switch)
    xline(delta_switch/R,':','Color',[0.6 0.6 0.6],'HandleVisibility','off');
end
if ~isempty(delta_rank_switch)
    xline(delta_rank_switch/R,'--','Color',[0.85 0.33 0.1],'HandleVisibility','off');
end
xlabel('\delta / R');
ylabel('aggregate magnitude');
legend('Location','best');
title('Cmap aggregate magnitude vs gap \delta');
grid on;

figure('Name','Enhancement-grid node count vs delta');
stairs(delta_vals/R,n_image(:,1),'DisplayName','body 1 source nodes');
hold on;
stairs(delta_vals/R,n_image(:,2),'DisplayName','body 2 source nodes');
set(gca,'XScale','log');
xlabel('\delta / R');
ylabel('# ellipse-segment source nodes kept');
legend('Location','best');
title('Discretisation switches responsible for the kinks above');
grid on;

figure('Name','Pseudoinverse rank vs delta');
stairs(delta_vals/R,rank_fine,'DisplayName','rank(fine pair block)');
hold on;
stairs(delta_vals/R,rank_peanut,'DisplayName','rank(coarse-to-peanut block)');
set(gca,'XScale','log');
xlabel('\delta / R');
ylabel('numerical rank (tol=1e-14)');
legend('Location','best');
title('Truncated-SVD rank changes: a second, independent source of kinks');
grid on;

fprintf('Inspect the figure: smooth curves in these semilogx plots indicate\n');
fprintf('that Cmap entries vary smoothly with log(delta/R); a visible kink\n');
fprintf('or discontinuity indicates a jump at that separation. Grey dotted\n');
fprintf('lines mark ellipse node-count switches; orange dashed lines mark\n');
fprintf('pseudoinverse rank switches (see the two figures above).\n\n');

%% Visualise every Cmap entry at once
% Each of the N*N entries is one row here, rescaled to [0,1] over its own
% range so entries of very different magnitude stay visible on one plot.
% A jump at some delta shows up as a sharp vertical edge across many rows.
V = reshape(C_all,N*N,n_delta);
V_range = max(V,[],2) - min(V,[],2);
V_range(V_range==0) = 1;
V_norm = (V - min(V,[],2))./V_range;

figure('Name','All Cmap entries vs delta');
imagesc(log10(delta_vals/R),1:size(V,1),V_norm);
set(gca,'YDir','normal');
xlabel('log_{10}(\delta/R)');
ylabel('flattened (row,col) entry index');
title('Every Cmap entry (rescaled to [0,1]) vs gap \delta');
colorbar;

%% Fourier spectrum in log(delta/R), split by discretisation regime
% The enhancement on/off switch (see delta_star above) is itself a large
% discontinuity, so one FFT over the whole delta range mixes that
% discretisation artifact with the smoothness of the underlying physics.
% Splitting the delta axis at the switch and computing the spectrum
% separately on each side isolates the two: within a side no further
% on/off switch occurs (though the smaller +/-1 ellipse-node changes seen
% above may still leave some residual high-frequency content on the
% enhanced side). x=log10(delta/R) stays uniformly spaced within each
% side because it is a sub-range of the original log-spaced grid.
is_enhanced = n_image(:,1) > 0 | n_image(:,2) > 0;
boundary = find(~is_enhanced,1,'first');
if isempty(boundary)
    regimes = struct('name',{'ellipse-enhanced'},'idx',{1:n_delta});
elseif boundary == 1
    regimes = struct('name',{'fine-proxy-only'},'idx',{1:n_delta});
else
    regimes = struct('name',{'ellipse-enhanced','fine-proxy-only'}, ...
        'idx',{1:boundary-1,boundary:n_delta});
end

min_pts = 8;
for r = 1:numel(regimes)
    plot_regime_spectrum(regimes(r).idx,regimes(r).name,entries,V,delta_vals,R, ...
        track_labels,min_pts);
end

fprintf(['A power spectrum decaying over many decades indicates a smooth\n' ...
    'function of log(delta/R) within that regime; a spectrum with a\n' ...
    'slowly decaying tail indicates a jump or kink inside the regime\n' ...
    '(the enhancement on/off switch itself is excluded by the split).\n\n']);

%% Smoothness test on every subinterval bounded by a known switch
% This refines the two-regime split above: known_switch_idx merges the
% ellipse node-count switches and the pseudoinverse rank switches, so each
% subinterval below has neither kind of discretisation change inside it.
% If Cmap is intrinsically smooth in log(delta/R), the power spectrum on
% every such subinterval should decay quickly; a slowly-decaying tail on a
% subinterval that contains no known switch would point to a further,
% still-unaccounted-for source of non-smoothness.
edges = [0; known_switch_idx(:); n_delta];
n_seg = numel(edges)-1;
fprintf('Smoothness test on %d subinterval(s) bounded by the %d known switches.\n\n', ...
    n_seg,numel(known_switch_idx));
for s = 1:n_seg
    idx = (edges(s)+1):edges(s+1);
    name = sprintf('subinterval %d/%d, delta/R in [%.3g, %.3g]', ...
        s,n_seg,delta_vals(idx(1))/R,delta_vals(idx(end))/R);
    plot_regime_spectrum(idx,name,entries,V,delta_vals,R,track_labels,min_pts);
end

function plot_regime_spectrum(idx,name,entries,V,delta_vals,R,track_labels,min_pts)
%PLOT_REGIME_SPECTRUM Detrended FFT power spectrum in log(delta/R) on one
% index range, for the tracked entries and for every Cmap entry at once.
if numel(idx) < min_pts
    fprintf('Skipping %s: only %d delta samples (need >= %d).\n', ...
        name,numel(idx),min_pts);
    return
end

x_r = log10(delta_vals(idx)/R);
dx_r = x_r(2)-x_r(1);
assert(max(abs(diff(x_r(:))-dx_r)) < 1e-10*abs(dx_r), ...
    'delta_vals must be log-spaced within each regime for this diagnostic.');
n_r = numel(idx);
freq_r = (0:floor(n_r/2))'/(n_r*dx_r); % cycles per unit log10(delta/R)

entries_r = detrend(entries(idx,:));
Y_track_r = fft(entries_r);
P_track_r = abs(Y_track_r(1:floor(n_r/2)+1,:)).^2/n_r;

V_r = detrend(V(:,idx).'); % n_r x (N*N), detrended along delta
Y_all_r = fft(V_r);
P_all_r = abs(Y_all_r(1:floor(n_r/2)+1,:)).^2/n_r;
spec_mean_r = mean(P_all_r,2);
spec_max_r = max(P_all_r,[],2);

figure('Name',sprintf('Cmap spectrum in log(delta/R): %s',name));
subplot(1,2,1);
semilogy(freq_r,P_track_r,'-o','MarkerSize',3);
xlabel('frequency (cycles per unit log_{10}(\delta/R))');
ylabel('power');
legend(track_labels,'Location','best');
title(sprintf('Tracked entries (%s, %d pts)',name,n_r));
grid on;

subplot(1,2,2);
semilogy(freq_r,spec_mean_r,'-','DisplayName','mean over all entries');
hold on;
semilogy(freq_r,spec_max_r,'--','DisplayName','max over all entries');
xlabel('frequency (cycles per unit log_{10}(\delta/R))');
ylabel('power');
legend('Location','best');
title(sprintf('All entries (%s, %d pts)',name,n_r));
grid on;
end
