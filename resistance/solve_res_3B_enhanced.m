function [FT,sol] = solve_res_3B_enhanced(q,U,W,opt)
%SOLVE_RES_3B_ENHANCED Solve a 2D Stokes resistance problem with triplet corrections.

if nargin == 0
    test_solve_res_3B;
    return
end

if nargin < 4 || ~isstruct(opt)
    error('solve_res_3B_enhanced requires q, U, W, and an options struct opt.');
end

q = q(:);
W = W(:);
P = numel(q);
assert(size(U,1)==P,'U must have one row per particle.');
assert(size(U,2)==2,'U must have two columns [Ux, Uy].');
assert(numel(W)==P,'W must have one entry per particle.');

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-10);
debug = logical(getOptField(opt,'debug',false));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
maxit = getOptField(opt,'maxit',800);
use_fmm = logical(getOptField(opt,'use_fmm',true));
triplet_mode = lower(char(getOptField(opt,'triplet_accumulation_mode','edge_weighted')));
triplet_lsq_mode = lower(char(getOptField(opt,'triplet_lsq_mode','full')));

fprintf('==== START: %s ====\n', mfilename);

opt.P = P;
opt.rad = 1;
opt.project_force = false;
opt.gmres_verbose = gmres_verbose;
opt.use_fmm = use_fmm;
opt.N_peanut = 0;
opt.cmap = 0;

N_c = getOptField(opt,'N_c',60);
N_f = getOptField(opt,'N_f',60);
a_c = getOptField(opt,'a_c',1.2);
a_f = getOptField(opt,'a_f',1.2);
tol_c = getOptField(opt,'tol_c',1e-12);
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',max([1-sep_f,0.01]));
delta_pair = getOptField(opt,'delta_pair',(1-Rp_c)^2/Rp_c);

tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1)';
tout_c = tout_c(1:end-1);
rbase_out_c = cos(tout_c)+1i*sin(tout_c);

tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
rbase_in_c = Rp_c*cos(tin_c)+1i*Rp_c*sin(tin_c);

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*cos(tin_f)+1i*Rp_f*sin(tin_f);

rout = zeros(P*numel(rbase_out_c),1);
for k = 1:P
    rout((k-1)*numel(rbase_out_c)+1:k*numel(rbase_out_c)) = q(k)+rbase_out_c;
end

opt.delta_pair = delta_pair;
[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
triplets = getCloseTripletsFromPairs(pairs,P);
pair_lookup = getPairRowLookup(pairs,P);
if strcmp(triplet_mode,'pair_subtracted')
    triplet_weights = [];
else
    triplet_weights = getTripletAccumulationWeights(triplets,P,triplet_mode);
end
[UU,YY] = getSelfPseudo(1,rbase_in_c,rbase_out_c);

pair_Upf = [];
pair_Ypf = [];
pair_cache = struct('enabled',false);
pair_test_data = [];
if strcmp(triplet_mode,'pair_subtracted')
    [pair_Upf,pair_Ypf,~,~,~,~,pair_cache] = ...
        getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,[],rbase_out_c);
end

need_triplet_pair_data = ~strcmp(triplet_lsq_mode,'full') || ...
    logical(getOptField(opt,'triplet_basis_debug',0));
if need_triplet_pair_data
    if isfield(pair_cache,'enabled') && pair_cache.enabled
        pair_test_data = struct('pairs',pairs,'pair_lookup',pair_lookup,'pair_cache',pair_cache);
    else
        opt_pair_trip = opt;
        opt_pair_trip.reuse_pair_basis_by_sep = true;
        opt_pair_trip.pair_basis_debug = 0;
        [~,~,~,~,~,~,pair_cache_trip] = getPairBasisStokes( ...
            q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_pair_trip,[],rbase_out_c);
        pair_test_data = struct('pairs',pairs,'pair_lookup',pair_lookup,'pair_cache',pair_cache_trip);
    end
end

triplet_cache = getTripletBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    triplets,opt,[],pair_test_data);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.refine = refine;
geom.rimage_vec = rimage_vec;
geom.opt = opt;
geom.rvec_out = rout;
geom.q = q;
geom.pairs = pairs;
geom.triplets = triplets;
geom.triplet_weights = triplet_weights;
geom.pair_lookup = pair_lookup;

basis = struct();
basis.U = UU;
basis.Y = YY;
basis.triplet_cache = triplet_cache;
basis.Upf = pair_Upf;
basis.Ypf = pair_Ypf;
basis.pair_cache = pair_cache;

foutx = [];
fouty = [];
for k = 1:P
    rhs_x = @(x) U(k,1)-W(k)*(imag(x)-imag(q(k)));
    rhs_y = @(x) U(k,2)+W(k)*(real(x-q(k)));
    seg = (k-1)*numel(rbase_out_c)+1:k*numel(rbase_out_c);
    foutx = [foutx; rhs_x(rout(seg))];
    fouty = [fouty; rhs_y(rout(seg))];
end
fout = [foutx; fouty];

n_bound = 2000;
tcheck = linspace(0,2*pi,n_bound)';
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k)+cos(tcheck)+1i*sin(tcheck);
end

if debug
    ncols = 2*numel(rout);
    CC = zeros(ncols,ncols);
    x = zeros(ncols,1);
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_res_3B_enhanced(x,geom,basis);
    end
    figure();
    imagesc(log10(abs(CC)));
    colorbar;
    title([mfilename ': log_{10} |matvec system matrix|'],'interpreter','none');
end

disp(' == Solving... == ');
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_res_3B_enhanced(x,geom,basis), ...
    fout,2*numel(rout),maxit,gmres_tol,opt.gmres_verbose,rout);

if visualise_sol
    figure();
    semilogy(resvec);
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual');
    axis tight
    grid on
    title('Convergence resistance with triplet corr','interpreter','latex');
end

disp(' == Postprocessing == ');
[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y,~,rimage_k,triplet_diag] = ...
    getTripletTransformationStokes(tau,geom,basis);
lambda = [tau_stokes_x; tau_stokes_y];

FT = zeros(3*P,1);
Kc = getKmat2D(rbase_in_c,0);
for k = 1:P
    FT((k-1)*3+1:3*k) = FT((k-1)*3+1:3*k) + ...
        Kc'*[tau_stokes_x(coarse_ind{k}); tau_stokes_y(coarse_ind{k})];
end

if strcmp(triplet_mode,'pair_subtracted')
    active_particles = sort(unique(pairs(:)));
else
    active_particles = sort(unique(triplets(:)));
end
offset = P*N_c;
Kf = getKmat2D(rbase_in_f,0);
for k = 1:numel(active_particles)
    pid = active_particles(k);
    seg = (pid-1)*3+1:3*pid;
    fine_idx = offset + (1:N_f);
    offset = offset + N_f;
    FT(seg) = FT(seg) + Kf'*[tau_stokes_x(fine_idx); tau_stokes_y(fine_idx)];

    if ~isempty(rimage_k{pid})
        n_extra = numel(rimage_k{pid});
        extra_idx = offset + (1:n_extra);
        offset = offset + n_extra;
        Kim = getKmat2D(rimage_k{pid},q(pid));
        FT(seg) = FT(seg) + Kim'*[tau_stokes_x(extra_idx); tau_stokes_y(extra_idx)];
    end
end

if offset ~= numel(tau_stokes_x)
    error('solve_res_3B_enhanced:IndexMismatch', ...
        'Unexpected fine-source indexing mismatch during postprocessing.');
end

ftest_b = getVelocityField(rvec_in,rcheck_b,tau_stokes_x,tau_stokes_y,use_fmm);
fbound_x = ftest_b(1:numel(rcheck_b));
fbound_y = ftest_b(numel(rcheck_b)+1:end);

fb_x = [];
fb_y = [];
for k = 1:P
    rhs_f = @(x) [U(k,1)-W(k)*(imag(x)-imag(q(k))); ...
                  U(k,2)+W(k)*(real(x-q(k)))];
    fb_true = rhs_f(rcheck_b((k-1)*n_bound+1:k*n_bound));
    fb_x = [fb_x; fb_true(1:n_bound)];
    fb_y = [fb_y; fb_true(n_bound+1:end)];
end

rel_res = max(hypot(fb_x-fbound_x,fb_y-fbound_y)) / max(hypot(fb_x,fb_y));
fprintf('Max surface rel residual at new nodes is %.3e\n', rel_res);

if visualise_sol
    figure();
    semilogy(abs(lambda));
    title('Triplet corr resistance');
    axis tight
end

sol = struct();
sol.lambda = lambda;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.rel_res = rel_res;
sol.resvec = resvec;
sol.triplet_diagnostics = triplet_diag;

end

function test_solve_res_3B
close all;
rng(9);

delta = 1e-3;
q = hexagonal_lattice(delta,1,1);
P = numel(q);
U = randn(P,2);
W = randn(P,1);

opt = get2Dparams(P,150,150);
opt.delta_pair = 0.2;
opt.gmres_tol = 1e-8;
opt.visualise_sol = 1;
opt.visualise_grid = 0;
opt.debug = 0;
opt.gmres_verbose = 0;
opt.show_counter = 0;
opt.use_fmm = false;
opt.reuse_pair_basis_by_sep = true;
opt.reuse_triplet_basis_by_shape = true;
opt.shared_triplet_tol = 1e-4;
opt.shared_sep_tol = 1e-4;
opt.triplet_basis_debug = 0;
opt.triplet_accumulation_mode = 'sum';
opt.N_peanut = 0;
opt.cmap = 0;

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
triplets = getCloseTripletsFromPairs(pairs,P);
t_c = linspace(0,2*pi,opt.N_c+1)';
t_c = t_c(1:end-1);
rbase_in_c = opt.Rp_c*cos(t_c)+1i*opt.Rp_c*sin(t_c);
tout_c = linspace(0,2*pi,ceil(opt.a_c*opt.N_c)+1)';
tout_c = tout_c(1:end-1);
rbase_out_c = cos(tout_c)+1i*sin(tout_c);
t_f = linspace(0,2*pi,opt.N_f+1)';
t_f = t_f(1:end-1);
rbase_in_f = opt.Rp_f*cos(t_f)+1i*opt.Rp_f*sin(t_f);

fprintf('=== solve_res_3B_enhanced self-test ===\n');
fprintf('P=%d, close pairs=%d, close triplets=%d\n', ...
    P,size(pairs,1),size(triplets,1));
disp('pairs ='); disp(pairs);
disp('triplets ='); disp(triplets);
if isempty(triplets)
    fprintf('WARNING: no close triplets were detected.\n');
end

pair_lookup = getPairRowLookup(pairs,P);
opt_debug = opt;
opt_debug.triplet_basis_debug = 1;
opt_debug.triplet_lsq_test_mode = 'compare_all';
opt_debug.pair_basis_debug = 0;
[~,~,~,~,~,~,pair_cache_debug] = ...
    getPairBasisStokes(q(:),rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_debug,[],rbase_out_c);
pair_test_data = struct('pairs',pairs,'pair_lookup',pair_lookup,'pair_cache',pair_cache_debug);
triplet_cache_debug = getTripletBasisStokes(q(:),rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    triplets,opt_debug,[],pair_test_data);
fprintf('debug compare_all: cached pair groups=%d, cached triplet groups=%d\n', ...
    pair_cache_debug.n_groups,triplet_cache_debug.n_groups);

for mode = {'edge_weighted','pair_subtracted'}
    opt_mode = opt;
    opt_mode.triplet_accumulation_mode = mode{1};
    if strcmp(mode{1},'pair_subtracted')
        opt_mode.triplet_lsq_mode = 'full';
    else
        opt_mode.triplet_lsq_mode = 'rhs_subtracted_pair_complement';
    end

    [~,~,~,~,~,~,pair_cache] = ...
        getPairBasisStokes(q(:),rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_mode,[],rbase_out_c);
    pair_groups = pair_cache.n_groups;
    if strcmp(opt_mode.triplet_lsq_mode,'full')
        pair_test_data_mode = [];
    else
        pair_test_data_mode = struct('pairs',pairs,'pair_lookup',pair_lookup,'pair_cache',pair_cache);
    end
    triplet_cache = getTripletBasisStokes(q(:),rbase_in_c,rbase_in_f,rimage_vec,refine, ...
        triplets,opt_mode,[],pair_test_data_mode);
    fprintf('\nmode: %s, triplet_lsq_mode: %s\n', mode{1},opt_mode.triplet_lsq_mode);
    fprintf('cached pair groups=%d, cached triplet groups=%d\n', ...
        pair_groups,triplet_cache.n_groups);

    [FT,sol] = solve_res_3B_enhanced(q,U,W,opt_mode);
    fprintf('resistance: it=%d, rel_res=%.3e, ||lambda||_2=%.3e, ||lambda||_inf=%.3e\n', ...
        sol.it,sol.rel_res,norm(sol.lambda),max(abs(sol.lambda)));

    if strcmp(mode{1},'pair_subtracted')
        print_triplet_diag(sol.triplet_diagnostics);
    end

    [F,T] = unpackFT_local(FT);
    [UW_back,mob_sol] = solve_mob_3B_enhanced(q,F,T,opt_mode);
    refUW = packUW_local(U,W);
    two_way = norm(UW_back-refUW,inf)/max(1,norm(refUW,inf));
    fprintf('two-way res->mob error: %.3e (mob it=%d, mob rel_res=%.3e)\n', ...
        two_way,mob_sol.it,mob_sol.rel_res);
end
end

function uw = packUW_local(U,W)
P = size(U,1);
uw = zeros(3*P,1);
uw(1:3:end) = U(:,1);
uw(2:3:end) = U(:,2);
uw(3:3:end) = W;
end

function [F,T] = unpackFT_local(FT)
P = numel(FT)/3;
F = [FT(1:3:end), FT(2:3:end)];
T = FT(3:3:end);
end

function print_triplet_diag(diag_info)
if isempty(diag_info) || ~isfield(diag_info,'used_pair_subtraction') || ~diag_info.used_pair_subtraction
    fprintf('pair-subtracted diagnostics: unavailable\n');
    return
end

fprintf(['pair-subtracted triplets: representative row %d (%d,%d,%d), ', ...
    '||beta_full||=%.3e, ||sum pair||=%.3e, ||beta_3only||=%.3e\n'], ...
    diag_info.rep_row,diag_info.rep_ids(1),diag_info.rep_ids(2),diag_info.rep_ids(3), ...
    diag_info.rep_norm_full,diag_info.rep_norm_pair_sum,diag_info.rep_norm_remainder);
fprintf('pair-subtracted triplets: max ||beta_full||=%.3e, max ||sum pair||=%.3e, max ||beta_3only||=%.3e\n', ...
    diag_info.max_norm_full,diag_info.max_norm_pair_sum,diag_info.max_norm_remainder);
fprintf('pair-subtracted triplets: max coarse-cross cancellation residual=%.3e\n', ...
    diag_info.max_cross_cancel_res);
end
