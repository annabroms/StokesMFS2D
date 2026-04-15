function [UW,sol] = solve_mob_3B_enhanced(q,F,T,opt)
%SOLVE_MOB_3B_ENHANCED Solve a 2D Stokes mobility problem with triplet corrections.

if nargin == 0
    test_solve_mob_3B;
    return
end

if nargin < 4 || ~isstruct(opt)
    error('solve_mob_3B_enhanced requires q, F, T, and an options struct opt.');
end

q = q(:);
T = T(:);
P = numel(q);
assert(size(F,1)==P,'F must have one row per particle.');
assert(size(F,2)==2,'F must have two columns [Fx, Fy].');
assert(numel(T)==P,'T must have one entry per particle.');

visualise_sol = logical(getOptField(opt,'visualise_sol',getOptField(opt,'visualise',0)));
gmres_tol = getOptField(opt,'gmres_tol',1e-10);
debug = logical(getOptField(opt,'debug',false));
surface_error_mode = lower(char(getOptField(opt,'surface_error_mode','abs')));
gmres_verbose = getOptField(opt,'gmres_verbose',0);
maxit = getOptField(opt,'maxit',800);
use_fmm = logical(getOptField(opt,'use_fmm',true));
triplet_mode = lower(char(getOptField(opt,'triplet_accumulation_mode','edge_weighted')));
triplet_lsq_mode = lower(char(getOptField(opt,'triplet_lsq_mode','full')));

if ~any(strcmp(surface_error_mode,{'abs','rel'}))
    error('surface_error_mode must be ''abs'' or ''rel''.');
end

fprintf('==== START: %s ====\n', mfilename);

opt.P = P;
opt.rad = 1;
opt.project_force = true;
opt.gmres_verbose = gmres_verbose;
opt.use_fmm = use_fmm;
opt.N_peanut = 0;
opt.cmap = 0;

N_c = getOptField(opt,'N_c',150);
N_f = getOptField(opt,'N_f',150);
a_c = getOptField(opt,'a_c',1.2);
tol_c = getOptField(opt,'tol_c',1e-12);
sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c);
Rp_c = getOptField(opt,'Rp_c',max([1-sep_c,0.01]));
Rp_f = getOptField(opt,'Rp_f',max([1-sep_f,0.01]));

tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1)';
tout_c = tout_c(1:end-1);
rbase_out_c = cos(tout_c)+1i*sin(tout_c);

tin_c = linspace(0,2*pi,N_c+1)';
tin_c = tin_c(1:end-1);
rbase_in_c = Rp_c*cos(tin_c)+1i*Rp_c*sin(tin_c);

tin_f = linspace(0,2*pi,N_f+1)';
tin_f = tin_f(1:end-1);
rbase_in_f = Rp_f*cos(tin_f)+1i*Rp_f*sin(tin_f);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);
triplets = getCloseTripletsFromPairs(pairs,P);
pair_lookup = getPairRowLookup(pairs,P);
if strcmp(triplet_mode,'pair_subtracted')
    triplet_weights = [];
else
    triplet_weights = getTripletAccumulationWeights(triplets,P,triplet_mode);
end

rin = zeros(P*N_c,1);
rout = zeros(P*numel(rbase_out_c),1);
for k = 1:P
    rin((k-1)*N_c+1:k*N_c) = q(k)+rbase_in_c;
    rout((k-1)*numel(rbase_out_c)+1:k*numel(rbase_out_c)) = q(k)+rbase_out_c;
end

[urhs,~,~] = getRecompletionFlow(rin,rout,q,F,T);
[U,Y,Lc] = getSelfPseudoMobilityStokes(1,q,rbase_in_c,rbase_out_c,[],[0,ceil(a_c*N_c)]);

pair_Upf = [];
pair_Ypf = [];
pair_cache = struct('enabled',false);
pair_test_data = [];
if strcmp(triplet_mode,'pair_subtracted')
    [pair_Upf,pair_Ypf,~,~,~,~,pair_cache] = ...
        getPairBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,Lc{1},rbase_out_c);
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
            q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_pair_trip,Lc{1},rbase_out_c);
        pair_test_data = struct('pairs',pairs,'pair_lookup',pair_lookup,'pair_cache',pair_cache_trip);
    end
end

triplet_cache = getTripletBasisStokes(q,rbase_in_c,rbase_in_f,rimage_vec,refine, ...
    triplets,opt,[],pair_test_data);

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.rimage_vec = rimage_vec;
geom.refine = refine;
geom.opt = opt;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.triplets = triplets;
geom.triplet_weights = triplet_weights;
geom.pair_lookup = pair_lookup;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Lc = Lc{1};
basis.triplet_cache = triplet_cache;
basis.Upf = pair_Upf;
basis.Ypf = pair_Ypf;
basis.pair_cache = pair_cache;

n_bound = 803;
tcheck = linspace(0,2*pi,n_bound)';
rcheck_b = zeros(P*n_bound,1);
for k = 1:P
    rcheck_b((k-1)*n_bound+1:k*n_bound) = q(k) + cos(tcheck) + 1i*sin(tcheck);
end

if debug
    ncols = numel(urhs);
    CC = zeros(ncols,ncols);
    x = zeros(size(urhs));
    fprintf('== Debug mode: building system matrix ==\n');
    for k = 1:ncols
        fprintf('build col nbr: %u/%u\n', k,ncols);
        x(:) = 0;
        x(k) = 1;
        CC(:,k) = matvec_mob_3B_enhanced(x,geom,basis);
    end
    figure();
    imagesc(log10(abs(CC)));
    colorbar;
    title([mfilename ': log_{10} |matvec system matrix|'],'interpreter','none');

    [V,D] = eig(CC);
    D = diag(D); 
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')

end

disp(' == Solving... == ');
[tau,it,resvec,~] = helsing_gmres(@(x) matvec_mob_3B_enhanced(x,geom,basis), ...
    urhs,2*numel(rout),maxit,gmres_tol,opt.gmres_verbose,rout);

if visualise_sol
    figure();
    semilogy(resvec);
    xlabel('iteration number','interpreter','latex');
    ylabel('Estimated relative residual');
    axis tight
    grid on
    title('GMRES convergence mobility, triplet corrections','interpreter','latex');
end

disp(' == Postprocessing == ');
[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx,tau_stokes_nonpy,tau_stokes_e_nonpx,tau_stokes_e_nonpy, ...
    rimage_k,triplet_diag] = getMobTripletTransformationStokes(tau,geom,basis);

lambda = [tau_stokes_x; tau_stokes_y];

Kc = getKmat2D(rbase_in_c,0);
UW = zeros(3*P,1);
for k = 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[tau_stokes_nonpx(coarse_ind{k}); tau_stokes_nonpy(coarse_ind{k})];
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
    fine_idx = offset + (1:N_f);
    offset = offset + N_f;
    UW((pid-1)*3+1:3*pid) = UW((pid-1)*3+1:3*pid) - ...
        Kf'*[tau_stokes_nonpx(fine_idx); tau_stokes_nonpy(fine_idx)];
    if ~isempty(rimage_k{pid})
        Kim = getKmat2D(rimage_k{pid},q(pid));
        UW((pid-1)*3+1:3*pid) = UW((pid-1)*3+1:3*pid) - ...
            Kim'*[tau_stokes_e_nonpx{pid}; tau_stokes_e_nonpy{pid}];
    end
end

if offset ~= numel(tau_stokes_nonpx)
    error('solve_mob_3B_enhanced:IndexMismatch', ...
        'Unexpected fine-source indexing mismatch during postprocessing.');
end

B = getKmat2D(rcheck_b(1:n_bound)-q(1),0);
u_lhs = zeros(2*P*n_bound,1);
for k = 1:P
    resk = B*UW(3*(k-1)+1:3*k);
    seg = (k-1)*n_bound+1:k*n_bound;
    u_lhs(seg) = resk(1:n_bound);
    u_lhs(P*n_bound + seg) = resk(n_bound+1:end);
end

u_rhs = getVelocityField(rvec_in,rcheck_b,tau_stokes_x,tau_stokes_y,use_fmm);
S_0 = getRecompletionFlow(rin,rcheck_b,q,F,T);
u_rhs = u_rhs - S_0;

diff_vec = u_rhs-u_lhs;
res = abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end));
abs_res = max(res);
max_abs = max(abs(S_0(1:end/2)+1i*S_0(end/2+1:end)));
if max_abs > 0
    err_vec = res/max_abs;
else
    err_vec = res;
end
rel_res = max(err_vec);
fprintf('Relative boundary error: %.3e\n', rel_res);
fprintf('Absolute boundary error: %.3e\n', abs_res);

if visualise_sol
    figure();
    semilogy(abs(lambda));
    title('Source strengths mobility with triplet corr');
    axis tight
end

sol = struct();
sol.lambda = lambda;
sol.it = it;
sol.gmres_tol = gmres_tol;
sol.rel_res = rel_res;
sol.abs_res = abs_res;
sol.resvec = resvec;
sol.triplet_diagnostics = triplet_diag;

end

function test_solve_mob_3B
close all;
rng(8);

delta = 1e-3;
q = hexagonal_lattice(delta,1,1);
q = q([4 6 7]); 
P = numel(q);
F = randn(P,2);
T = randn(P,1);

opt = get2Dparams(P,60,60);
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
[~,~,Lc_test] = getSelfPseudoMobilityStokes(1,q(:),rbase_in_c,rbase_out_c,[],[0,ceil(opt.a_c*opt.N_c)]);

fprintf('=== solve_mob_3B_enhanced self-test ===\n');
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

% [~,~,~,~,~,~,pair_cache] = ...
%     getPairBasisStokes(q(:),rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_debug,Lc_test{1},rbase_out_c);
% pair_test_data = struct('pairs',pairs,'pair_lookup',pair_lookup,'pair_cache',pair_cache);
% triplet_cache_debug = getTripletBasisStokes(q(:),rbase_in_c,rbase_in_f,rimage_vec,refine, ...
%     triplets,opt_debug,[],pair_test_data);
% fprintf('debug compare_all: cached pair groups=%d, cached triplet groups=%d\n', ...
%     pair_cache.n_groups,triplet_cache_debug.n_groups);


mode = 'pair_subtracted';
%mode = 'sum';
opt_mode = opt;
opt_mode.triplet_accumulation_mode = mode;
opt_mode.triplet_lsq_mode = 'full';
%opt_mode.triplet_lsq_mode = 'rhs_subtracted_pair_complement';

fprintf('\nmode: %s, triplet_lsq_mode: %s\n', mode,opt_mode.triplet_lsq_mode);

[UW,sol] = solve_mob_3B_enhanced(q,F,T,opt_mode);
fprintf('mobility: it=%d, rel_res=%.3e, abs_res=%.3e, ||lambda||_2=%.3e, ||lambda||_inf=%.3e\n', ...
    sol.it,sol.rel_res,sol.abs_res,norm(sol.lambda),max(abs(sol.lambda)));

[UW2,sol2] = solve_mob_2B_enhanced(q,F,T,opt_mode);
fprintf('mobility: it=%d, rel_res=%.3e, abs_res=%.3e, ||lambda||_2=%.3e, ||lambda||_inf=%.3e\n', ...
    sol2.it,sol2.rel_res,sol2.abs_res,norm(sol2.lambda),max(abs(sol2.lambda)));

if strcmp(mode,'pair_subtracted')
    print_triplet_diag(sol.triplet_diagnostics);
end

[U,W] = unpackUW_local(UW);
[FT_back,res_sol] = solve_res_3B_enhanced(q,U,W,opt_mode);
refFT = packFT_local(F,T);
two_way = norm(FT_back-refFT,inf)/max(1,norm(refFT,inf));
fprintf('two-way mob->res error: %.3e (res it=%d, res rel_res=%.3e)\n', ...
    two_way,res_sol.it,res_sol.rel_res);

end

function ft = packFT_local(F,T)
    P = size(F,1);
    ft = zeros(3*P,1);
    ft(1:3:end) = F(:,1);
    ft(2:3:end) = F(:,2);
    ft(3:3:end) = T;
end

function [U,W] = unpackUW_local(UW)
    U = [UW(1:3:end), UW(2:3:end)];
    W = UW(3:3:end);
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
end
