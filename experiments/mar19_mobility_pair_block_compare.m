clear;
close all;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Mobility Pair/Peanut Projection Benchmark (Apr 24, 2026) ===\n');

cfg = struct();
cfg.rad = 1;
cfg.delta_pair = 0.2;
cfg.rel_tol = 1e-13;
cfg.timeit_repeats = 3;
cfg.gaps = [0.02 0.05 0.10 0.15];
cfg.case_specs = [60 150 400; 80 60 400];

results = repmat(init_result_row(),0,1);

for case_id = 1:size(cfg.case_specs,1)
    N_c = cfg.case_specs(case_id,1);
    N_f = cfg.case_specs(case_id,2);
    N_peanut = cfg.case_specs(case_id,3);

    fprintf('\nCase %d: N_c=%d, N_f=%d, N_peanut=%d\n', ...
        case_id,N_c,N_f,N_peanut);

    for gap_id = 1:numel(cfg.gaps)
        gap = cfg.gaps(gap_id);
        case_data = build_case_data(cfg.rad,N_c,N_f,N_peanut,gap,cfg.delta_pair);
        result = benchmark_case(case_data,cfg);

        results(end+1,1) = result; %#ok<SAGROW>
        print_case_result(result);
    end
end

fprintf('\n=== Summary Table ===\n');
fprintf(['%4s %4s %8s %8s %6s %6s %10s %10s %10s %10s %10s %10s\n'], ...
    'N_c','N_f','N_peanut','gap','nsrc','ntar', ...
    'stok x','pair x','dense x','peanut x','dense x','map err');
for k = 1:numel(results)
    result = results(k);
    fprintf(['%4d %4d %8d %8.3f %6d %6d %10.2f %10.2f %10.2f ', ...
        '%10.2f %10.2f %10.2e\n'], ...
        result.N_c,result.N_f,result.N_peanut,result.gap, ...
        result.nsrc,result.ntar,result.stok_speedup, ...
        result.pair_speedup,result.pair_dense_speedup, ...
        result.peanut_speedup,result.peanut_dense_speedup, ...
        result.peanut_map_relerr);
end

function result = benchmark_case(case_data,cfg)
warmup_case(case_data);

S_new = case_data.S_new;
S_old = case_data.S_old;

Atot_old = assemble_pair_legacy(S_new,case_data.Lf_pair,case_data.Lr_pair);
Atot_dense = assemble_pair_dense_projection(S_new,case_data.Ppair,case_data.Lr_pair);
Atot_new = assemble_pair_lowrank(S_new,case_data.pair_moment_map, ...
    case_data.pair_rbm_map,case_data.pair_moment_gram, ...
    case_data.pair_target_rbm_map);

Ntot_old = project_peanut_legacy(case_data.Nf_peanut,case_data.Lf_pair);
Ntot_dense = project_peanut_dense(case_data.Nf_peanut,case_data.Ppair);
Ntot_new = project_peanut_lowrank(case_data.Nf_peanut, ...
    case_data.pair_moment_map,case_data.pair_rbm_map, ...
    case_data.pair_moment_gram);

assert_relerr(relerr(S_new,S_old),cfg.rel_tol, ...
    'stokSLPmat new/legacy mismatch');
assert_relerr(relerr(Atot_new,Atot_old),cfg.rel_tol, ...
    'low-rank pair assembly mismatch');
assert_relerr(relerr(Ntot_new,Ntot_old),cfg.rel_tol, ...
    'low-rank peanut projection mismatch');
assert_relerr(relerr(Atot_dense,Atot_old),cfg.rel_tol, ...
    'dense projected pair assembly mismatch');
assert_relerr(relerr(Ntot_dense,Ntot_old),cfg.rel_tol, ...
    'dense projected peanut projection mismatch');

[~,~,Atot_prod] = getPairBlockStokes(case_data.rin_pair_f,case_data.rout_f, ...
    case_data.pair_moment_map,case_data.pair_rbm_map, ...
    case_data.pair_moment_gram,case_data.pair_target_rbm_map);
assert_relerr(relerr(Atot_prod,Atot_new),cfg.rel_tol, ...
    'getPairBlockStokes returned the wrong pair operator');

[DC_prod,Y_prod] = getPeanutBlockStokes(case_data.rin_pair_c, ...
    case_data.rin_pair_f,case_data.rout_peanut,case_data.Lc_pair, ...
    case_data.pair_moment_map,case_data.pair_rbm_map, ...
    case_data.pair_moment_gram);
[Y_ref,U_ref] = getPseudoFactors(case_data.Npeanut_proj,1e-14,0);
peanut_map_prod = Y_prod*DC_prod;
peanut_map_ref = Y_ref*(U_ref'*Ntot_new);
peanut_map_relerr = relerr(peanut_map_prod,peanut_map_ref);
assert_relerr(peanut_map_relerr,cfg.rel_tol, ...
    'getPeanutBlockStokes returned the wrong peanut map');

result = init_result_row();
result.N_c = case_data.N_c;
result.N_f = case_data.N_f;
result.N_peanut = case_data.N_peanut;
result.gap = case_data.gap;
result.nsrc = case_data.nsrc;
result.ntar = case_data.ntar;
result.relerr_stok = relerr(S_new,S_old);
result.relerr_pair = relerr(Atot_new,Atot_old);
result.relerr_pair_dense = relerr(Atot_dense,Atot_old);
result.relerr_peanut = relerr(Ntot_new,Ntot_old);
result.relerr_peanut_dense = relerr(Ntot_dense,Ntot_old);
result.peanut_map_relerr = peanut_map_relerr;

result.t_stok_old = median_timeit(@() stokSLPmat_legacy( ...
    case_data.rin_pair_f,case_data.rout_f,1),cfg.timeit_repeats);
result.t_stok_new = median_timeit(@() stokSLPmat( ...
    case_data.rin_pair_f,case_data.rout_f,1),cfg.timeit_repeats);

result.t_pair_old = median_timeit(@() assemble_pair_legacy( ...
    S_new,case_data.Lf_pair,case_data.Lr_pair),cfg.timeit_repeats);
result.t_pair_dense = median_timeit(@() assemble_pair_dense_projection( ...
    S_new,case_data.Ppair,case_data.Lr_pair),cfg.timeit_repeats);
result.t_pair_new = median_timeit(@() assemble_pair_lowrank( ...
    S_new,case_data.pair_moment_map,case_data.pair_rbm_map, ...
    case_data.pair_moment_gram,case_data.pair_target_rbm_map), ...
    cfg.timeit_repeats);

result.t_peanut_old = median_timeit(@() project_peanut_legacy( ...
    case_data.Nf_peanut,case_data.Lf_pair),cfg.timeit_repeats);
result.t_peanut_dense = median_timeit(@() project_peanut_dense( ...
    case_data.Nf_peanut,case_data.Ppair),cfg.timeit_repeats);
result.t_peanut_new = median_timeit(@() project_peanut_lowrank( ...
    case_data.Nf_peanut,case_data.pair_moment_map, ...
    case_data.pair_rbm_map,case_data.pair_moment_gram), ...
    cfg.timeit_repeats);

result.stok_speedup = result.t_stok_old/max(result.t_stok_new,eps);
result.pair_speedup = result.t_pair_old/max(result.t_pair_new,eps);
result.pair_dense_speedup = result.t_pair_dense/max(result.t_pair_new,eps);
result.peanut_speedup = result.t_peanut_old/max(result.t_peanut_new,eps);
result.peanut_dense_speedup = result.t_peanut_dense/max(result.t_peanut_new,eps);
end

function case_data = build_case_data(rad,N_c,N_f,N_peanut,gap,delta_pair)
q_pair = [0; 2*rad + gap];

opt = get2Dparams(2,N_c,N_f);
opt.rad = rad;
opt.delta_pair = delta_pair;
opt.N_peanut = N_peanut;

[~,~,~,rimage_pairs,refine,pairs] = getEnhancedGrid(q_pair,opt);
assert(size(pairs,1) == 1 && all(pairs(1,:) == [1 2]), ...
    'Expected exactly one close pair in the benchmark geometry.');

rbase_in_c = build_circle_nodes(opt.Rp_c,N_c);
rbase_in_f = build_circle_nodes(opt.Rp_f,N_f);
rout_base_f = build_circle_nodes(rad,ceil(opt.a_f*N_f));

refine_i = refine{1,2};
refine_j = refine{2,1};
rimage_i = rimage_pairs{1,2};
rimage_j = rimage_pairs{2,1};

rin_pair_f = [q_pair(1)+rbase_in_f; rimage_i; q_pair(2)+rbase_in_f; rimage_j];
rout_f = [q_pair(1)+rout_base_f; refine_i; q_pair(2)+rout_base_f; refine_j];
rin_pair_c = [q_pair(1)+rbase_in_c; q_pair(2)+rbase_in_c];
rout_peanut = createPeanut(q_pair(1),q_pair(2),N_peanut,0,rad);

assert(mod(numel(rin_pair_f),2) == 0,'Expected symmetric fine-source layout.');
assert(mod(numel(rout_f),2) == 0,'Expected symmetric pair-target layout.');

Kf1 = getKmat2D(rin_pair_f(1:end/2),q_pair(1));
Kf2 = getKmat2D(rin_pair_f(end/2+1:end),q_pair(2));
B1 = getKmat2D(rout_f(1:end/2),q_pair(1));
B2 = getKmat2D(rout_f(end/2+1:end),q_pair(2));
pair_moment_map = getKftPair(Kf1,Kf2);
pair_rbm_map = pair_moment_map';
pair_moment_gram = pair_moment_map*pair_rbm_map;
pair_target_rbm_map = getKftPair(B1,B2)';
Lf_pair = getLfPair(Kf1,Kf2);
Lr_pair = getLrPair(B1,B2,Kf1,Kf2);
Ppair = eye(size(Lf_pair)) - Lf_pair;

Kc = getKmat2D(rbase_in_c,0);
Lc_body = Kc*((Kc'*Kc)\Kc');
Lc_pair = getILpair(Lc_body);

S_old = stokSLPmat_legacy(rin_pair_f,rout_f,1);
S_new = stokSLPmat(rin_pair_f,rout_f,1);
Nf_peanut = stokSLPmat(rin_pair_f,rout_peanut,1);
Npeanut = stokSLPmat(rin_pair_c,rout_peanut,1);

case_data = struct();
case_data.N_c = N_c;
case_data.N_f = N_f;
case_data.N_peanut = N_peanut;
case_data.gap = gap;
case_data.q_pair = q_pair;
case_data.rin_pair_f = rin_pair_f;
case_data.rin_pair_c = rin_pair_c;
case_data.rout_f = rout_f;
case_data.rout_peanut = rout_peanut;
case_data.Lf_pair = Lf_pair;
case_data.Lr_pair = Lr_pair;
case_data.Ppair = Ppair;
case_data.Lc_pair = Lc_pair;
case_data.pair_moment_map = pair_moment_map;
case_data.pair_rbm_map = pair_rbm_map;
case_data.pair_moment_gram = pair_moment_gram;
case_data.pair_target_rbm_map = pair_target_rbm_map;
case_data.S_old = S_old;
case_data.S_new = S_new;
case_data.Nf_peanut = Nf_peanut;
case_data.Npeanut_proj = Npeanut*Lc_pair;
case_data.nsrc = numel(rin_pair_f)/2;
case_data.ntar = numel(rout_f)/2;
end

function nodes = build_circle_nodes(radius,N)
t = linspace(0,2*pi,N+1)';
t = t(1:end-1);
nodes = radius*(cos(t) + 1i*sin(t));
end

function warmup_case(case_data)
stokSLPmat_legacy(case_data.rin_pair_f,case_data.rout_f,1);
stokSLPmat(case_data.rin_pair_f,case_data.rout_f,1);
assemble_pair_legacy(case_data.S_new,case_data.Lf_pair,case_data.Lr_pair);
assemble_pair_dense_projection(case_data.S_new,case_data.Ppair,case_data.Lr_pair);
assemble_pair_lowrank(case_data.S_new,case_data.pair_moment_map, ...
    case_data.pair_rbm_map,case_data.pair_moment_gram, ...
    case_data.pair_target_rbm_map);
project_peanut_legacy(case_data.Nf_peanut,case_data.Lf_pair);
project_peanut_dense(case_data.Nf_peanut,case_data.Ppair);
project_peanut_lowrank(case_data.Nf_peanut,case_data.pair_moment_map, ...
    case_data.pair_rbm_map,case_data.pair_moment_gram);
end

function t = median_timeit(fun,n_repeats)
timings = zeros(n_repeats,1);
for k = 1:n_repeats
    timings(k) = timeit(fun);
end
t = median(timings);
end

function S = stokSLPmat_legacy(rin,rout,mu)
r = bsxfun(@minus, rout, rin.');
irr = 1./(conj(r).*r);
d1 = real(r);
d2 = imag(r);
c = 1/(4*pi*mu);
logir = -log(abs(r));
A12 = d1.*d2.*irr;
S = c*[logir + d1.^2.*irr, A12; ...
       A12,                logir + d2.^2.*irr];
end

function Atot = assemble_pair_legacy(S,Lf_pair,Lr_pair)
Atot = S - S*Lf_pair + Lr_pair;
end

function Atot = assemble_pair_dense_projection(S,Ppair,Lr_pair)
Atot = S*Ppair + Lr_pair;
end

function Atot = assemble_pair_lowrank(S,pair_moment_map,pair_rbm_map, ...
    pair_moment_gram,pair_target_rbm_map)
Atot = S + (pair_target_rbm_map - (S*pair_rbm_map)/pair_moment_gram) * ...
    pair_moment_map;
end

function Ntot = project_peanut_legacy(Nf,Lf_pair)
Ntot = Nf - Nf*Lf_pair;
end

function Ntot = project_peanut_dense(Nf,Ppair)
Ntot = Nf*Ppair;
end

function Ntot = project_peanut_lowrank(Nf,pair_moment_map,pair_rbm_map, ...
    pair_moment_gram)
Ntot = Nf - (Nf*pair_rbm_map)/pair_moment_gram * pair_moment_map;
end

function e = relerr(A,B)
e = norm(A-B,'fro') / max(1,norm(B,'fro'));
end

function assert_relerr(err_value,tol,msg)
assert(err_value < tol,'%s (relerr %.3e, tol %.3e).',msg,err_value,tol);
end

function print_case_result(result)
fprintf(['  gap=%.3f  nsrc=%d  ntar=%d  rel(S)=%.3e  ', ...
    'rel(pair)=%.3e  rel(peanut)=%.3e\n'], ...
    result.gap,result.nsrc,result.ntar,result.relerr_stok, ...
    result.relerr_pair,result.relerr_peanut);
fprintf(['    stok old/new      : %.4gs / %.4gs  speedup %.2fx\n'], ...
    result.t_stok_old,result.t_stok_new,result.stok_speedup);
fprintf(['    pair old/dense/new: %.4gs / %.4gs / %.4gs  ', ...
    'new vs old %.2fx  new vs dense %.2fx\n'], ...
    result.t_pair_old,result.t_pair_dense,result.t_pair_new, ...
    result.pair_speedup,result.pair_dense_speedup);
fprintf(['    peanut old/dense/new: %.4gs / %.4gs / %.4gs  ', ...
    'new vs old %.2fx  new vs dense %.2fx  map rel %.3e\n'], ...
    result.t_peanut_old,result.t_peanut_dense,result.t_peanut_new, ...
    result.peanut_speedup,result.peanut_dense_speedup, ...
    result.peanut_map_relerr);
end

function result = init_result_row()
result = struct( ...
    'N_c',0,'N_f',0,'N_peanut',0,'gap',0, ...
    'nsrc',0,'ntar',0, ...
    'relerr_stok',0,'relerr_pair',0,'relerr_pair_dense',0, ...
    'relerr_peanut',0,'relerr_peanut_dense',0, ...
    'peanut_map_relerr',0, ...
    't_stok_old',0,'t_stok_new',0, ...
    't_pair_old',0,'t_pair_dense',0,'t_pair_new',0, ...
    't_peanut_old',0,'t_peanut_dense',0,'t_peanut_new',0, ...
    'stok_speedup',0,'pair_speedup',0,'pair_dense_speedup',0, ...
    'peanut_speedup',0,'peanut_dense_speedup',0);
end
