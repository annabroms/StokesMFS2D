function report = test_laplace_peanut_big_sparse()
%TEST_LAPLACE_PEANUT_BIG_SPARSE Compare serial and big-sparse matvecs.

fprintf('test_laplace_peanut_big_sparse: use_fmm=true matvec comparison\n');
rng(220422);

use_fmm = true;
problems = {'capacitance','elastance'};
geometries = {'axis_shared','rotated_shared'};
build_modes = {'precomputed','streaming'};
reuse_modes = [false true];

case_name = {};
relerr = [];
tol_used = [];

for ip = 1:numel(problems)
    for ig = 1:numel(geometries)
        for ir = 1:numel(reuse_modes)
            ref = buildLaplaceMatvecCase(problems{ip},geometries{ig}, ...
                reuse_modes(ir),'serial',use_fmm);
            tau = randn(numel(ref.geom.rvec_out),1);
            res_ref = matvec_lap_peanut_enhanced(tau,ref.geom,ref.basis);

            for im = 1:numel(build_modes)
                big = buildLaplaceMatvecCase(problems{ip},geometries{ig}, ...
                    reuse_modes(ir),build_modes{im},use_fmm);
                res_big = matvec_lap_peanut_big_sparse(tau,big.geom,big.basis);
                err = relativeError(res_big,res_ref);
                tol = matvecTolerance(reuse_modes(ir));

                label = sprintf('%s %s reuse=%d %s',problems{ip}, ...
                    geometries{ig},reuse_modes(ir),build_modes{im});
                fprintf('  %-49s relerr %.3e (tol %.3e)\n',label,err,tol);
                assert(err <= tol, ...
                    'test_laplace_peanut_big_sparse:MatvecMismatch', ...
                    'Matvec mismatch for %s: %.3e exceeds %.3e.',label,err,tol);

                case_name{end+1,1} = label; %#ok<AGROW>
                relerr(end+1,1) = err; %#ok<AGROW>
                tol_used(end+1,1) = tol; %#ok<AGROW>
            end
        end
    end
end

report = table(case_name,relerr,tol_used, ...
    'VariableNames',{'case','relerr','tolerance'});
fprintf('test_laplace_peanut_big_sparse: passed\n');
end

function ctx = buildLaplaceMatvecCase(problem,geometry,reuse_pair_basis, ...
    mode,use_fmm)
R = 2;
q = buildLaplaceTestGeometry(geometry,R);
P = numel(q);
N_c = 20;
N_f = 28;

opt = getLaplace2Dparams(P,R,N_c,N_f);
opt.N_peanut = 64;
opt.delta_pair = 0.2*R;
opt.cmap = 1;
opt.compress_cmap = 0;
opt.reuse_pair_basis_by_sep = reuse_pair_basis;
opt.shared_sep_tol = 1e-8;
opt.show_counter = 0;
opt.visualise_sol = 0;
opt.visualise_grid = 0;
opt.gmres_verbose = 0;
opt.get_bndry_field = 0;
opt.use_fmm = use_fmm;
opt.project_charge = strcmp(problem,'elastance');
opt.lap_big_sparse_build_mode = mode;
opt.use_big_sparse = ~strcmp(mode,'serial');

[rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] = ...
    buildLaplaceCircleGrids(R,opt);

rvec_in_c = zeros(P*N_c,1);
rout = zeros(P*nout,1);
for k = 1:P
    rvec_in_c((k-1)*N_c+1:k*N_c) = q(k)+rbase_in_c;
    rout((k-1)*nout+1:k*nout) = q(k)+rbase_out_c;
end

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt);

if strcmp(mode,'streaming')
    UB_all = [];
    YB_all = [];
    UC_all = [];
    YC_all = [];
    Cmap = [];
    Cmap_QV = [];
    pair_cache = initLaplacePairCache();
else
    [UB_all,YB_all,UC_all,YC_all,Cmap,Cmap_QV,pair_cache] = ...
        getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f, ...
        rbase_out_c,rimage_vec,refine,pairs,opt);
end

if opt.project_charge
    [U,Y] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout],true);
else
    [U,Y] = getSelfPseudoLaplace(1,rbase_in_c,rbase_out_c,[0 nout]);
end

geom = struct();
geom.rbase_in_c = rbase_in_c;
geom.rbase_in_f = rbase_in_f;
geom.rout_base_f = rout_base_f;
geom.refine = refine;
geom.opt = opt;
geom.rvec_out = rout;
geom.rcheck = rout;
geom.q = q;
geom.pairs = pairs;
geom.rimage_vec = rimage_vec;
geom.rvec_in = rvec_in_c;
geom.pair_cache = pair_cache;

basis = struct();
basis.U = U;
basis.Y = Y;
basis.Upf = UB_all;
basis.Ypf = YB_all;
basis.DC_all = UC_all;
basis.YC_all = YC_all;
basis.Cmap = Cmap;
basis.Cmap_QV = Cmap_QV;
basis.pair_cache = pair_cache;
basis.Nii = lapSLPmat(rbase_in_c,rbase_out_c);

if opt.use_big_sparse
    [big_sparse,~,pair_cache] = buildLaplacePeanutBigSparse(geom,basis);
    basis.big_sparse = big_sparse;
    basis.pair_cache = pair_cache;
    geom.pair_cache = pair_cache;
end

ctx = struct('geom',geom,'basis',basis);
end

function q = buildLaplaceTestGeometry(geometry,R)
gap = 2e-2;
sep = 2*R + gap;
switch geometry
    case 'axis_shared'
        q = [0; sep; 2*sep];
    case 'rotated_shared'
        theta = 1.3;
        q = [0; sep; sep*exp(1i*theta)];
    otherwise
        error('test_laplace_peanut_big_sparse:BadGeometry', ...
            'Unknown geometry "%s".',geometry);
end
q = q(:);
end

function [rbase_in_c,rbase_out_c,rbase_in_f,rout_base_f,nout] = ...
    buildLaplaceCircleGrids(R,opt)
N_c = opt.N_c;
N_f = opt.N_f;
nout = ceil(opt.a_c*N_c);

t = linspace(0,2*pi,N_c+1)';
t = t(1:end-1);
rbase_in_c = opt.Rp_c*(cos(t)+1i*sin(t));

t = linspace(0,2*pi,nout+1)';
t = t(1:end-1);
rbase_out_c = R*(cos(t)+1i*sin(t));

t = linspace(0,2*pi,N_f+1)';
t = t(1:end-1);
rbase_in_f = opt.Rp_f*(cos(t)+1i*sin(t));

t = linspace(0,2*pi,ceil(opt.a_f*N_f)+1)';
t = t(1:end-1);
rout_base_f = R*(cos(t)+1i*sin(t));
end

function tol = matvecTolerance(reuse_pair_basis)
if reuse_pair_basis
    tol = 2e-6;
else
    tol = 5e-7;
end
end

function e = relativeError(a,b)
e = norm(a-b,inf)/max(1,norm(b,inf));
end
