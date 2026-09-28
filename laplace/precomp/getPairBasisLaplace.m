function [Uf,Yf,Up,Yp,Cmap,Cmap_QV,pair_cache] = getPairBasisLaplace( ...
    q,rbase_in_c,rbase_in_f,rout_base_f,rout_base_c,rimage_vec,refine, ...
    pairs,opt,rbase_in_solver)
%GETPAIRBASISLAPLACE Build pair-basis pseudoinverse factors for Laplace.
%
% Syntax:
%   [Uf,Yf,Up,Yp,Cmap,Cmap_QV,pair_cache] = ...
%       getPairBasisLaplace(q,rbase_in_c,rbase_in_f,rout_base_f,rout_base_c, ...
%       rimage_vec,refine,pairs,opt,rbase_in_solver)
%
% rbase_in_c is the grid on which the square Cmap is constructed.  The
% optional rbase_in_solver is the global interaction grid.  If opt.refit is
% true, charge-preserving per-particle field refits between the two grids
% are stored in pair_cache.  Otherwise the transform uses Fourier transfer.
%
% See also: getPairBlockLaplace, getPeanutBlockLaplace, ...
%   evaluateCoarseOnPairLaplace, getPairTransformationLaplace.
%
% Anna Broms, Mar 2026

q = q(:);
if nargin < 10 || isempty(rbase_in_solver)
    rbase_in_solver = rbase_in_c;
end
rbase_in_c = rbase_in_c(:);
rbase_in_solver = rbase_in_solver(:);
P = numel(q);
N_peanut = opt.N_peanut;
R = opt.rad;
if numel(R) > 1
    R = R(1);
end
reuse_pair_basis = opt.reuse_pair_basis_by_sep;
parallel_precomp = opt.parallel_precomp;
check_rotations = opt.check_rotations;
shared_sep_tol = opt.shared_sep_tol;
project_charge = opt.project_charge;
show_counter = opt.show_counter;
use_pair_map = opt.cmap;
use_refit = use_pair_map && logical(getOptField(opt,'refit',false));
interpolation_mode = getLaplaceInterpolationMode(opt);
N_cmap = getOptField(opt,'N_cmap',numel(rbase_in_c));
if use_pair_map && numel(rbase_in_c) ~= N_cmap
    error('getPairBasisLaplace:BadCmapGridSize', ...
        ['The Cmap source grid has %d nodes, but opt.N_cmap=%d. ', ...
         'Pass the N_cmap-node grid as rbase_in_c.'], ...
        numel(rbase_in_c),N_cmap);
end
if use_pair_map && numel(rbase_in_solver) ~= N_cmap && ...
        (~reuse_pair_basis || ~use_refit)
    error('getPairBasisLaplace:MixedCoarseRequiresReusedRefit', ...
        ['N_cmap may differ from N_c only when ', ...
         'opt.reuse_pair_basis_by_sep=true and opt.refit=true.']);
end

if reuse_pair_basis
    Uf = [];
    Yf = [];
    Up = [];
    Yp = [];
    Cmap = [];
    Cmap_QV = [];
else
    Uf = cell(P);
    Yf = cell(P);

    if N_peanut
        if use_pair_map
            Up = [];
            Yp = [];
            Cmap = cell(P);
            Cmap_QV = cell(P);
        else
            Up = cell(P);
            Yp = cell(P);
            Cmap = [];
            Cmap_QV = [];
        end
    else
        Up = [];
        Yp = [];
        Cmap = [];
        Cmap_QV = [];
    end
end

pair_cache = initLaplacePairCache();
pair_cache.enabled = reuse_pair_basis;
pair_cache.interpolation_mode = interpolation_mode;
pair_cache.check_rotations = check_rotations;
pair_cache.shared_sep_tol = shared_sep_tol;
pair_cache.rout_base_f = rout_base_f(:);
pair_cache.stats.requested_parallel = parallel_precomp;

if isempty(pairs)
    return
end

pair_cache.meta = buildLaplacePairMeta(q,pairs,numel(rbase_in_c), ...
    numel(rbase_in_f),opt);
pair_cache.stats.n_pairs = size(pairs,1);

if use_refit
    if reuse_pair_basis
        rotations = reshape([pair_cache.meta.rot],[],1);
    else
        % Per-pair maps use the same global angular origin as the solver
        % grid, so only the identity rectangular refit is needed.
        rotations = 1;
    end
    transfer_svd_opts = struct( ...
        'use_tikhonov',logical(getOptField(opt,'use_tikhonov',false)), ...
        'tikhonov_tol',getOptField(opt,'tikhonov_tol',[]));
    pair_cache.coarse_transfer = buildLaplaceCoarseGridTransfers( ...
        rbase_in_solver,rbase_in_c,rotations,R, ...
        getOptField(opt,'a_c',1.2),1e-14,transfer_svd_opts);
end

if ~strcmp(interpolation_mode,'none')
    if logical(project_charge)
        error('getPairBasisLaplace:InterpolationCapacitanceOnly', ...
            ['opt.use_interpolation is currently implemented only for ', ...
             'Laplace capacitance (opt.project_charge=false).']);
    end
    if ~reuse_pair_basis || ~use_pair_map || ~logical(N_peanut)
        error('getPairBasisLaplace:InterpolationRequiresCanonicalCmap', ...
            ['Interpolation requires opt.reuse_pair_basis_by_sep=true, ', ...
             'opt.cmap=true, and opt.N_peanut>0.']);
    end
    if check_rotations
        error('getPairBasisLaplace:InterpolationRotationCheckUnsupported', ...
            'opt.check_rotations is not supported with interpolated maps.');
    end
    if ~logical(getOptField(opt,'ellipse_constant',false)) || ...
            getOptField(opt,'Nclust',[]) ~= 150 || ...
            ~logical(getOptField(opt,'use_tikhonov',false)) || ...
            getOptField(opt,'tikhonov_tol',[]) ~= 1e-11 || ...
            logical(getOptField(opt,'compress_cmap',false))
        error('getPairBasisLaplace:InterpolationProfileMismatch', ...
            ['Interpolation requires ellipse_constant=true, Nclust=150, ', ...
             'Tikhonov tolerance 1e-11, and compress_cmap=false.']);
    end
    pair_cache = buildLaplaceInterpolatedPairCache(pair_cache,q,pairs, ...
        rbase_in_c,rbase_in_f,pair_cache.rout_base_f,rbase_in_solver, ...
        rout_base_c,opt);
    if show_counter
        fprintf(['getPairBasisLaplace: evaluated %s interpolation for ', ...
            '%d separation groups covering %d pairs\n'], ...
            interpolation_mode,pair_cache.n_groups,size(pairs,1));
    end
    return
end

if ~reuse_pair_basis
    total_pairs = size(pairs,1);
    pair_cache.stats.branch = 'per_pair';
    use_parallel_pairs = parallel_precomp && total_pairs > 1;
    [use_parallel_pairs,pool_size] = resolve_parallel_precomp(use_parallel_pairs, ...
        'getPairBasisLaplace');
    pair_cache.stats.used_parallel = use_parallel_pairs;
    pair_cache.stats.pool_size = pool_size;

    if use_parallel_pairs
        if show_counter
            fprintf('getPairBasisLaplace: parallel pair build for %d pairs\n', total_pairs);
        end
        pair_entries = cell(total_pairs,1);
        parfor ii = 1:total_pairs
            pair_entries{ii} = buildLaplacePairGroup([],ii, ...
                pair_cache.meta(ii).sep,q,rbase_in_c,rbase_in_f, ...
                rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
                project_charge,false);
        end
        pair_cache.stats.pool_size = get_parallel_pool_size();
        for ii = 1:total_pairs
            [Uf,Yf,Up,Yp,Cmap,Cmap_QV] = assign_pair_entry_outputs( ...
                Uf,Yf,Up,Yp,Cmap,Cmap_QV,pairs,ii,pair_entries{ii},N_peanut,use_pair_map);
        end
        if show_counter
            fprintf('getPairBasisLaplace: finished parallel pair build for %d pairs\n', ...
                total_pairs);
        end
    else
        processed_pairs = 0;
        for ii = 1:total_pairs
            pair_entry = buildLaplacePairGroup([],ii, ...
                pair_cache.meta(ii).sep,q,rbase_in_c,rbase_in_f, ...
                rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
                project_charge,false);
            [Uf,Yf,Up,Yp,Cmap,Cmap_QV] = assign_pair_entry_outputs( ...
                Uf,Yf,Up,Yp,Cmap,Cmap_QV,pairs,ii,pair_entry,N_peanut,use_pair_map);

            processed_pairs = processed_pairs + 1;
            if show_counter
                i = pairs(ii,1);
                p2 = pairs(ii,2);
                fprintf('getPairBasisLaplace: processed pair %d/%d (%d,%d)\n', ...
                    processed_pairs,total_pairs,i,p2);
            end
        end
    end

    if use_pair_map && N_peanut
        i = pairs(1,1);
        j = pairs(1,2);
        report_cmap_size(Cmap{i,j},N_cmap);
    end
    return
end

[group_id,group_sep,rep_rows] = groupLaplacePairSeparations( ...
    pair_cache.meta,shared_sep_tol);
n_groups = numel(group_sep);
pair_cache.n_groups = n_groups;
pair_cache.group_id = group_id;
pair_cache.group_sep = group_sep;
pair_cache.representative_rows = rep_rows;
total_pairs = size(pairs,1);
pair_cache.stats.branch = 'canonical_group';
pair_cache.stats.n_groups = n_groups;
use_parallel_groups = parallel_precomp && n_groups > 1;
[use_parallel_groups,pool_size] = resolve_parallel_precomp(use_parallel_groups, ...
    'getPairBasisLaplace');
pair_cache.stats.used_parallel = use_parallel_groups;
pair_cache.stats.pool_size = pool_size;

if use_parallel_groups
    if show_counter
        fprintf(['getPairBasisLaplace: parallel canonical build for %d groups ', ...
            'covering %d pairs\n'], n_groups,total_pairs);
    end
    group_cells = cell(n_groups,1);
    parfor gg = 1:n_groups
        group_cells{gg} = buildLaplacePairGroup(gg,rep_rows(gg),group_sep(gg), ...
            q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
            project_charge);
    end
    pair_cache.groups = vertcat(group_cells{:});
    pair_cache.stats.pool_size = get_parallel_pool_size();
else
    group_cells = cell(n_groups,1);
    covered_pairs = 0;
    for gg = 1:n_groups
        group_cells{gg} = buildLaplacePairGroup(gg,rep_rows(gg),group_sep(gg), ...
            q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
            project_charge);

        if show_counter
            rep_i = pairs(rep_rows(gg),1);
            rep_j = pairs(rep_rows(gg),2);
            n_rep = sum(group_id == gg);
            covered_pairs = covered_pairs + n_rep;
            fprintf(['getPairBasisLaplace: processed canonical group %d/%d ', ...
                'from pair (%d,%d), pp_sep = %.3g, covers %d pairs -> %d/%d pairs covered\n'], ...
                gg,n_groups,rep_i,rep_j,group_sep(gg)-2*R,n_rep,covered_pairs,total_pairs);
        end
    end
    pair_cache.groups = vertcat(group_cells{:});
end

if use_pair_map && N_peanut && ~isempty(pair_cache.groups)
    report_cmap_size(pair_cache.groups(1).Cmap,N_cmap);
end

if check_rotations
    % Keep a per-pair copy in the solve geometry so the debug path can
    % compare against the canonical cached group later.
    check_pair_cells = cell(size(pairs,1),1);
    for row = 1:size(pairs,1)
        gid = group_id(row);
        check_pair_cells{row} = buildLaplacePairGroup(gid,row, ...
            pair_cache.meta(row).sep,q,rbase_in_c,rbase_in_f, ...
            rimage_vec,refine,pairs,opt,pair_cache.rout_base_f, ...
            project_charge,false);
    end
    pair_cache.check_pairs = vertcat(check_pair_cells{:});
end

for row = 1:size(pairs,1)
    gid = group_id(row);
    pair_cache.meta(row).group_id = gid;
    pair_cache.meta(row).sep = group_sep(gid);
    if use_pair_map
        [Ucross_actual,Ec_actual,Lr_actual] = buildLaplaceActualPairCollocFactors( ...
            pair_cache.meta(row),q,rbase_in_solver,rout_base_c);
        pair_cache.meta(row).Ucross_colloc_actual = Ucross_actual;
        pair_cache.meta(row).Ec_colloc_actual = Ec_actual;
        pair_cache.meta(row).Lr_colloc_actual = Lr_actual;
    end
end

if show_counter
    fprintf('getPairBasisLaplace: reused %d canonical groups for %d pairs\n', ...
        n_groups,size(pairs,1));
end

end

function [use_parallel,pool_size] = resolve_parallel_precomp(requested_parallel,caller_name)
use_parallel = requested_parallel;
pool_size = 0;
if ~requested_parallel
    return
end

if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox') || ...
        exist('gcp','file') ~= 2
    error([caller_name ':ParallelToolboxRequired'], ...
        ['opt.parallel_precomp requires Parallel Computing Toolbox. ', ...
         'Open a pool before benchmarking warm-pool speedups if you want ', ...
         'to exclude startup overhead.']);
end

pool_size = get_parallel_pool_size();
end

function pool_size = get_parallel_pool_size()
pool_size = 0;
pool = gcp('nocreate');
if ~isempty(pool)
    pool_size = pool.NumWorkers;
end
end

function report_cmap_size(C,N_cmap)
if isempty(C)
    return
end
fprintf(['getPairBasisLaplace: Cmap size = %d x %d ', ...
    '(N_cmap=%d, square aligned pair map)\n'], ...
    size(C,1),size(C,2),N_cmap);
end

function [Uf,Yf,Up,Yp,Cmap,Cmap_QV] = assign_pair_entry_outputs( ...
    Uf,Yf,Up,Yp,Cmap,Cmap_QV,pairs,row,pair_entry,N_peanut,use_pair_map)
i = pairs(row,1);
j = pairs(row,2);
Uf{i,j} = pair_entry.Upf;
Yf{i,j} = pair_entry.Ypf;

if N_peanut
    if use_pair_map
        Cmap{i,j} = pair_entry.Cmap;
        Cmap_QV{i,j} = pair_entry.Cmap_QV;
    else
        Up{i,j} = pair_entry.DC;
        Yp{i,j} = pair_entry.YC;
    end
end
end
