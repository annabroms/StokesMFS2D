function snapshot = buildCanonicalLaplacePairMap(delta,opt,grids,problem)
%BUILDCANONICALLAPLACEPAIRMAP Exact aligned Cmap and auxiliary map at one gap.
%
% This is the shared production snapshot constructor used by both the
% runtime interpolator and the Sep. 25 experiments.

if nargin < 4
    problem = getOptField(opt,'interpolation_problem','capacitance');
end
problem = resolveLaplaceInterpolationProblem(problem,mfilename);

R = opt.rad;
sep = 2*R+delta;
q = [0;sep];
opt_exact = opt;
opt_exact.use_interpolation = 'none';
opt_exact.reuse_pair_basis_by_sep = false;
opt_exact.check_rotations = false;
opt_exact.show_counter = 0;
opt_exact.delta_pair = max(opt.delta_pair,delta)*(1+1e-12);

[~,~,~,rimage_vec,refine,pairs] = getEnhancedGrid(q,opt_exact);
if ~isequal(pairs,[1 2])
    error('buildCanonicalLaplacePairMap:BadPair', ...
        'Expected the canonical pair (1,2) at delta=%.16g.',delta);
end

group = buildLaplacePairGroup(1,1,sep,q,grids.rbase_in_c, ...
    grids.rbase_in_f,rimage_vec,refine,pairs,opt_exact, ...
    grids.rout_base_f,strcmp(problem,'elastance'),true);
snapshot = struct('delta',delta,'Cmap',group.Cmap, ...
    'Cmap_QV',group.Cmap_QV,'problem',problem);
end
