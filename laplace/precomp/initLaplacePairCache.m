function pair_cache = initLaplacePairCache()
%INITLAPLACEPAIRCACHE Empty cache struct returned by getPairBasisLaplace.

empty_group = struct('group_id',{},'sep',{},'q_pair',{}, ...
    'rimage_canon',{},'refine_canon',{},'Upf',{},'Ypf',{}, ...
    'DC',{},'YC',{},'Cmap',{},'Cmap_QV',{},'nout_f',{}, ...
    'nsrc_f',{},'ntar_f',{},'rep_pair',{});

pair_cache = struct();
pair_cache.enabled = false;
pair_cache.check_rotations = false;
pair_cache.shared_sep_tol = [];
pair_cache.rout_base_f = [];
pair_cache.meta = repmat(struct('i',[],'j',[],'group_id',[],'sep',[], ...
    'mid',[],'rot',[],'phase_c',[],'phase_c_inv',[], ...
    'phase_f',[],'phase_f_inv',[], ...
    'Ucross_colloc_actual',[],'Ec_colloc_actual',[], ...
    'Lr_colloc_actual',[]),0,1);
pair_cache.groups = empty_group;
pair_cache.check_pairs = empty_group;
pair_cache.group_id = zeros(0,1);
pair_cache.group_sep = zeros(0,1);
pair_cache.representative_rows = zeros(0,1);
pair_cache.n_groups = 0;
pair_cache.stats = struct('requested_parallel',false,'used_parallel',false, ...
    'branch','','n_pairs',0,'n_groups',0,'pool_size',0);
end
