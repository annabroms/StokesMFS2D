function [pair, debug_data] = buildStokesMobilityPairData(...
    q, rbase_in_c, rbase_in_f, rimage_pairs, refine, pairs, opt, Lc, row,...
    debug, svd_opts, use_canonical, payload_mode)
%BUILDSTOKESMOBILITYPAIRDATA Build one Stokes mobility pair payload.
%
% This helper contains the pair solve used by getPairBasisStokes and by the
% streamed Stokes big-sparse builder. The caller decides whether to keep the
% full fine-source payload or only the compressed maps.

    % --- Default arguments ---
    if nargin < 10 || isempty(debug);         debug         = false;      end
    if nargin < 11 || isempty(svd_opts);      svd_opts      = struct();   end
    if nargin < 12 || isempty(use_canonical); use_canonical = false;      end
    if nargin < 13 || isempty(payload_mode);  payload_mode  = 'full';     end

    store_full = strcmp(payload_mode, 'full');

    % --- Cache frequently used opt fields up front ---
    N_c             = opt.N_c;
    a_c             = getOptField(opt, 'a_c', 1.2);
    a_f             = opt.a_f;
    N_f             = opt.N_f;
    project_force   = logical(getOptField(opt, 'project_force', false));
    use_cmap        = logical(getOptField(opt, 'cmap',          false));
    N_peanut        = getOptField(opt, 'N_peanut', 0);
    use_left_weight = logical(getOptField(svd_opts, 'left_weight', false));

    % --- Ensure column vectors ---
    q          = q(:);
    rbase_in_c = rbase_in_c(:);
    rbase_in_f = rbase_in_f(:);

    % --- Pair indices ---
    i = pairs(row, 1);
    j = pairs(row, 2);
    pair = init_pair_payload(i, j);

    % --- Geometry ---
    delta = q(j) - q(i);
    sep   = abs(delta);

    if sep == 0
        rot = 1;
    else
        rot = delta / sep;
    end
    mid = 0.5 * (q(i) + q(j));

    % --- Canonical or physical frame ---
    if use_canonical
        q_pair   = [-sep/2; sep/2];
        rimage_i = map_points_to_canonical(rimage_pairs{i,j}, mid, rot);
        rimage_j = map_points_to_canonical(rimage_pairs{j,i}, mid, rot);
        refine_i = map_points_to_canonical(refine{i,j},       mid, rot);
        refine_j = map_points_to_canonical(refine{j,i},       mid, rot);
    else
        q_pair   = [q(i); q(j)];
        rimage_i = rimage_pairs{i,j};
        rimage_j = rimage_pairs{j,i};
        refine_i = refine{i,j};
        refine_j = refine{j,i};
    end

    pair.q_pair       = q_pair;
    pair.rimage_canon = {rimage_i; rimage_j};
    pair.refine_canon = {refine_i; refine_j};

    % --- Quadrature nodes via exp(1i*t) — single trig call each ---
    nout_c      = ceil(a_c * N_c);
    rout_base_c = exp(1i * linspace(0, 2*pi*(1-1/nout_c), nout_c)');

    nout_f    = ceil(a_f * N_f);
    rout_base = exp(1i * linspace(0, 2*pi*(1-1/nout_f), nout_f)');

    % --- Source and target point arrays ---
    rin_1_f  = q_pair(1) + rbase_in_f;
    rin_2_f  = q_pair(2) + rbase_in_f;
    rout_1   = q_pair(1) + rout_base;
    rout_2   = q_pair(2) + rout_base;
    rout_f   = [rout_1; refine_i; rout_2; refine_j];
    rin_pair = [rin_1_f; rimage_i; rin_2_f; rimage_j];

    % Precompute half-index for rin_pair (used multiple times)
    rin_half = numel(rin_pair) / 2;

    if store_full
        pair.rin_pair = rin_pair;
    end

    % --- SVD weights ---
    svd_pair = svd_opts;
    if use_left_weight
        svd_pair.row_weights = [...
            getPeriodicCurveWeights([rout_1; refine_i], q_pair(1));...
            getPeriodicCurveWeights([rout_2; refine_j], q_pair(2))];
    end

    % --- Debug geometry plot ---
    if debug
        plot_pair_geometry(i, j, q_pair, rin_1_f, rin_2_f, rout_base,...
            refine_i, refine_j, rimage_i, rimage_j);
    end

    % --- Coarse correction operator ---
    Lc_pair = [];
    if ~isempty(Lc)
        Lc_pair = getILpair(Lc);
    end

    % --- Moment / RBM maps (shared by project_force and cmap) ---
    if project_force || use_cmap
        Kf1              = getKmat2D(rin_pair(1:rin_half),     q_pair(1));
        Kf2              = getKmat2D(rin_pair(rin_half+1:end), q_pair(2));
        pair_moment_map  = getKftPair(Kf1, Kf2);
        pair_rbm_map     = pair_moment_map';
        pair_moment_gram = pair_moment_map * pair_rbm_map;
    else
        Kf1 = []; Kf2 = [];
        pair_moment_map  = [];
        pair_rbm_map     = [];
        pair_moment_gram = [];
    end

    % --- Projection operators ---
    if project_force
        B1 = getKmat2D([rout_1; refine_i], q_pair(1));
        B2 = getKmat2D([rout_2; refine_j], q_pair(2));
        Lf_pair               = getLfPair(Kf1, Kf2);
        pair_target_rbm_map   = getKftPair(B1, B2)';
        pair_proj_moment_map  = pair_moment_map;
        pair_proj_rbm_map     = pair_rbm_map;
        pair_proj_moment_gram = pair_moment_gram;
    else
        Lf_pair               = [];
        pair_target_rbm_map   = [];
        pair_proj_moment_map  = [];
        pair_proj_rbm_map     = [];
        pair_proj_moment_gram = [];
    end

    % --- Core pair solve ---
    [Uf_pair, Yf_pair] = getPairBlockStokes(rin_pair, rout_f,...
        pair_proj_moment_map, pair_proj_rbm_map,...
        pair_proj_moment_gram, pair_target_rbm_map, svd_pair);

    % Separate the unsigned intermediate UfTN = Uf_pair'*Npair
    % from the stored Upf = -UfTN, to avoid sign confusion in Cmap
    % expressions. UfTN is a direct substitute for Uf_pair'*Npair everywhere.
    Npair = evaluateCoarseOnPair(q_pair, rbase_in_c, rout_f);
    UfTN  = Uf_pair' * Npair;    % = Uf_pair'*Npair, used in Cmap expressions
    Upf   = -UfTN;               % stored quantity, used in colloc block

    % --- Coarse (peanut) block ---
    DC = []; YC = [];
    rin_pair_c = [q_pair(1) + rbase_in_c; q_pair(2) + rbase_in_c];

    if N_peanut > 0
        rout_peanut = createPeanut(q_pair(1), q_pair(2), N_peanut, false);
        [DC, YC] = getPeanutBlockStokes(rin_pair_c, rin_pair, rout_peanut,...
            Lc_pair, pair_proj_moment_map, pair_proj_rbm_map,...
            pair_proj_moment_gram, svd_opts);
        if use_cmap           
            pair.Cmap = -YC * (DC * Yf_pair * UfTN);
        end
    end

    if use_cmap
        pair.Cmap_FU = -pair_moment_map * Yf_pair * UfTN;
    end

    % --- Store full payload ---
    if store_full
        pair.Upf     = Upf;
        pair.Ypf     = Yf_pair;
        pair.Lf_pair = Lf_pair;
        pair.Lc_pair = Lc_pair;
        pair.DC      = DC;
        pair.YC      = YC;
    end

    % --- Collocation maps (only when full payload + project_force) ---
    if store_full && project_force
        rout_pair_c = [q_pair(1) + rout_base_c; q_pair(2) + rout_base_c];
        Ppair       = eye(size(Lf_pair)) - Lf_pair;

        % Upf = -UfTN is correct here — matches original definition
        pair.Upair_colloc  = stokSLPmat(rin_pair, rout_pair_c, 1) *...
                             Ppair * Yf_pair * Upf;
        pair.Ecolloc       = stokSLPmat(rin_pair_c, rout_pair_c, 1);
        pair.Ucross_colloc = buildStokesCrossPairVelocityMap(...
                             pair.Ecolloc, N_c, numel(rout_base_c));
    end

    if debug
        debug_data = struct(...
            'q_pair',              q_pair,...
            'rin_pair',            rin_pair,...
            'Upf',                 Upf,...
            'Yf_pair',             Yf_pair,...
            'Lf_pair',             Lf_pair,...
            'Kf1',                 Kf1,...
            'Kf2',                 Kf2,...
            'pair_moment_map',     pair_moment_map,...
            'pair_rbm_map',        pair_rbm_map,...
            'pair_moment_gram',    pair_moment_gram,...
            'pair_target_rbm_map', pair_target_rbm_map,...
            'Lc_pair',             Lc_pair,...
            'DC',                  DC,...
            'YC',                  YC,...
            'rout_f',              rout_f,...
            'Npair',               Npair...
        );
    else
        debug_data = [];
    end

end

% =========================================================================

function pair = init_pair_payload(i, j)
    pair = struct(...
        'rep_pair',      [i j],...
        'q_pair',        [],...
        'rimage_canon',  {{[]; []}},...
        'refine_canon',  {{[]; []}},...
        'rin_pair',      [],...
        'Upf',           [],...
        'Ypf',           [],...
        'Lf_pair',       [],...
        'Lc_pair',       [],...
        'DC',            [],...
        'YC',            [],...
        'Cmap',          [],...
        'Cmap_FU',       [],...
        'Upair_colloc',  [],...
        'Ucross_colloc', [],...
        'Ecolloc',       []...
    );
end

% =========================================================================

function plot_pair_geometry(i, j, q_pair, rin_1_f, rin_2_f, rout_base,...
        refine_i, refine_j, rimage_i, rimage_j)
    figure(801); clf;
    plot(real(rin_1_f), imag(rin_1_f), 'r.', 'MarkerSize', 10); hold on;
    plot(real(rin_2_f), imag(rin_2_f), 'b.', 'MarkerSize', 10);
    plot(real(q_pair(1)+rout_base), imag(q_pair(1)+rout_base), 'ro', 'MarkerSize', 4);
    plot(real(q_pair(2)+rout_base), imag(q_pair(2)+rout_base), 'bo', 'MarkerSize', 4);
    plot(real(refine_i), imag(refine_i), 'r+', 'MarkerSize', 6);
    plot(real(refine_j), imag(refine_j), 'b+', 'MarkerSize', 6);
    if ~isempty(rimage_i)
        plot(real(rimage_i), imag(rimage_i), 'ks', 'MarkerSize', 5);
    end
    if ~isempty(rimage_j)
        plot(real(rimage_j), imag(rimage_j), 'kd', 'MarkerSize', 5);
    end
    axis equal; grid on;
    title(sprintf('getPairBasisStokes pair (%d,%d)', i, j), 'Interpreter', 'none');
    drawnow;
end

% =========================================================================

function z_canon = map_points_to_canonical(z, mid, rot)
    z = z(:);
    if isempty(z)
        z_canon = zeros(0, 1);
        return
    end
    z_canon = conj(rot) * (z - mid);
end