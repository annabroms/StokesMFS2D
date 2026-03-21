function [rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx, tau_stokes_nonpy,tau_stokes_e_nonpx, tau_stokes_e_nonpy, rimage_k] = getMobPairTransformationStokes(tau,geom,basis)
%GETMOBPAIRTRANSFORMATIONSTOKES Map coarse boundary data to coarse and fine Stokes source strengths.
% This implementation caches coarse lambda for all particles first, then loops over pairs.
%
% Syntax:
%   [rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y,...
%    tau_stokes_nonpx,tau_stokes_nonpy,tau_stokes_e_nonpx,tau_stokes_e_nonpy,...
%    rimage_k] = getMobPairTransformationStokes(tau,geom,basis)
%
% Inputs:
%   tau   - Coarse data at collocation points [tau_x; tau_y].
%   geom  - Struct with fields:
%           rbase_in_c, rbase_in_f, rimage_vec, opt, rvec_out, q, pairs.
%   basis - Struct with fields:
%           U, Y, Lc, Upf, Ypf. (precomputed SVD factors)

rvec_out = geom.rvec_out;
q = geom.q;
pairs = geom.pairs;

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
rimage_vec = geom.rimage_vec;
opt = geom.opt;

U = basis.U;
Y = basis.Y;
Lc = basis.Lc;
Upf = basis.Upf;
Ypf = basis.Ypf;
if isfield(basis,'pair_cache')
    pair_cache = basis.pair_cache;
else
    pair_cache = struct('enabled',false);
end
use_pair_cache = isfield(pair_cache,'enabled') && pair_cache.enabled;
project_coarse_in_reference_frame = logical( ...
    getOptField(opt,'project_coarse_in_reference_frame',false));
project_fine_in_reference_frame = logical( ...
    getOptField(opt,'project_fine_in_reference_frame',false));
use_matrix_free_Lc_pair = logical(getOptField(opt,'use_matrix_free_Lc_pair',true));

P = length(q);
N_coarse = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;
PM = length(rvec_out);

if project_coarse_in_reference_frame
    if isempty(Lc)
        error('getMobPairTransformationStokes:MissingLc', ...
            ['project_coarse_in_reference_frame requires the single-body ', ...
             'coarse projection matrix Lc.']);
    end
    Lc_pair = getILpair(Lc);
else
    Lc_pair = [];
end

% Preallocate coarse contributions (same size for all particles).
n_coarse_tot = P*N_coarse;
tau_stokes_x_coarse = zeros(n_coarse_tot,1);
tau_stokes_y_coarse = zeros(n_coarse_tot,1);
tau_stokes_nonpx_coarse = zeros(n_coarse_tot,1);
tau_stokes_nonpy_coarse = zeros(n_coarse_tot,1);
rvec_in_coarse = zeros(n_coarse_tot,1);
coarse_ind = cell(P,1);

% Prepaare to store fine sources per body
tau_stokes_f_x = cell(P,1);
tau_stokes_f_y = cell(P,1);
tau_stokes_e_x = cell(P,1);
tau_stokes_e_y = cell(P,1);
tau_stokes_f_nonpx = cell(P,1);
tau_stokes_f_nonpy = cell(P,1);

tau_stokes_e_nonpx = cell(P,1);
tau_stokes_e_nonpy = cell(P,1);

rimage_k = cell(P,1);
rimage_k_chunks = repmat({cell(0,1)},P,1);
tau_stokes_e_x_chunks = repmat({cell(0,1)},P,1);
tau_stokes_e_y_chunks = repmat({cell(0,1)},P,1);
tau_stokes_e_nonpx_chunks = repmat({cell(0,1)},P,1);
tau_stokes_e_nonpy_chunks = repmat({cell(0,1)},P,1);

%% Phase 1: map all bodies once on the coarse grid to recover proxy source
% strengths from data at the boundary.
for i = 1:P
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    step1 = U{1}*[tau_particle_x; tau_particle_y];
    lambda_coarse_nonproj_i = Y{1}*step1;

    tau_i_x = lambda_coarse_nonproj_i(1:N_coarse);
    tau_i_y = lambda_coarse_nonproj_i(N_coarse+1:end);
    lambda_coarse_i = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y];

    coarse_ind_i = (i-1)*N_coarse+1:i*N_coarse;
    coarse_ind{i} = coarse_ind_i;

    tau_stokes_nonpx_coarse(coarse_ind_i) = lambda_coarse_nonproj_i(1:N_coarse);
    tau_stokes_nonpy_coarse(coarse_ind_i) = lambda_coarse_nonproj_i(N_coarse+1:end);
    tau_stokes_x_coarse(coarse_ind_i) = lambda_coarse_i(1:N_coarse);
    tau_stokes_y_coarse(coarse_ind_i) = lambda_coarse_i(N_coarse+1:end);
    rvec_in_coarse(coarse_ind_i) = q(i)+rbase_in_c;
end

%% Phase 2: loop over close pairs and build fine-grid corrections.
for pair_row = 1:size(pairs,1)
    i = pairs(pair_row,1);
    p2 = pairs(pair_row,2);

    if project_coarse_in_reference_frame
        % Experimental path: rotate the non-projected coarse densities to the
        % pair frame first, then apply the coarse rigid-mode projection there.
        lambda_coarse_i = [tau_stokes_nonpx_coarse(coarse_ind{i}); ...
            tau_stokes_nonpy_coarse(coarse_ind{i})];
        lambda_coarse_p2 = [tau_stokes_nonpx_coarse(coarse_ind{p2}); ...
            tau_stokes_nonpy_coarse(coarse_ind{p2})];
    else
        % Default/original path: project each body's coarse sources before
        % assembling the pair-local coarse data.
        lambda_coarse_i = [tau_stokes_x_coarse(coarse_ind{i}); ...
            tau_stokes_y_coarse(coarse_ind{i})];
        lambda_coarse_p2 = [tau_stokes_x_coarse(coarse_ind{p2}); ...
            tau_stokes_y_coarse(coarse_ind{p2})];
    end

    rimage_i = rimage_vec{i,p2};
    rimage_p2 = rimage_vec{p2,i};

    rimage_k_chunks{i}{end+1,1} = rimage_i;
    rimage_k_chunks{p2}{end+1,1} = rimage_p2;

    % Keep track of local ordering of source vector for the pair.
    im_nr = length(rimage_i);
    s_ind1_x = 1:opt.N_f;
    s_ind2_x = opt.N_f+im_nr+1:2*opt.N_f+im_nr;
    s_ind1_y = 2*opt.N_f+2*im_nr+1:3*opt.N_f+2*im_nr;
    s_ind2_y = 3*opt.N_f+3*im_nr+1:4*opt.N_f+3*im_nr;
    e_ind1_x = opt.N_f+1:opt.N_f+im_nr;
    e_ind2_x = 2*opt.N_f+im_nr+1:2*opt.N_f+2*im_nr;
    e_ind1_y = 3*opt.N_f+2*im_nr+1:3*opt.N_f+3*im_nr;
    e_ind2_y = 4*opt.N_f+3*im_nr+1:4*opt.N_f+4*im_nr;

    % Fine source locations for the pair.
    rin_pair = [rbase_in_f+q(i); rimage_i; rbase_in_f+q(p2); rimage_p2];

    if use_pair_cache
        pair = getStokesPairInstance(pair_cache,pair_row);

        %rotate source to frame of reference
        coarse_to_fine_tot = rotatePairOrderedStokesData( ...
            [lambda_coarse_i(1:N_coarse) lambda_coarse_p2(1:N_coarse) ...
             lambda_coarse_i(N_coarse+1:end) lambda_coarse_p2(N_coarse+1:end)], ...
            N_coarse,pair.meta.phase_c,conj(pair.meta.rot));
        coarse_to_fine_tot = coarse_to_fine_tot(:);
    else
        %same thing, but without rotating to reference frame and back
        coarse_to_fine_tot = [lambda_coarse_i(1:end/2); lambda_coarse_p2(1:end/2); ...
            lambda_coarse_i(end/2+1:end); lambda_coarse_p2(end/2+1:end)];
    end

    if project_coarse_in_reference_frame
        if use_matrix_free_Lc_pair
            if use_pair_cache
                coarse_to_fine_tot = projectOutRigidPair2D( ...
                    coarse_to_fine_tot, ...
                    rbase_in_c + pair.group.q_pair(1), pair.group.q_pair(1), ...
                    rbase_in_c + pair.group.q_pair(2), pair.group.q_pair(2));
            else
                coarse_to_fine_tot = projectOutRigidPair2D( ...
                    coarse_to_fine_tot, ...
                    rbase_in_c+q(i), q(i), rbase_in_c+q(p2), q(p2));
            end
        else
            coarse_to_fine_tot = Lc_pair*coarse_to_fine_tot;
        end
    end

    if use_pair_cache
        % Take pseudoinverse of the fine representation to determine fine
        % sources for both pair corrections related to the pair.
        pair_mapped = pair.group.Upf*coarse_to_fine_tot;
        beta_tot_ref = pair.group.Ypf*pair_mapped;

        if project_fine_in_reference_frame
            rin_pair_ref_i = [rbase_in_f+pair.group.q_pair(1); ...
                pair.group.rimage_canon{1}];
            rin_pair_ref_p2 = [rbase_in_f+pair.group.q_pair(2); ...
                pair.group.rimage_canon{2}];

            tau_fine_ref_i = extract_pair_body_sources(beta_tot_ref, ...
                s_ind1_x,e_ind1_x,s_ind1_y,e_ind1_y);
            tau_fine_ref_p2 = extract_pair_body_sources(beta_tot_ref, ...
                s_ind2_x,e_ind2_x,s_ind2_y,e_ind2_y);

            tau_mapped_proj_ref_i = projectOutRigid2D( ...
                tau_fine_ref_i, rin_pair_ref_i, pair.group.q_pair(1));
            tau_mapped_proj_ref_p2 = projectOutRigid2D( ...
                tau_fine_ref_p2, rin_pair_ref_p2, pair.group.q_pair(2));

            beta_tot_proj_ref = beta_tot_ref;
            beta_tot_proj_ref = insert_pair_body_sources(beta_tot_proj_ref, ...
                tau_mapped_proj_ref_i,s_ind1_x,e_ind1_x,s_ind1_y,e_ind1_y);
            beta_tot_proj_ref = insert_pair_body_sources(beta_tot_proj_ref, ...
                tau_mapped_proj_ref_p2,s_ind2_x,e_ind2_x,s_ind2_y,e_ind2_y);

            beta_tot_proj = rotateStokesPairSourceVector( ...
                beta_tot_proj_ref,N_f,length(rimage_i),length(rimage_p2), ...
                pair.meta.phase_f_inv,pair.meta.rot);
        else
            beta_tot_proj = [];
        end

        beta_tot = rotateStokesPairSourceVector( ...
            beta_tot_ref,N_f,length(rimage_i),length(rimage_p2), ...
            pair.meta.phase_f_inv,pair.meta.rot);
    else
        pair_mapped = Upf{i,p2}*coarse_to_fine_tot;
        beta_tot = Ypf{i,p2}*pair_mapped;
        beta_tot_proj = [];
    end
        
    %% Project to remove force/torque-producing modes in fine source beta_tot

    % Reorder unknowns
    tau_mapped_f_xi = beta_tot(s_ind1_x);
    tau_mapped_f_yi = beta_tot(s_ind1_y);
    tau_mapped_f_xp2 = beta_tot(s_ind2_x);
    tau_mapped_f_yp2 = beta_tot(s_ind2_y);

    tau_mapped_e_xi = beta_tot(e_ind1_x);
    tau_mapped_e_yi = beta_tot(e_ind1_y);
    tau_mapped_e_xp2 = beta_tot(e_ind2_x);
    tau_mapped_e_yp2 = beta_tot(e_ind2_y);
    tau_fine_i = [tau_mapped_f_xi; tau_mapped_e_xi; tau_mapped_f_yi; tau_mapped_e_yi];
    tau_fine_p2 = [tau_mapped_f_xp2; tau_mapped_e_xp2; tau_mapped_f_yp2; tau_mapped_e_yp2];

    if use_pair_cache && project_fine_in_reference_frame
        tau_mapped_proj_i = extract_pair_body_sources(beta_tot_proj, ...
            s_ind1_x,e_ind1_x,s_ind1_y,e_ind1_y);
        tau_mapped_proj_p2 = extract_pair_body_sources(beta_tot_proj, ...
            s_ind2_x,e_ind2_x,s_ind2_y,e_ind2_y);
    else
        % Default/original path: rotate back first, then project on the
        % physical pair geometry.
        tau_mapped_proj_i = projectOutRigid2D( ...
            tau_fine_i, rin_pair(1:end/2), q(i));
        tau_mapped_proj_p2 = projectOutRigid2D( ...
            tau_fine_p2, rin_pair(end/2+1:end), q(p2));
    end


    %% Store fine source strengths for later evaluation.
    pair_idx = [i; p2];
    pair_proj = {tau_mapped_proj_i, tau_mapped_proj_p2};
    pair_f_nonpx = {tau_mapped_f_xi, tau_mapped_f_xp2};
    pair_f_nonpy = {tau_mapped_f_yi, tau_mapped_f_yp2};
    pair_e_nonpx = {tau_mapped_e_xi, tau_mapped_e_xp2};
    pair_e_nonpy = {tau_mapped_e_yi, tau_mapped_e_yp2};

    for pair_it = 1:2
        idx = pair_idx(pair_it);
        tau_proj = pair_proj{pair_it};

        tau_f_x = tau_proj(s_ind1_x);
        tau_f_y = tau_proj(opt.N_f+im_nr+1:2*opt.N_f+im_nr);
        tau_e_x = tau_proj(e_ind1_x);
        tau_e_y = tau_proj(2*opt.N_f+im_nr+1:2*opt.N_f+2*im_nr);

        if isempty(tau_stokes_f_x{idx})
            tau_stokes_f_x{idx} = tau_f_x;
            tau_stokes_f_y{idx} = tau_f_y;
            tau_stokes_f_nonpx{idx} = pair_f_nonpx{pair_it};
            tau_stokes_f_nonpy{idx} = pair_f_nonpy{pair_it};
        else
            tau_stokes_f_x{idx} = tau_stokes_f_x{idx} + tau_f_x;
            tau_stokes_f_y{idx} = tau_stokes_f_y{idx} + tau_f_y;
            tau_stokes_f_nonpx{idx} = tau_stokes_f_nonpx{idx} + pair_f_nonpx{pair_it};
            tau_stokes_f_nonpy{idx} = tau_stokes_f_nonpy{idx} + pair_f_nonpy{pair_it};
        end

        tau_stokes_e_x_chunks{idx}{end+1,1} = tau_e_x;
        tau_stokes_e_y_chunks{idx}{end+1,1} = tau_e_y;
        tau_stokes_e_nonpx_chunks{idx}{end+1,1} = pair_e_nonpx{pair_it};
        tau_stokes_e_nonpy_chunks{idx}{end+1,1} = pair_e_nonpy{pair_it};
    end
end

%% Recover arrays of all Stokeslet source locations and their strengths.
for k = 1:P
    if isempty(rimage_k_chunks{k})
        rimage_k{k} = zeros(0,1);
    else
        rimage_k{k} = vertcat(rimage_k_chunks{k}{:});
    end

    if isempty(tau_stokes_e_x_chunks{k})
        tau_stokes_e_x{k} = zeros(0,1);
        tau_stokes_e_y{k} = zeros(0,1);
        tau_stokes_e_nonpx{k} = zeros(0,1);
        tau_stokes_e_nonpy{k} = zeros(0,1);
    else
        tau_stokes_e_x{k} = vertcat(tau_stokes_e_x_chunks{k}{:});
        tau_stokes_e_y{k} = vertcat(tau_stokes_e_y_chunks{k}{:});
        tau_stokes_e_nonpx{k} = vertcat(tau_stokes_e_nonpx_chunks{k}{:});
        tau_stokes_e_nonpy{k} = vertcat(tau_stokes_e_nonpy_chunks{k}{:});
    end
end

has_neigh = sort(unique(pairs(:)));
n_has = numel(has_neigh);

% Collect all source points and source data using one-time concatenation.
x_chunks = cell(1+2*n_has,1);
y_chunks = cell(1+2*n_has,1);
nonpx_chunks = cell(1+n_has,1);
nonpy_chunks = cell(1+n_has,1);
rvec_chunks = cell(1+2*n_has,1);

x_chunks{1} = tau_stokes_x_coarse;
y_chunks{1} = tau_stokes_y_coarse;
nonpx_chunks{1} = tau_stokes_nonpx_coarse;
nonpy_chunks{1} = tau_stokes_nonpy_coarse;
rvec_chunks{1} = rvec_in_coarse;

for i = 1:n_has
    k = has_neigh(i);
    x_chunks{2*i} = tau_stokes_f_x{k};
    x_chunks{2*i+1} = tau_stokes_e_x{k};
    y_chunks{2*i} = tau_stokes_f_y{k};
    y_chunks{2*i+1} = tau_stokes_e_y{k};
    nonpx_chunks{i+1} = tau_stokes_f_nonpx{k};
    nonpy_chunks{i+1} = tau_stokes_f_nonpy{k};
    rvec_chunks{2*i} = q(k)+rbase_in_f;
    rvec_chunks{2*i+1} = rimage_k{k};
end

tau_stokes_x = vertcat(x_chunks{:});
tau_stokes_y = vertcat(y_chunks{:});
tau_stokes_nonpx = vertcat(nonpx_chunks{:});
tau_stokes_nonpy = vertcat(nonpy_chunks{:});
rvec_in = vertcat(rvec_chunks{:});

end

function tau_body = extract_pair_body_sources(beta_tot,idx_fx,idx_ex,idx_fy,idx_ey)
tau_body = [beta_tot(idx_fx); beta_tot(idx_ex); beta_tot(idx_fy); beta_tot(idx_ey)];
end

function beta_tot = insert_pair_body_sources(beta_tot,tau_body,idx_fx,idx_ex,idx_fy,idx_ey)
n_f = numel(idx_fx);
n_e = numel(idx_ex);

beta_tot(idx_fx) = tau_body(1:n_f);
beta_tot(idx_ex) = tau_body(n_f+1:n_f+n_e);
beta_tot(idx_fy) = tau_body(n_f+n_e+1:2*n_f+n_e);
beta_tot(idx_ey) = tau_body(2*n_f+n_e+1:2*(n_f+n_e));
end
