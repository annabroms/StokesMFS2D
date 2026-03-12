function [rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y,u_corr] = getPairTransformationStokes(tau,geom,basis)
%GETPAIRTRANSFORMATIONSTOKES Map coarse boundary data to interior coarse and fine Stokes source strengths.
%
% Syntax:
%   [rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y,u_corr] = ...
%       getPairTransformationStokes(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] at outer collocation points.
%   geom  - Struct with fields:
%           rbase_in_c, rbase_in_f, refine, rimage_vec, opt, rvec_out, q, pairs.
%   basis - Struct with fields:
%           U, Y, Upf, Ypf. (precomputed SVD factors)
%
% Outputs:
%   rvec_in      - All source locations used in the matvec.
%   coarse_ind   - Coarse-source index set per particle.
%   tau_stokes_x - x-strengths at rvec_in.
%   tau_stokes_y - y-strengths at rvec_in.
%   u_corr       - Pair/self correction used to create identity blocks in the system matrix.

rvec_out = geom.rvec_out;
q = geom.q;
pairs = geom.pairs;

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
refine = geom.refine;
rimage_vec = geom.rimage_vec;
opt = geom.opt;

U = basis.U;
Y = basis.Y;
Upf = basis.Upf;
Ypf = basis.Ypf;

P = length(q);
N_coarse = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;
PM = length(rvec_out);
precomp = opt.precomp;

% Preallocate coarse contributions (same size for all particles).
n_coarse_tot = P*N_coarse;
tau_stokes_x_coarse = zeros(n_coarse_tot,1);
tau_stokes_y_coarse = zeros(n_coarse_tot,1);
rvec_in_coarse = zeros(n_coarse_tot,1);

coarse_ind = cell(P,1);

% Fine/proxy and extra/image strengths, accumulated per particle.
tau_stokes_f_x = cell(P,1);
tau_stokes_f_y = cell(P,1);
tau_stokes_e_x = cell(P,1);
tau_stokes_e_y = cell(P,1);
rimage_k = cell(P,1);

rimage_k_chunks = repmat({cell(0,1)},P,1);
tau_stokes_e_x_chunks = repmat({cell(0,1)},P,1);
tau_stokes_e_y_chunks = repmat({cell(0,1)},P,1);

u_corr = zeros(2*PM,1);

%% Phase 1: map all bodies once to recover coarse source strenghts from 
% coarse velocity data
for i = 1:P
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    step1 = U{1}*[tau_particle_x; tau_particle_y];
    lambda_coarse_i = Y{1}*step1;

    coarse_ind_i = (i-1)*N_coarse+1:i*N_coarse;
    coarse_ind{i} = coarse_ind_i;
    tau_stokes_x_coarse(coarse_ind_i) = lambda_coarse_i(1:N_coarse);
    tau_stokes_y_coarse(coarse_ind_i) = lambda_coarse_i(N_coarse+1:end);
    rvec_in_coarse(coarse_ind_i) = q(i)+rbase_in_c;
end

%% Phase 2: loop once over each close pair
for pair_row = 1:size(pairs,1)
    i = pairs(pair_row,1);
    p2 = pairs(pair_row,2);

    lambda_i = [tau_stokes_x_coarse(coarse_ind{i}); tau_stokes_y_coarse(coarse_ind{i})];
    lambda_p2 = [tau_stokes_x_coarse(coarse_ind{p2}); tau_stokes_y_coarse(coarse_ind{p2})];

    % Keep per-particle image source locations for final assembly.
    rimage_k_chunks{i}{end+1,1} = rimage_vec{i,p2};
    rimage_k_chunks{p2}{end+1,1} = rimage_vec{p2,i};

    % Local ordering in pair source vector.
    im_nr = length(rimage_vec{i,p2});
    s_ind1_x = 1:N_f;
    s_ind2_x = N_f+im_nr+1:2*N_f+im_nr;
    s_ind1_y = 2*N_f+2*im_nr+1:3*N_f+2*im_nr;
    s_ind2_y = 3*N_f+3*im_nr+1:4*N_f+3*im_nr;
    e_ind1_x = N_f+1:N_f+im_nr;
    e_ind2_x = 2*N_f+im_nr+1:2*N_f+2*im_nr;
    e_ind1_y = 3*N_f+2*im_nr+1:3*N_f+3*im_nr;
    e_ind2_y = 4*N_f+3*im_nr+1:4*N_f+4*im_nr;

    if ~precomp
        rout_fine_other2 = getFineOther(opt.a_f,opt.N_f,opt.rad,refine,q,i,p2);
        [u2,v2] = StokesletDirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
            real(rout_fine_other2),imag(rout_fine_other2),...
            lambda_i(1:N_coarse),lambda_i(N_coarse+1:2*N_coarse),N_coarse);
        R2 = -[u2; v2];

        rout_fine_other1 = getFineOther(opt.a_f,opt.N_f,opt.rad,refine,q,p2,i);
        [u1,v1] = StokesletDirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
            real(rout_fine_other1),imag(rout_fine_other1),...
            lambda_p2(1:N_coarse),lambda_p2(N_coarse+1:2*N_coarse),N_coarse);
        R1 = -[u1; v1];

        coarse_to_fine_tot = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)];
    else
        coarse_to_fine_tot = [lambda_i(1:end/2); lambda_p2(1:end/2); ...
            lambda_i(end/2+1:end); lambda_p2(end/2+1:end)];
    end

    pair_mapped = Upf{i,p2}*coarse_to_fine_tot;
    beta_tot = Ypf{i,p2}*pair_mapped;

    %% Evaluate flow field on the pair itself.
    rout_pair = [rvec_out((i-1)*N_large+1:i*N_large,:); ...
        rvec_out((p2-1)*N_large+1:p2*N_large,:)];
    rin_pair = [rbase_in_f+q(i); rimage_vec{i,p2}; rbase_in_f+q(p2); rimage_vec{p2,i}];
    st_all = length(rin_pair);
    [u1,v1] = stokSLPdirect(real(rin_pair),imag(rin_pair),...
        real(rout_pair),imag(rout_pair),...
        beta_tot(1:st_all),beta_tot(st_all+1:2*st_all),st_all);
    u_pair = [u1; v1];

    ind1x = (i-1)*N_large+1:i*N_large;
    ind2x = (p2-1)*N_large+1:p2*N_large;
    ind1y = (i-1)*N_large+PM+1:i*N_large+PM;
    ind2y = (p2-1)*N_large+PM+1:p2*N_large+PM;
    pair_ind = [ind1x ind2x ind1y ind2y]';

    u_corr(pair_ind) = u_corr(pair_ind)+u_pair;

    % Subtract self-interaction on the neighbour block.
    N2 = stokSLPmat(rbase_in_c+q(i),rvec_out((p2-1)*N_large+1:p2*N_large,:),1);
    u2 = N2*lambda_i;

    N1 = stokSLPmat(rbase_in_c+q(p2),rvec_out((i-1)*N_large+1:i*N_large,:),1);
    u1 = N1*lambda_p2;

    u_corr(ind1x) = u_corr(ind1x)+u1(1:end/2);
    u_corr(ind2x) = u_corr(ind2x)+u2(1:end/2);
    u_corr(ind1y) = u_corr(ind1y)+u1(end/2+1:end);
    u_corr(ind2y) = u_corr(ind2y)+u2(end/2+1:end);

    % Store fine strengths for later global evaluation.
    pair_idx = [i; p2];
    pair_f_x = {beta_tot(s_ind1_x), beta_tot(s_ind2_x)};
    pair_f_y = {beta_tot(s_ind1_y), beta_tot(s_ind2_y)};
    pair_e_x = {beta_tot(e_ind1_x), beta_tot(e_ind2_x)};
    pair_e_y = {beta_tot(e_ind1_y), beta_tot(e_ind2_y)};

    for pair_it = 1:2
        idx = pair_idx(pair_it);

        if isempty(tau_stokes_f_x{idx})
            tau_stokes_f_x{idx} = pair_f_x{pair_it};
            tau_stokes_f_y{idx} = pair_f_y{pair_it};
        else
            tau_stokes_f_x{idx} = tau_stokes_f_x{idx} + pair_f_x{pair_it};
            tau_stokes_f_y{idx} = tau_stokes_f_y{idx} + pair_f_y{pair_it};
        end

        tau_stokes_e_x_chunks{idx}{end+1,1} = pair_e_x{pair_it};
        tau_stokes_e_y_chunks{idx}{end+1,1} = pair_e_y{pair_it};
    end
end

%% Recover per-particle extra source arrays from chunked storage.
for k = 1:P
    if isempty(tau_stokes_f_x{k})
        tau_stokes_f_x{k} = zeros(0,1);
        tau_stokes_f_y{k} = zeros(0,1);
    end

    if isempty(rimage_k_chunks{k})
        rimage_k{k} = zeros(0,1);
    else
        rimage_k{k} = vertcat(rimage_k_chunks{k}{:});
    end

    if isempty(tau_stokes_e_x_chunks{k})
        tau_stokes_e_x{k} = zeros(0,1);
        tau_stokes_e_y{k} = zeros(0,1);
    else
        tau_stokes_e_x{k} = vertcat(tau_stokes_e_x_chunks{k}{:});
        tau_stokes_e_y{k} = vertcat(tau_stokes_e_y_chunks{k}{:});
    end
end

%% Collect all source locations and strengths with one-time concatenation.
has_neigh = sort(unique(pairs(:)));
n_has = numel(has_neigh);

x_chunks = cell(1+2*n_has,1);
y_chunks = cell(1+2*n_has,1);
rvec_chunks = cell(1+2*n_has,1);

x_chunks{1} = tau_stokes_x_coarse;
y_chunks{1} = tau_stokes_y_coarse;
rvec_chunks{1} = rvec_in_coarse;

for i = 1:n_has
    k = has_neigh(i);
    x_chunks{2*i} = tau_stokes_f_x{k};
    x_chunks{2*i+1} = tau_stokes_e_x{k};
    y_chunks{2*i} = tau_stokes_f_y{k};
    y_chunks{2*i+1} = tau_stokes_e_y{k};
    rvec_chunks{2*i} = q(k)+rbase_in_f;
    rvec_chunks{2*i+1} = rimage_k{k};
end

tau_stokes_x = vertcat(x_chunks{:});
tau_stokes_y = vertcat(y_chunks{:});
rvec_in = vertcat(rvec_chunks{:});

end
