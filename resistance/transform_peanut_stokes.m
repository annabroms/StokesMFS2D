function [tau_stokes_x, tau_self_x, tau_beta_x, tau_cf_x, ...
    tau_stokes_y, tau_self_y, tau_beta_y, tau_cf_y, u_corr] = transform_peanut_stokes(tau,geom,basis)
%TRANSFORM_PEANUT_STOKES Map coarse boundary data to coarse and fine Stokes source strengths.
%
% Syntax:
%   [tau_stokes_x,tau_self_x,tau_beta_x,tau_cf_x,...
%    tau_stokes_y,tau_self_y,tau_beta_y,tau_cf_y,u_corr] = ...
%       transform_peanut_stokes(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] at outer collocation points.
%   geom  - Struct with fields:
%           rbase_in_c, rbase_in_f, refine, rimage_vec, opt, rvec_out, rcheck, q, pairs.
%   basis - Struct with fields:
%           U, Y, Upf, Ypf, DC_all, YC_all, Cmap.
%
% Outputs:
%   tau_stokes_x - Coarse x-source strengths used in the compressed representation.
%   tau_self_x   - One-body coarse x-source strengths (used for self-block subtraction).
%   tau_beta_x   - Fine x-source strengths from pair-local solves (for postprocessing).
%   tau_cf_x     - Reserved coarse x-storage for alternative compressed maps.
%   tau_stokes_y - Coarse y-source strengths used in the compressed representation.
%   tau_self_y   - One-body coarse y-source strengths (used for self-block subtraction).
%   tau_beta_y   - Fine y-source strengths from pair-local solves (for postprocessing).
%   tau_cf_y     - Reserved coarse y-storage for alternative compressed maps.
%   u_corr       - Pair correction that replaces compressed pair blocks by fine pair evaluation.

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
refine = geom.refine;
rimage_vec = geom.rimage_vec;
opt = geom.opt;
rvec_out = geom.rvec_out;
rcheck_out = geom.rcheck;
q = geom.q;
pairs = geom.pairs;

U = basis.U;
Y = basis.Y;
Upf = basis.Upf;
Ypf = basis.Ypf;
DC_all = basis.DC_all;
YC_all = basis.YC_all;
Cmap = basis.Cmap;

P = length(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P;
PM = length(rvec_out);
N_check = length(rcheck_out)/P;
PM2 = length(rcheck_out);

precomp = opt.precomp;

% Preallocate coarse/fine source storage.
tau_stokes_x = zeros(N_c*P,1);
tau_stokes_y = zeros(N_c*P,1);
tau_beta_x = zeros(N_f*P,1);
tau_beta_y = zeros(N_f*P,1);
tau_cf_x = zeros(N_c*P,1);
tau_cf_y = zeros(N_c*P,1);

% Pair correction evaluated on rcheck_out.
u_corr = zeros(2*N_check*P,1);

% Phase 1: map all bodies once from boundary data to coarse source strengths.
for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c;

    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    step1 = U{1}*[tau_particle_x; tau_particle_y];
    lambda_coarse_i = Y{1}*step1;

    tau_stokes_x(coarse_ind) = lambda_coarse_i(1:N_c);
    tau_stokes_y(coarse_ind) = lambda_coarse_i(N_c+1:2*N_c);
end

% One-body part is needed when correcting diagonal blocks in the matvec.
tau_self_x = tau_stokes_x;
tau_self_y = tau_stokes_y;

if isempty(pairs)
    return;
end

% Phase 2: process each pair once.
for pair_row = 1:size(pairs,1)
    i = pairs(pair_row,1);
    p2 = pairs(pair_row,2);

    coarse_i = (i-1)*N_c+1:i*N_c;
    coarse_p2 = (p2-1)*N_c+1:p2*N_c;

    lambda_i = [tau_self_x(coarse_i); tau_self_y(coarse_i)];
    lambda_p2 = [tau_self_x(coarse_p2); tau_self_y(coarse_p2)];

    if precomp
        rhs = [lambda_i(1:N_c); lambda_p2(1:N_c); ...
               lambda_i(N_c+1:2*N_c); lambda_p2(N_c+1:2*N_c)];
    else
        rout_fine_other2 = getFineOther(opt.a_f,opt.N_f,refine,q,i,p2);
        [u2,v2] = stokesletDirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
            real(rout_fine_other2),imag(rout_fine_other2),...
            lambda_i(1:N_c),lambda_i(N_c+1:2*N_c),N_c);
        R2 = -[u2; v2];

        rout_fine_other1 = getFineOther(opt.a_f,opt.N_f,refine,q,p2,i);
        [u1,v1] = stokesletDirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
            real(rout_fine_other1),imag(rout_fine_other1),...
            lambda_p2(1:N_c),lambda_p2(N_c+1:2*N_c),N_c);
        R1 = -[u1; v1];

        rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)];
    end

    pair_mapped = Upf{i,p2}*rhs;
    tau_mapped_tot = Ypf{i,p2}*pair_mapped;

    if opt.cmap
        tau_peanut_tot = Cmap{i,p2}*rhs;
    else
        tau_peanut_temp = DC_all{i,p2}*tau_mapped_tot;
        tau_peanut_tot = YC_all{i,p2}*tau_peanut_temp;
    end

    % Store fine source strengths for postprocessing.
    fine_i = (i-1)*N_f+1:i*N_f;
    fine_p2 = (p2-1)*N_f+1:p2*N_f;
    tau_beta_x(fine_i) = tau_beta_x(fine_i) + tau_mapped_tot(1:N_f);
    tau_beta_x(fine_p2) = tau_beta_x(fine_p2) + tau_mapped_tot(N_f+1:2*N_f);
    tau_beta_y(fine_i) = tau_beta_y(fine_i) + tau_mapped_tot(2*N_f+1:3*N_f);
    tau_beta_y(fine_p2) = tau_beta_y(fine_p2) + tau_mapped_tot(3*N_f+1:4*N_f);

    % Evaluate pair correction on check grid.
    rout_pair = [rcheck_out((i-1)*N_check+1:i*N_check,:); ...
                 rcheck_out((p2-1)*N_check+1:p2*N_check,:)];

    rin_pair_f = [rbase_in_f+q(i); rimage_vec{i,p2}; rbase_in_f+q(p2); rimage_vec{p2,i}];
    n_fpair = length(rin_pair_f);
    [u1,v1] = stokesletDirect(real(rin_pair_f),imag(rin_pair_f),...
        real(rout_pair),imag(rout_pair),...
        tau_mapped_tot(1:n_fpair),tau_mapped_tot(n_fpair+1:2*n_fpair),n_fpair);
    u_pair = [u1; v1];

    rin_pair_c = [rbase_in_c+q(i); rbase_in_c+q(p2)];
    [u1,v1] = stokesletDirect(real(rin_pair_c),imag(rin_pair_c),...
        real(rout_pair),imag(rout_pair),...
        tau_peanut_tot(1:2*N_c),tau_peanut_tot(2*N_c+1:4*N_c),2*N_c);
    u_peanut_corr = [u1; v1];

    pair_ind = [(i-1)*N_check+1:i*N_check ...
                (p2-1)*N_check+1:p2*N_check ...
                (i-1)*N_check+PM2+1:i*N_check+PM2 ...
                (p2-1)*N_check+PM2+1:p2*N_check+PM2]';
    u_corr(pair_ind) = u_corr(pair_ind) + u_pair - u_peanut_corr;

    % Add pair-compressed coarse strengths.
    tau_stokes_x(coarse_i) = tau_stokes_x(coarse_i) + tau_peanut_tot(1:N_c);
    tau_stokes_x(coarse_p2) = tau_stokes_x(coarse_p2) + tau_peanut_tot(N_c+1:2*N_c);
    tau_stokes_y(coarse_i) = tau_stokes_y(coarse_i) + tau_peanut_tot(2*N_c+1:3*N_c);
    tau_stokes_y(coarse_p2) = tau_stokes_y(coarse_p2) + tau_peanut_tot(3*N_c+1:4*N_c);
end

end
