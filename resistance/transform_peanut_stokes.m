function [lam_c_x, lam_self_x, lam_f_x, ...
    lam_c_y, lam_self_y, lam_f_y, u_corr, rimage_k] = transform_peanut_stokes(tau,geom,basis)
%TRANSFORM_PEANUT_STOKES Map coarse boundary data to interior coarse and fine Stokes source strengths.
%
% Syntax:
%   [lam_c_x,lam_self_x,lam_f_x,...
%    lam_c_y,lam_self_y,lam_f_y,u_corr,rimage_k] = ...
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
%   lam_c_x - Coarse x-source strengths used in the compressed representation.
%   lam_self_x   - One-body coarse x-source strengths (used for self-block subtraction).
%   lam_f_x   - Cell array with per-particle pair-source strengths in x,
%                  ordered as [fine-body; fine-image].
%   lam_c_y - Coarse y-source strengths used in the compressed representation.
%   lam_self_y   - One-body coarse y-source strengths (used for self-block subtraction).
%   lam_f_y   - Cell array with per-particle pair-source strengths in y,
%                  ordered as [fine-body; fine-image].
%   u_corr       - Pair correction that replaces compressed pair blocks by fine pair evaluation.
%   rimage_k     - Per-particle concatenated image/source nodes from all close pairs.

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
lam_c_x = zeros(N_c*P,1);
lam_c_y = zeros(N_c*P,1);
lam_beta_f_x = zeros(N_f*P,1);
lam_beta_f_y = zeros(N_f*P,1);
lam_beta_e_x_chunks = repmat({cell(0,1)},P,1);
lam_beta_e_y_chunks = repmat({cell(0,1)},P,1);
rimage_k = cell(P,1);

% Pair correction evaluated on rcheck_out.
u_corr = zeros(2*N_check*P,1);

% Phase 1: map all bodies once from boundary data to coarse source strengths.
for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c;

    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    step1 = U{1}*[tau_particle_x; tau_particle_y];
    lambda_coarse_i = Y{1}*step1;

    lam_c_x(coarse_ind) = lambda_coarse_i(1:N_c);
    lam_c_y(coarse_ind) = lambda_coarse_i(N_c+1:2*N_c);
end

% One-body part is needed when correcting diagonal blocks in the matvec.
lam_self_x = lam_c_x;
lam_self_y = lam_c_y;

% Phase 2: process each pair once.
for pair_row = 1:size(pairs,1)
    i = pairs(pair_row,1);
    p2 = pairs(pair_row,2);

    coarse_i = (i-1)*N_c+1:i*N_c;
    coarse_p2 = (p2-1)*N_c+1:p2*N_c;

    lambda_i = [lam_self_x(coarse_i); lam_self_y(coarse_i)];
    lambda_p2 = [lam_self_x(coarse_p2); lam_self_y(coarse_p2)];
    rimage_k{i} = [rimage_k{i}; rimage_vec{i,p2}];
    rimage_k{p2} = [rimage_k{p2}; rimage_vec{p2,i}];

    if precomp
        rhs = [lambda_i(1:N_c); lambda_p2(1:N_c); ...
               lambda_i(N_c+1:2*N_c); lambda_p2(N_c+1:2*N_c)];
    else
        rout_fine_other2 = getFineOther(opt.a_f,opt.N_f,refine,q,i,p2);
        [u2,v2] = stokSLPdirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
            real(rout_fine_other2),imag(rout_fine_other2),...
            lambda_i(1:N_c),lambda_i(N_c+1:2*N_c),N_c);
        R2 = -[u2; v2];

        rout_fine_other1 = getFineOther(opt.a_f,opt.N_f,refine,q,p2,i);
        [u1,v1] = stokSLPdirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
            real(rout_fine_other1),imag(rout_fine_other1),...
            lambda_p2(1:N_c),lambda_p2(N_c+1:2*N_c),N_c);
        R1 = -[u1; v1];

        rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)];
    end

    % 1) Compute peanut compressed coarse source strengths
    % 2) Determine flow field on pair itself due to fine representation.
    % With compression matrices concatenated (opt.cmap = 1), we can do this using the BVP
    % rhs for the fine solve

    rcheck_i = rcheck_out((i-1)*N_check+1:i*N_check,:);
    rcheck_p2 = rcheck_out((p2-1)*N_check+1:p2*N_check,:);
    rout_pair = [rcheck_i; rcheck_p2];
    rin_ci = rbase_in_c+q(i);
    rin_cp2 = rbase_in_c+q(p2);
    rin_pair_c = [rin_ci; rin_cp2];

    if opt.cmap
        tau_peanut_tot = Cmap{i,p2}*rhs;
        [ui,vi] = stokSLPdirect(real(rin_cp2),imag(rin_cp2),...
            real(rcheck_i),imag(rcheck_i), lam_self_x(coarse_p2),...
            lam_self_y(coarse_p2),N_c);
        [up2,vp2] = stokSLPdirect(real(rin_ci),imag(rin_ci),...
            real(rcheck_p2),imag(rcheck_p2),lam_self_x(coarse_i),...
            lam_self_y(coarse_i),N_c);
        u_fine = -[ui; up2; vi; vp2];
    else
        pair_mapped = Upf{i,p2}*rhs;
        tau_mapped_tot = Ypf{i,p2}*pair_mapped;
        tau_peanut_temp = DC_all{i,p2}*tau_mapped_tot;
        tau_peanut_tot = YC_all{i,p2}*tau_peanut_temp;

        % Pair-local indexing in tau_mapped_tot:
        % [f_i_x; e_i_x; f_p2_x; e_p2_x; f_i_y; e_i_y; f_p2_y; e_p2_y]
        im_nr = length(rimage_vec{i,p2});
        f_ind1_x = 1:N_f;
        e_ind1_x = N_f+1:N_f+im_nr;
        f_ind2_x = N_f+im_nr+1:2*N_f+im_nr;
        e_ind2_x = 2*N_f+im_nr+1:2*N_f+2*im_nr;
        f_ind1_y = 2*N_f+2*im_nr+1:3*N_f+2*im_nr;
        e_ind1_y = 3*N_f+2*im_nr+1:3*N_f+3*im_nr;
        f_ind2_y = 3*N_f+3*im_nr+1:4*N_f+3*im_nr;
        e_ind2_y = 4*N_f+3*im_nr+1:4*N_f+4*im_nr;
    
        % Store fine/body and fine/image strengths for postprocessing.
        fine_i = (i-1)*N_f+1:i*N_f;
        fine_p2 = (p2-1)*N_f+1:p2*N_f;
        lam_beta_f_x(fine_i) = lam_beta_f_x(fine_i) + tau_mapped_tot(f_ind1_x);
        lam_beta_f_x(fine_p2) = lam_beta_f_x(fine_p2) + tau_mapped_tot(f_ind2_x);
        lam_beta_f_y(fine_i) = lam_beta_f_y(fine_i) + tau_mapped_tot(f_ind1_y);
        lam_beta_f_y(fine_p2) = lam_beta_f_y(fine_p2) + tau_mapped_tot(f_ind2_y);

        lam_beta_e_x_chunks{i}{end+1,1} = tau_mapped_tot(e_ind1_x);
        lam_beta_e_x_chunks{p2}{end+1,1} = tau_mapped_tot(e_ind2_x);
        lam_beta_e_y_chunks{i}{end+1,1} = tau_mapped_tot(e_ind1_y);
        lam_beta_e_y_chunks{p2}{end+1,1} = tau_mapped_tot(e_ind2_y);


         % Evaluate fine grid pair correction on check grid.
        rin_pair_f = [rbase_in_f+q(i); rimage_vec{i,p2}; rbase_in_f+q(p2); rimage_vec{p2,i}];
        n_fpair = length(rin_pair_f);
        [u1,v1] = stokSLPdirect(real(rin_pair_f),imag(rin_pair_f),...
            real(rout_pair),imag(rout_pair),...
            tau_mapped_tot(1:n_fpair),tau_mapped_tot(n_fpair+1:2*n_fpair),n_fpair);
        u_fine = [u1; v1];
   
    
    end

    % The velocity field correction for the pair happens here! 
    [u1,v1] = stokSLPdirect(real(rin_pair_c),imag(rin_pair_c),...
        real(rout_pair),imag(rout_pair),...
        tau_peanut_tot(1:2*N_c),tau_peanut_tot(2*N_c+1:4*N_c),2*N_c);
    u_peanut = [u1; v1];

    pair_ind = [(i-1)*N_check+1:i*N_check ...
                (p2-1)*N_check+1:p2*N_check ...
                (i-1)*N_check+PM2+1:i*N_check+PM2 ...
                (p2-1)*N_check+PM2+1:p2*N_check+PM2]';
    u_corr(pair_ind) = u_corr(pair_ind) + u_fine - u_peanut;

   

    % Add pair-compressed coarse strengths.
    lam_c_x(coarse_i) = lam_c_x(coarse_i) + tau_peanut_tot(1:N_c);
    lam_c_x(coarse_p2) = lam_c_x(coarse_p2) + tau_peanut_tot(N_c+1:2*N_c);
    lam_c_y(coarse_i) = lam_c_y(coarse_i) + tau_peanut_tot(2*N_c+1:3*N_c);
    lam_c_y(coarse_p2) = lam_c_y(coarse_p2) + tau_peanut_tot(3*N_c+1:4*N_c);
end


if ~opt.cmap
    % Assemble per-particle pair-source vectors [f; e]. Only needed if the
    % effective mapping Cmap_FU is not in use
    lam_f_x = cell(P,1);
    lam_f_y = cell(P,1);
    for k = 1:P
        fine_ind = (k-1)*N_f+1:k*N_f;
        if isempty(lam_beta_e_x_chunks{k})
            lam_f_x{k} = lam_beta_f_x(fine_ind);
            lam_f_y{k} = lam_beta_f_y(fine_ind);
        else
            lam_f_x{k} = [lam_beta_f_x(fine_ind); vertcat(lam_beta_e_x_chunks{k}{:})];
            lam_f_y{k} = [lam_beta_f_y(fine_ind); vertcat(lam_beta_e_y_chunks{k}{:})];
        end
    end
else
    lam_f_x = [];
    lam_f_y = [];
end

end
