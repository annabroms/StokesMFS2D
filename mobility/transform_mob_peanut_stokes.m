function [tau_stokes_x, tau_stokes_nonpx,tau_self_x, tau_beta_x, ...
    tau_stokes_y,tau_stokes_nonpy,tau_self_y,tau_beta_y,u_corr,rimage_k] = ...
    transform_mob_peanut_stokes(tau,geom,basis)
%TRANSFORM_MOB_PEANUT_STOKES Map coarse boundary data to compressed Stokes sources.
%
% Syntax:
%   [tau_stokes_x,tau_stokes_nonpx,tau_self_x,tau_beta_x,...
%    tau_stokes_y,tau_stokes_nonpy,tau_self_y,tau_beta_y,u_corr,rimage_k] = ...
%       transform_mob_peanut_stokes(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] at outer collocation points.
%   geom  - Struct with fields:
%           rbase_in_c, rbase_in_f, refine, rimage_vec, opt, rvec_out, rcheck, q, pairs.
%   basis - Struct with fields:
%           U, Y, Lc, Upf, Ypf, DC_all, YC_all, Cmap, Lc_pair.
%
% Outputs:
%   tau_stokes_x, tau_stokes_y     - Projected coarse source strengths.
%   tau_stokes_nonpx, tau_stokes_nonpy - Unprojected coarse source strengths.
%   tau_self_x, tau_self_y         - One-body coarse projected strengths.
%   tau_beta_x, tau_beta_y         - Cell arrays with per-particle pair-source strengths
%                                    ordered as [fine-body; fine-image].
%   u_corr                         - Pair-local correction on the check grid.
%   rimage_k                       - Per-particle image-source locations accumulated from pair corrections.

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
Lc = basis.Lc;
Upf = basis.Upf;
Ypf = basis.Ypf;
DC_all = basis.DC_all;
YC_all = basis.YC_all;
Cmap = basis.Cmap;
Lc_pair = basis.Lc_pair;
                                                                              
P = length(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P; 
PM = length(rvec_out);
N_check = length(rcheck_out)/P;
PM2 = length(rcheck_out);

%map densities back

%Preallocate
tau_stokes_x = zeros(N_c*P,1); 
tau_stokes_y = zeros(N_c*P,1);
tau_stokes_nonpx = zeros(N_c*P,1);
tau_stokes_nonpy = zeros(N_c*P,1);

tau_beta_f_x = zeros(N_f*P,1);
tau_beta_f_y = zeros(N_f*P,1);
tau_beta_e_x_chunks = repmat({cell(0,1)},P,1);
tau_beta_e_y_chunks = repmat({cell(0,1)},P,1);

rimage_k = cell(P,1);


u_corr = zeros(2*N_check*P,1); 
precomp = opt.precomp;  
use_matrix_free_projection = true; % set false to use the original K-based projector

% Phase 1 
for i = 1:P

    coarse_ind = (i-1)*N_c+1:i*N_c;  

    %Retrieve self evaluation blocks    
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    %U{i}, Y{i}
    step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other?
    tau_mapped = Y{1}*step1; %this is the mapped density for this particle to throw in to the kernel

    tau_stokes_nonpx(coarse_ind) = tau_mapped(1:N_c);
    tau_stokes_nonpy(coarse_ind) = tau_mapped(N_c+1:end);

    %Project
    tau_i_x = tau_mapped(1:N_c);
    tau_i_y = tau_mapped(N_c+1:end);
    tau_mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y]; 

    tau_stokes_x(coarse_ind) = tau_mapped(1:N_c);
    tau_stokes_y(coarse_ind) = tau_mapped(N_c+1:end);

end


%For subtracting off self-contribution in matvec
tau_self_x = tau_stokes_x;
tau_self_y = tau_stokes_y; 

s_ind1_x = 1:opt.N_c; 
s_ind2_x = opt.N_c+1:2*opt.N_c;
s_ind1_y = 2*opt.N_c+1:3*opt.N_c;
s_ind2_y = 3*opt.N_c+1:4*opt.N_c;

% Phase 2
for i = 1:P

    % Reuse phase-1 coarse projected mapping (lambda_coarse) for this body.
    coarse_ind_i = (i-1)*N_c+1:i*N_c;
    lambda_coarse_i = [tau_self_x(coarse_ind_i); tau_self_y(coarse_ind_i)];

    %check if particle is in pair.
    if ~isempty(pairs)
        %neigh = find(pairs(:,1)==i);
        [I,J] = find(pairs==i);
        neigh = [I J];
    else
        neigh = [];
    end
    
    if ~isempty(neigh) %if isempty(neigh) has no neighbours 
        
        for k = 1:size(neigh,1)

            if neigh(k,2) == 2 %has already been counted in 
                break
            end
            %send in two particles with the fine grid 
            p2 = pairs(neigh(k),2); 

            rimage_k{i} = [rimage_k{i}; rimage_vec{i,p2}];
            rimage_k{p2} = [rimage_k{p2}; rimage_vec{p2,i}];

            coarse_ind_p2 = (p2-1)*N_c+1:p2*N_c;
            lambda_coarse_p2 = [tau_self_x(coarse_ind_p2); tau_self_y(coarse_ind_p2)];

            %% First, determine beta, the fine sources
            if ~precomp                
                
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,refine,q,i,p2); 
                %Nother = singleLayer(rbase_in_c+q(i),rout_fine_other,mu);
                %R2 = -Nother*tau_mapped; %read off on particle 2
    
                [u2,v2] = stokesletDirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
                    real(rout_fine_other),imag(rout_fine_other),...
                    lambda_coarse_i(1:N_c),lambda_coarse_i(N_c+1:2*N_c),N_c);
                R2  = -[u2; v2];
    
                %% Do a similar thing for the other order of the particles in the pair
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,refine,q,p2,i); 
                %Nother2 = singleLayer(rbase_in_c+q(p2),rout_fine_other,mu);
                %R1 = -Nother2*mapped; %read off on particle 1
    
                %To be replaced with C implementation?
                [u1,v1] = stokesletDirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
                    real(rout_fine_other),imag(rout_fine_other),...
                    lambda_coarse_p2(1:N_c),lambda_coarse_p2(N_c+1:2*N_c),N_c);
                R1 = -[u1; v1];

            
                rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)]; 

            else              
                rhs = [lambda_coarse_i(1:end/2); lambda_coarse_p2(1:end/2); ...
                    lambda_coarse_i(end/2+1:end); lambda_coarse_p2(end/2+1:end)];
            end
   
            pair_mapped = Upf{i,p2}*rhs; 
            tau_mapped_tot = Ypf{i,p2}*pair_mapped;

            %now we have beta. Lets do peanut compression from
            %here
            
            if opt.cmap
                tau_peanut_ntot = Cmap{i,p2}*rhs;
            else
                tau_peanut_temp = DC_all{i,p2}*tau_mapped_tot;
                tau_peanut_ntot = YC_all{i,p2}*tau_peanut_temp;
            end

            tau_peanut_tot= Lc_pair*tau_peanut_ntot; 

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

            %Store
            tau_stokes_x((i-1)*N_c+1:N_c*i) = tau_stokes_x((i-1)*N_c+1:N_c*i)+...
                    tau_peanut_tot(s_ind1_x);
            tau_stokes_y((i-1)*N_c+1:N_c*i) = tau_stokes_y((i-1)*N_c+1:N_c*i)+...
                tau_peanut_tot(s_ind1_y);

            tau_stokes_x((p2-1)*N_c+1:N_c*p2) = tau_stokes_x((p2-1)*N_c+1:N_c*p2)+...
                tau_peanut_tot(s_ind2_x);
            tau_stokes_y((p2-1)*N_c+1:N_c*p2) = tau_stokes_y((p2-1)*N_c+1:N_c*p2)+...
                tau_peanut_tot(s_ind2_y);

            % Store pair-source strengths for BC correction:
            % - body-fine (f) terms are summed on the shared fine grid
            % - image-fine (e) terms are concatenated per pair
            tau_beta_f_x((i-1)*N_f+1:N_f*i) = tau_beta_f_x((i-1)*N_f+1:N_f*i)+tau_mapped_tot(f_ind1_x);
            tau_beta_f_x((p2-1)*N_f+1:N_f*p2) = tau_beta_f_x((p2-1)*N_f+1:N_f*p2)+tau_mapped_tot(f_ind2_x);
            tau_beta_f_y((i-1)*N_f+1:N_f*i) = tau_beta_f_y((i-1)*N_f+1:N_f*i)+tau_mapped_tot(f_ind1_y);
            tau_beta_f_y((p2-1)*N_f+1:N_f*p2) = tau_beta_f_y((p2-1)*N_f+1:N_f*p2)+tau_mapped_tot(f_ind2_y);

            tau_beta_e_x_chunks{i}{end+1,1} = tau_mapped_tot(e_ind1_x);
            tau_beta_e_x_chunks{p2}{end+1,1} = tau_mapped_tot(e_ind2_x);
            tau_beta_e_y_chunks{i}{end+1,1} = tau_mapped_tot(e_ind1_y);
            tau_beta_e_y_chunks{p2}{end+1,1} = tau_mapped_tot(e_ind2_y);

            %% Evaluate flow field on pair and subtract this contribution. 
            % Replace with fine grid on the pair
            rout_pair = [rcheck_out((i-1)*N_check+1:i*N_check,:); rcheck_out((p2-1)*N_check+1:p2*N_check,:)];

            rin_pair_f = [rbase_in_f+q(i); rimage_vec{i,p2}; rbase_in_f+q(p2); rimage_vec{p2,i}];
            N_im = length(rin_pair_f);

            % Split beta into body-wise vectors and project each body separately
            % (same strategy as getMobPairTransformationStokesCached).
            tau_fine_i = [tau_mapped_tot(f_ind1_x); tau_mapped_tot(e_ind1_x); ...
                tau_mapped_tot(f_ind1_y); tau_mapped_tot(e_ind1_y)];
            tau_fine_p2 = [tau_mapped_tot(f_ind2_x); tau_mapped_tot(e_ind2_x); ...
                tau_mapped_tot(f_ind2_y); tau_mapped_tot(e_ind2_y)];

            rin_i = [rbase_in_f+q(i); rimage_vec{i,p2}];
            rin_p2 = [rbase_in_f+q(p2); rimage_vec{p2,i}];

            if use_matrix_free_projection
                tau_mapped_proj_i = projectOutRigid2D(tau_fine_i,rin_i,q(i));
                tau_mapped_proj_p2 = projectOutRigid2D(tau_fine_p2,rin_p2,q(p2));
            else
                Kf1 = getKmat2D(rin_i,q(i));
                Kf2 = getKmat2D(rin_p2,q(p2));
                Lf1 = Kf1*((Kf1'*Kf1)\Kf1');
                Lf2 = Kf2*((Kf2'*Kf2)\Kf2');
                tau_mapped_proj_i = tau_fine_i-Lf1*tau_fine_i;
                tau_mapped_proj_p2 = tau_fine_p2-Lf2*tau_fine_p2;
            end

            nsrc_i = N_f + im_nr;
            nsrc_p2 = N_f + im_nr;
            tau_mapped_proj = [tau_mapped_proj_i(1:nsrc_i); tau_mapped_proj_p2(1:nsrc_p2); ...
                tau_mapped_proj_i(nsrc_i+1:end); tau_mapped_proj_p2(nsrc_p2+1:end)];

            [u1,v1] = stokesletDirect(real(rin_pair_f),imag(rin_pair_f),...
                real(rout_pair),imag(rout_pair),...
                tau_mapped_proj(1:N_im),tau_mapped_proj(N_im+1:2*N_im),N_im);
            u_pair = [u1; v1];

            %Determine contribution to be subtracted

            rin_pair_c = [rbase_in_c+q(i); rbase_in_c+q(p2)];               

            [u1,v1] = stokesletDirect(real(rin_pair_c),imag(rin_pair_c),...
                real(rout_pair),imag(rout_pair),...
                tau_peanut_tot(1:2*N_c),tau_peanut_tot(2*N_c+1:4*N_c),2*N_c);
            u_peanut_corr = [u1; v1];
 
            pair_ind = [(i-1)*N_check+1:i*N_check (p2-1)*N_check+1:p2*N_check ...
                (i-1)*N_check+PM2+1:i*N_check+PM2 ...
                (p2-1)*N_check+PM2+1:p2*N_check+PM2]'; 
            u_corr(pair_ind) = u_corr(pair_ind)+u_pair-u_peanut_corr;


            
        end
        


    end
       
end

% Assemble per-particle pair-source vectors [f; e] for BC correction.
tau_beta_x = cell(P,1);
tau_beta_y = cell(P,1);
for k = 1:P
    fine_ind = (k-1)*N_f+1:k*N_f;
    if isempty(tau_beta_e_x_chunks{k})
        tau_beta_x{k} = tau_beta_f_x(fine_ind);
        tau_beta_y{k} = tau_beta_f_y(fine_ind);
    else
        tau_beta_x{k} = [tau_beta_f_x(fine_ind); vertcat(tau_beta_e_x_chunks{k}{:})];
        tau_beta_y{k} = [tau_beta_f_y(fine_ind); vertcat(tau_beta_e_y_chunks{k}{:})];
    end
end

end
