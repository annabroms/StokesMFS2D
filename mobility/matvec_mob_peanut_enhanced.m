function res = matvec_mob_peanut_enhanced(tau,geom,basis)
%MATVEC_MOB_PEANUT_ENHANCED Matrix-vector action for mobility peanut preconditioner.
%
% Syntax:
%   res = matvec_mob_peanut_enhanced(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] at outer collocation points.
%   geom  - Struct with fields:
%           rbase_in_c, rbase_in_f, rvec_in, opt, rvec_out, rcheck, q, pairs.
%   basis - Struct with fields:
%           U, Y, Lc, Upf, Ypf, DC_all, YC_all, Cmap, Lc_pair.
%
% Output:
%   res   - Velocity data [u_x; u_y] on the check/collocation grid.

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
rvec_in = geom.rvec_in;
opt = geom.opt;
rvec_out = geom.rvec_out;
rcheck = geom.rcheck;
q = geom.q;
pairs = geom.pairs;
if isfield(geom,'rimage_vec')
    rimage_vec = geom.rimage_vec;
else
    rimage_vec = [];
end

Cmap_FU = basis.Cmap_FU;

P = length(q);
N_large = length(rvec_out)/P;
mu = 1; 
PM = length(rvec_out);
N_c = opt.N_c; 


%Transform coarse \mu -> coarse \lambda
[tau_stokes_x, tau_stokes_nonpx,tau_self_x, tau_beta_x,tau_stokes_y,...
    tau_stokes_nonpy,tau_self_y,tau_beta_y,u_corr,rimage_k] = ...
    transform_mob_peanut_stokes(tau,geom,basis);

res = getVelocityField(rvec_in,rcheck,tau_stokes_x,tau_stokes_y);

res = res+u_corr; %Subtraction of the contribution from the peanut compressed basis on the pair itself, 
% adding the fine representation instead

%Imposing boundary conditions, if solving
if isequal(rcheck,rvec_out)

    rbase_out_rel = rvec_out(1:N_large)-q(1);
    
    %This part is already taken care of... i.e. we never add and subtract the same thing. 
    % for k= 1:P
    %     bcvec = B*K'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
    %     res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/2);
    %     res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:end);
    % end
    
    %Need contribution from every pair to the velocity vector that sets
    %boundary conditions.
    if opt.cmap
        if isempty(Cmap_FU)
            error('opt.cmap=true but basis.Cmap_FU is empty.');
        end
        if isempty(rimage_vec)
            error('opt.cmap=true requires geom.rimage_vec for pair-local BK'' corrections.');
        end
        % Loop over all pairs in contact
        for pair_it = 1:size(pairs,1)
            i = pairs(pair_it,1);
            p2 = pairs(pair_it,2);

            coarse_i = (i-1)*N_c+1:i*N_c;
            coarse_p2 = (p2-1)*N_c+1:p2*N_c;

            %Here the projected sources are used, as (I-L) is not yet
            %applied for Cmap_FU.
            rhs_pair = [tau_self_x(coarse_i); tau_self_x(coarse_p2); ...
                        tau_self_y(coarse_i); tau_self_y(coarse_p2)];

            % Determine rigid body motion for the pair, using ansatz
            pair_vel = -Cmap_FU{i,p2}*rhs_pair;

            %These velocities now must be mapped to boundary velocities at the collocation points on each particle, and added to the result.
            pair_vel = pair_vel(:); % [Ux_i; Uy_i; W_i; Ux_p2; Uy_p2; W_p2]
            rot_x = -imag(rbase_out_rel);
            rot_y =  real(rbase_out_rel);
            bvec_i = [pair_vel(1) + pair_vel(3)*rot_x; ...
                      pair_vel(2) + pair_vel(3)*rot_y];
            bvec_p2 = [pair_vel(4) + pair_vel(6)*rot_x; ...
                       pair_vel(5) + pair_vel(6)*rot_y];

            res((i-1)*N_large+1:i*N_large) = res((i-1)*N_large+1:i*N_large) + bvec_i(1:end/2);
            res((i-1)*N_large+PM+1:i*N_large+PM) = res((i-1)*N_large+PM+1:i*N_large+PM) + bvec_i(end/2+1:end);

            res((p2-1)*N_large+1:p2*N_large) = res((p2-1)*N_large+1:p2*N_large) + bvec_p2(1:end/2);
            res((p2-1)*N_large+PM+1:p2*N_large+PM) = res((p2-1)*N_large+PM+1:p2*N_large+PM) + bvec_p2(end/2+1:end);
        end
    else
        has_neigh = sort(unique(pairs(:)));
        for i = 1:length(has_neigh)
            k = has_neigh(i); 

            % tau_beta_{x,y}{k} stores all pair-source strengths on particle k
            % as [fine-body; fine-image], matching source points [rbase_in_f+q(k); rimage_k{k}].
            rsrc_k = [rbase_in_f+q(k); rimage_k{k}];
            bcvec = applyBKt2D(rbase_out_rel,0,rsrc_k,q(k),tau_beta_x{k},tau_beta_y{k});

            res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/2);
            res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:end);
        end
    end

    %Need to subract off the 1-body basis contribution computed twice
    rout = rvec_out(1:N_large)-q(1);
    Nii = singleLayer(rbase_in_c,rout,mu);
    
    for i = 1:P
        % Get contribution on this particle from 1.body basis
        coarse_ind = (i-1)*N_c+1:i*N_c;    
        tau_xy = [tau_self_x(coarse_ind); tau_self_y(coarse_ind)];
        
        uii = Nii*tau_xy; %This is G*(I-L)
    
        %subract contribution in x
        res((i-1)*N_large+1:i*N_large) = res((i-1)*N_large+1:i*N_large)-uii(1:end/2)+tau((i-1)*N_large+1:i*N_large);
    
        %subract contribution in y
        res((i-1)*N_large+1+PM:i*N_large+PM) = res((i-1)*N_large+1+PM:i*N_large+PM)-...
            uii(end/2+1:end)+tau((i-1)*N_large+PM+1:i*N_large+PM);
        
    end

end



 




end
