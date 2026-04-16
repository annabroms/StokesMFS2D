function res = matvec_mob_2B_enhanced(tau,geom,basis)
%MATVEC_MOB_2B_ENHANCED Matrix-vector action for mobility 2-body preconditioner.
%
% Syntax:
%   res = matvec_mob_2B_enhanced(tau,geom,basis)
%
% Inputs:
%   tau   - Stacked boundary data [tau_x; tau_y] on outer collocation points.
%   geom  - Geometry/problem struct with fields:
%           rbase_in_c, rbase_in_f, opt, rvec_out, rcheck, q, pairs.
%   basis - Precomputed basis struct with fields:
%           U, Y, Lc, Upf, Ypf.
%
% Output:
%   res   - Matvec result at check/collocation points [u_x; u_y].
%

rbase_in_c = geom.rbase_in_c;
rbase_in_f = geom.rbase_in_f;
opt = geom.opt;
rvec_out = geom.rvec_out;
rcheck = geom.rcheck;
q = geom.q;
pairs = geom.pairs;
use_fmm = logical(getOptField(opt,'use_fmm',true));
is_collocation = isequal(rcheck,rvec_out);
Cmap_FU = basis.Cmap_FU;
if isfield(basis,'pair_cache')
    pair_cache = basis.pair_cache;
else
    pair_cache = struct('enabled',false);
end
use_pair_cache = isfield(pair_cache,'enabled') && pair_cache.enabled;
[rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx, tau_stokes_nonpy,tau_stokes_e_nonpx, tau_stokes_e_nonpy, rimage_k] = ...
    getMobPairTransformationStokes(tau,geom,basis);
 
P = opt.P; 
PM = length(rvec_out);
N_large = PM/P;
mu = 1; 
N_c = opt.N_c;
N_f = opt.N_f;



%% Get flow field from all source types

res = getVelocityField(rvec_in,rcheck,tau_stokes_x,tau_stokes_y,use_fmm);

% Away from the collocation grid we just evaluate the recovered source
% representation directly. The boundary-only BK' and identity terms below
% are part of the solve operator, not the field evaluation.
if ~is_collocation
    return
end

rbase_out_rel = rvec_out(1:N_large)-q(1);

    
    %This part is already taken care of... 
    % for k= 1:P
    %     bcvec = B*K'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
    %     res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/2);
    %     res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:end);
    % end
    
% Add action of Lr = BK' for close pairs
if  (~isempty(Cmap_FU) || use_pair_cache) && logical(getOptField(opt,'cmap',false))
    for row = 1:size(pairs,1)
        i = pairs(row,1);
        p2 = pairs(row,2);

        rhs_pair = [tau_stokes_x(coarse_ind{i}); tau_stokes_x(coarse_ind{p2}); ...
                    tau_stokes_y(coarse_ind{i}); tau_stokes_y(coarse_ind{p2})];

        if use_pair_cache
            pair = getStokesPairInstance(pair_cache,row);

            % rotate to frame of reference
            rhs_pair = rotatePairOrderedStokesData(rhs_pair,opt.N_c,pair.meta.phase_c,conj(pair.meta.rot));
            pair_vel = pair.group.Cmap_FU*rhs_pair;

            %apply rotation back to frame of pair
            vel_i = pair.meta.rot*(pair_vel(1) + 1i*pair_vel(2));
            vel_p2 = pair.meta.rot*(pair_vel(4) + 1i*pair_vel(5));
            pair_vel = [real(vel_i); imag(vel_i); pair_vel(3); ...
                        real(vel_p2); imag(vel_p2); pair_vel(6)];
        else
            pair_vel = Cmap_FU{i,p2}*rhs_pair;
        end

        % get boundary velocity from translational and angular velocity
        rot_x = -imag(rbase_out_rel);
        rot_y = real(rbase_out_rel);
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

        % Add BK' part for the uniformly sampled sources...
        fcx = tau_stokes_nonpx((k-1)*N_f+1+P*N_c:k*N_f+P*N_c);
        fcy = tau_stokes_nonpy((k-1)*N_f+1+P*N_c:k*N_f+P*N_c);
        
        bcvec_f = applyBKt2D(rbase_out_rel,0,rbase_in_f,0,fcx,fcy);

        %and for the enhancing extra sources
        if isempty(rimage_k{k})
            bcvec_e = zeros(2*N_large,1);
        else    
            bcvec_e = applyBKt2D(rbase_out_rel,0,rimage_k{k},q(k),...
                tau_stokes_e_nonpx{k},tau_stokes_e_nonpy{k}); 
        end

        res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec_f(1:end/2)+bcvec_e(1:end/2);
        res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec_f(end/2+1:end) + bcvec_e(end/2+1:end);

    end
end


%% Correct idenity blocks 
rout = rvec_out(1:N_large)-q(1);
rin = rbase_in_c;
Nii = stokSLPmat(rin,rout,mu);
%Nii_fine = stokSLPmat(rbase_in_f,rout,mu);

for i = 1:P
    % Get sources on this particle from single layer evaluation.
    tau_xy = [tau_stokes_x(coarse_ind{i}); tau_stokes_y(coarse_ind{i})];
    
    %Not needed to redo for all! 
    %rout = rvec_out((i-1)*N_large+1:i*N_large,:);
    %rin = rbase_in_c+q(i);

    %Nii = stokSLPmat(rin,rout,mu);
    uii = Nii*tau_xy;

    %subract contribution in x
    res((i-1)*N_large+1:i*N_large) = res((i-1)*N_large+1:i*N_large)-uii(1:end/2)+tau((i-1)*N_large+1:i*N_large);

    %subract contribution in y
    res((i-1)*N_large+1+PM:i*N_large+PM) = res((i-1)*N_large+1+PM:i*N_large+PM)-...
        uii(end/2+1:end)+tau((i-1)*N_large+PM+1:i*N_large+PM);
end

end
