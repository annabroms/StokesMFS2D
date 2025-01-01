function res = matvec_mob_pairprecond_peanut(tau,rbase_in_c,rbase_in_f,rvec_in,refine,rimage_vec,nimage,opt,rvec_out,rcheck,q,Us,Ys,Lc,pairs,Ucf,Ycf,Up, Yp,Lc_pair,Lf_pair,debug)
%Resistance matvec with peanut compressions for the pair corrections

% Us, Ys SVD factors for self-interaction block
% Ucf, Ycf cell array of SVD factors for pair-corrections using fine grid
% (N_f stokeslets + image points)
% Up, Yp factors for peanut compression

P = length(q);
N_large = length(rvec_out)/P;
mu = 1; 
PM = length(rvec_out);
N_f = opt.N_f;
N_c = opt.N_c; 


%Transform coarse \mu -> coarse \lambda
[tau_stokes_x, tau_stokes_nonpx,tau_self_x, tau_beta_x,tau_stokes_y,tau_stokes_nonpy,tau_self_y,tau_beta_y,u_corr] = transform_mob_peanut(tau,rbase_in_c,...
    rbase_in_f,refine,rimage_vec,nimage,opt,rvec_out,rcheck,q,Us,Ys,Lc,pairs,Ucf,Ycf,Up,Yp,Lc_pair,Lf_pair);


[ufmm,vfmm] = stokesSLPfmm(tau_stokes_x,tau_stokes_y,real(rvec_in),imag(rvec_in),real(rcheck),imag(rcheck),...
         0,5); %the last number here determines the tolerance. Higher number: higher tolerance. 

  
%NOTE! ERRORS ACCUMULATE WITH FMM IF LARGE SOURCE STRENGTHS 

%[udirect,vdirect] = StokesletDirect(real(rvec_in),imag(rvec_in),real(rcheck),imag(rcheck),tau_stokes_x,tau_stokes_y,length(rvec_in));
% See the SLP FMM test in the fast tools folder

%res = [udirect; vdirect];
res = [ufmm; vfmm];

res = res+u_corr; %Subtraction of the contribution from the peanut compressed basis on the pair itself, 
% adding the fine representation instead

%Imposing boundary conditions, if solving
if isequal(rcheck,rvec_out)

    B = getKmat2D(rvec_out(1:N_large)-q(1),0);
    K = getKmat2D(rbase_in_f,0);
    
    %This part is already taken care of... i.e. we never add and subtract the same thing. 
    % for k= 1:P
    %     bcvec = B*K'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
    %     res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/2);
    %     res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:end);
    % end
    
    %Need contribution from every pair to the velocity vector that sets
    %boundary conditions.
    has_neigh = sort(unique(pairs(:)));
    for i = 1:length(has_neigh)
        k = has_neigh(i); 
        bcvec = B*K'*[tau_beta_x((k-1)*N_f+1:k*N_f); 
            tau_beta_y((k-1)*N_f+1:k*N_f)];
    
        res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/2);
        res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:end);

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

