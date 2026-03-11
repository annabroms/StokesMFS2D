function res = matvec_res_peanut(tau,rbase_in_c,rbase_in_f,rvec_in,rbase_out_f,refine,rimage_vec,nimage,opt,rvec_out,q,Us,Ys,pairs,Ucf,Ycf,Up, Yp,Cmap,debug)
%Resistance matvec with peanut compressions for the pair corrections

% Us, Ys SVD factors for self-interaction block
% Ucf, Ycf cell array of SVD factors for pair-corrections using fine grid
% (N_f stokeslets + image points)
% Up, Yp factors for peanut compression

P = length(q);
N_large = length(rvec_out)/P;
mu = 1; 
PM = length(rvec_out);
N_c = opt.N_c;

%Transform coarse \mu -> coarse \lambda
[tau_stokes_x, tau_self_x, ~,tau_stokes_y,tau_self_y,~,u_corr] = transform_peanut(tau,rbase_in_c,...
    rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rvec_out,rvec_out,q,Us,Ys,pairs,Ucf,Ycf,Up,Yp,Cmap,debug);


%[ufmm,vfmm] = stokesSLPfmm(tau_stokes_x,tau_stokes_y,real(rvec_in),imag(rvec_in),real(rvec_out),imag(rvec_out),...
       %  0,5); %the last number here determines the tolerance. Higher number: higher tolerance. 

  
%NOTE! ERRORS ACCUMULATE WITH FMM IF LARGE SOURCE STRENGTHS 

%[udirect,vdirect] = StokesletDirect(real(rvec_in),imag(rvec_in),real(rvec_out),imag(rvec_out),tau_stokes_x,tau_stokes_y,length(rvec_out));
% See the SLP FMM test in the fast tools folder


%res = [udirect'/4/pi; vdirect'/4/pi];
%res = [ufmm; vfmm];

res = getVelocityField(rvec_in,rvec_out,tau_stokes_x,tau_stokes_y);

res = res+u_corr; %Subtraction of the contribution from the peanut compressed basis on the pair itself, 
% adding the fine representation instead


%Need to subract off the 1-body basis contribution computed twice
rout = rvec_out(1:N_large)-q(1);
Nii = singleLayer(rbase_in_c,rout,mu);

for i = 1:P
    % Get sources on this particle from single layer
        %tau_xy = [tau_stokes_x(coarse_ind{i});
        %tau_stokes_y(coarse_ind{i})]; %Probably the wrong thing to do as
        %this includes the effect of the peanut compressed pair
        %corrections. 

        coarse_ind = (i-1)*N_c+1:i*N_c;

      %  tau_xy = [tau_self_x(coarse_ind{i}); tau_self_y(coarse_ind{i})];
        tau_xy = [tau_self_x(coarse_ind); tau_self_y(coarse_ind)];
        
        uii = Nii*tau_xy;

        %subract contribution in x
        res((i-1)*N_large+1:i*N_large) = res((i-1)*N_large+1:i*N_large)-uii(1:end/2)+tau((i-1)*N_large+1:i*N_large);

        %subract contribution in y
        res((i-1)*N_large+1+PM:i*N_large+PM) = res((i-1)*N_large+1+PM:i*N_large+PM)-...
            uii(end/2+1:end)+tau((i-1)*N_large+PM+1:i*N_large+PM);
    
end



 




end
