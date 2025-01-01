function res = matvec_2D_Stokes(tau,rvec_in,rvec_out,rimage,nimage,q,Uii,Yii,pairs,s)
%MATVEC_2D_STOKES applies the matrix corresponding to one body
%preconditioning for the 2D Stokes resistance problem using a grid of
%coarse sources in rvec_in enhance with images stored in rimage.

P = size(q,1);
N_small = size(rvec_in,1)/P;
PM = length(rvec_out);
mu = 1; 

%Transform density mu -> lambda
[tau_stokes_x,tau_stokes_y,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(tau,rimage,Uii,Yii,P,N_small,PM,pairs,s);

if isempty(rimage)
    images = 0;
else
    images = 1; 
end

res = tau; 

%Apply all stokeslets

%For debugging... 
%NS = singleLayer(rvec_in,rvec_out,mu);
%U_proxy = NS*[tau_stokes_x; tau_stokes_y]; %I think this vector has to be ordered differently


[ufmm,vfmm] = stokesSLPfmm(tau_stokes_x,tau_stokes_y,real(rvec_in),imag(rvec_in),real(rvec_out),imag(rvec_out),...
        0,5); %the last number here determines the tolerance. Higher number: higher tolerance. 
  
%NOTE! Errors might accumulate in the FMM if density norm is large!  

%[udirect,vdirect] = StokesletDirect(real(rvec_in),imag(rvec_in),real(rvec_out),imag(rvec_out),tau_stokes_x,tau_stokes_y,length(rvec_out));
% See the SLP FMM test in the fast tools folder


if ~images
    res = res+[ufmm; vfmm];
    %res = res+[udirect'/4/pi; vdirect'/4/pi];
else

    u_stress = getStresslets(tau_stress_x,tau_stress_y,rimage,rvec_out,real(nimage),imag(nimage));
    u_pot = getPotdip(tau_pot_x, tau_pot_y,rimage,rvec_out);
   
    %res = res + [udirect'/4/pi; vdirect'/4/pi]+U_image;

    res = res + [ufmm; vfmm]+u_stress + u_pot; 
end

%Need to subract off the part computed twice for the one-body interactions
start_colloc = 0; 
start_im= 0; 

if ~images
    rout = rvec_out(start_colloc+1:start_colloc+pairs(1,2));
    rin = rvec_in(1:N_small);
    Nii = singleLayer(rin,rout,mu);
end

for i = 1:P
    % Get sources on this particle from single layer
    if images 
        rout = rvec_out(start_colloc+1:start_colloc+pairs(i,2));
        rin = rvec_in((i-1)*N_small+1:i*N_small);
        Nii = singleLayer(rin,rout,mu);
    end
    tau_xy = [tau_stokes_x((i-1)*N_small+1:i*N_small); tau_stokes_y((i-1)*N_small+1:i*N_small)];
    
    
    uii = Nii*tau_xy;
    

    %% Do the same thing with the images
    if images 
        rim = rimage(start_im+1:start_im+pairs(i,1)); 
        nim = nimage(start_im+1:start_im+pairs(i,1)); 
    
        Nim =  getImageKernels2D(rim,nim,rout,mu,s);
        stress_x = tau_stress_x(start_im+1:start_im+pairs(i,1));
        stress_y = tau_stress_y(start_im+1:start_im+pairs(i,1));
        pot_x = tau_pot_x(start_im+1:start_im+pairs(i,1));
        pot_y = tau_pot_y(start_im+1:start_im+pairs(i,1));
        tau_xy = [stress_x; stress_y; pot_x; pot_y];
        
        uii = uii+Nim*tau_xy;
        
    end

    %Take care of contribution in x
    res(start_colloc+1:start_colloc+pairs(i,2)) = res(start_colloc+1:start_colloc+pairs(i,2))-uii(1:end/2);

    %subract contribution in y
    start = start_colloc+PM;
    res(start+1:start+pairs(i,2)) = res(start+1:start+pairs(i,2))-uii(end/2+1:end);
    
    
    start_colloc = start_colloc+pairs(i,2);
    start_im = start_im+pairs(i,1); 


end


end