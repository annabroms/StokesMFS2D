function res = matvec_res_Stokes(tau,rvec_in,rvec_out,rimage,nimage,q,Uii,Yii,pairs,s)
%MATVEC_RES_STOKES applies the matrix corresponding to one-body
%preconditioning for the 2D Stokes resistance problem using a grid of
%coarse sources in rvec_in enhance with images stored in rimage.

P = size(q,1);
N_small = size(rvec_in,1)/P;
PM = length(rvec_out);
mu = 1; 

%Transform density mu -> lambda
[tau_stokes_x,tau_stokes_y,rot,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(tau,rimage,Uii,Yii,P,N_small,PM,pairs,s);

if isempty(rimage)
    images = 0;
else
    images = 1; 
end

%get velocity field from all source types
res = getVelocityField(rvec_in,rvec_out,tau_stokes_x,tau_stokes_y,rimage,nimage,rot,...
    tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y);

res = res+tau; %add identity

%Need to subract off the part computed twice for the one-body interactions
start_colloc = 0; 
start_im= 0; 

if ~images
    rout = rvec_out(start_colloc+1:start_colloc+pairs(1,2));
    rin = rvec_in(1:N_small);
    Nii = stokSLPmat(rin,rout,mu);
end

for i = 1:P
    % Get sources on this particle from single layer
    if images 
        rout = rvec_out(start_colloc+1:start_colloc+pairs(i,2));
        rin = rvec_in((i-1)*N_small+1:i*N_small);
        Nii = stokSLPmat(rin,rout,mu);
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
