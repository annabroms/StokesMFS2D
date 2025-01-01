function res = matvec_2D_mobility(tau,rvec_in,rvec_out,rvec_check,rimage,nimage,q,Uii,Yii,L,pairs,s,solve,ver1)
%MATVEC_2D_STOKES applies the matrix corresponding to one body
%preconditioning for the 2D Stokes mobility problem using a grid of
%coarse sources in rvec_in enhanced with images stored in rimage.

P = size(q,1);
N_small = size(rvec_in,1)/P;
PM = length(rvec_out);
mu = 1; 



%Transform density mu -> lambda
[tau_stokes_x,tau_stokes_y,tau_stress_x,tau_stress_y,tau_pot_x,tau_pot_y] = getTransformedDensity(tau,rimage,Uii,Yii,P,N_small,PM,pairs,s,ver1);

%Transform with L matrix
if ver1
    L = L{1};
else
    start_ind = 0; 
    start_colloc = 0; 
end

for i = 1:P
    if ver1
        tau_x = tau_stokes_x((i-1)*N_small+1:i*N_small);
        tau_y = tau_stokes_y((i-1)*N_small+1:i*N_small);
        tau_stokes = [tau_x; tau_y]-L*[tau_x; tau_y]; 
    
        tau_stokes_x((i-1)*N_small+1:i*N_small) = tau_stokes(1:end/2);
        tau_stokes_y((i-1)*N_small+1:i*N_small) = tau_stokes(end/2+1:end);
    else
        tau_x = tau_stokes_x((i-1)*N_small+1:i*N_small);
        tau_y = tau_stokes_y((i-1)*N_small+1:i*N_small);
        inds = start_ind+1:start_ind+pairs(i,1);
        tau_st_x = tau_stress_x(inds);
        tau_st_y = tau_stress_y(inds);
        tau_p_x = tau_pot_x(inds);
        tau_p_y = tau_pot_y(inds);


        density = [tau_x; tau_st_x; tau_p_x; tau_y; tau_st_y; tau_p_y];

        projected = density-L{i}*density;
        tau_stokes_x((i-1)*N_small+1:i*N_small) = projected(1:N_small);
        tau_stokes_y((i-1)*N_small+1:i*N_small) = projected(N_small+2*pairs(i,1)+1:2*N_small+2*pairs(i,1));
        
        
        tau_stress_x(inds) = projected(N_small+1:N_small+pairs(i,1));
        tau_pot_x(inds) = projected(N_small+pairs(i,1)+1:N_small+2*pairs(i,1));
        
        tau_stress_y(inds) = projected(2*N_small+2*pairs(i,1)+1:2*N_small+3*pairs(i,1));
        tau_pot_y(inds) = projected(2*N_small+3*pairs(i,1)+1:end);
        
        start_ind = start_ind+pairs(i,1);
        start_colloc = start_colloc+pairs(i,2);
    end
end


if isempty(rimage)
    images = 0;
else
    images = 1; 
end


%Apply all stokeslets

%For debugging... 
%NS = singleLayer(rvec_in,rvec_out,mu);
%U_proxy = NS*[tau_stokes_x; tau_stokes_y]; %I think this vector has to be ordered differently


[ufmm,vfmm] = stokesSLPfmm(tau_stokes_x,tau_stokes_y,real(rvec_in),imag(rvec_in),real(rvec_check),imag(rvec_check),...
        0,5); %the last number here determines the tolerance. Higher number: higher tolerance. 
  
%NOTE! Errors might accumulate in the FMM if density norm is large!  

%[udirect,vdirect] = StokesletDirect(real(rvec_in),imag(rvec_in),real(rvec_out),imag(rvec_out),tau_stokes_x,tau_stokes_y,length(rvec_out));
% See the SLP FMM test in the fast tools folder


if ~images
    %no image points are used
    %res = res + U_proxy;
    res = [ufmm; vfmm];
    %res = res+[udirect'/4/pi; vdirect'/4/pi];
else    
    %want to do this with a matrix vector multiply instead!
    %Nim =  getImageKernels2D(rimage,nimage,rvec_check,mu,s); 
    u_stress = getStresslets(tau_stress_x,tau_stress_y,rimage,rvec_check,real(nimage),imag(nimage));
    %We assume a small number of images and do this with a dense matrix for
    %now. 
    %U_image = Nim*[tau_stress_x; tau_stress_y; tau_pot_x; tau_pot_y];
    u_pot = getPotdip(tau_pot_x, tau_pot_y,rimage,rvec_check);

    %for debugg purposes
    %Nim_temp =  getImageKernels2D(rimage,nimage,rvec_check,mu,[0 0 0 1]); 
    %U_image_temp = Nim_temp*[tau_pot_x; tau_pot_y]; 
    %res = [ufmm; vfmm]+U_image;
    res = [ufmm; vfmm]+u_stress+u_pot; 
    %res = res + [udirect'/4/pi; vdirect'/4/pi]+U_image;
end
    
if solve
    res = res + tau;
end

%Need to subract off the part computed twice for the one-body interactions
start_colloc = 0; 
start_im= 0; 

if ~images
    rout = rvec_out(start_colloc+1:start_colloc+pairs(1,2));
    rin = rvec_in(1:N_small);
    Nii = singleLayer(rin,rout,mu);
end

if solve

for i = 1:P
    % Get sources on this particle from single layer
    if images 
        rout = rvec_out(start_colloc+1:start_colloc+pairs(i,2));
        rin = rvec_in((i-1)*N_small+1:i*N_small);
        Nii = singleLayer(rin,rout,mu);
    end
    tau_xy = [tau_stokes_x((i-1)*N_small+1:i*N_small); tau_stokes_y((i-1)*N_small+1:i*N_small)];
    uii = Nii*tau_xy;

    if ver1 || (~images)
        uii = Nii*tau_xy;
    end

    %% Do the same thing with the images
    if images 
        rim = rimage(start_im+1:start_im+pairs(i,1)); 
        nim = nimage(start_im+1:start_im+pairs(i,1)); 
    
        Nim =  getImageKernels2D(rim,nim,rout,mu,s);
        stress_x = tau_stress_x(start_im+1:start_im+pairs(i,1));
        stress_y = tau_stress_y(start_im+1:start_im+pairs(i,1));
        pot_x = tau_pot_x(start_im+1:start_im+pairs(i,1));
        pot_y = tau_pot_y(start_im+1:start_im+pairs(i,1));
        
        if ver1
            tau_xy = [stress_x; stress_y; pot_x; pot_y];
            uii = uii+Nim*tau_xy;
        else
            Nimod = [Nii(:,1:end/2) Nim(:,1:pairs(i,1)) Nim(:,2*pairs(i,1)+1:3*pairs(i,1)) ... 
            Nii(:,end/2+1:end) Nim(:,pairs(i,1)+1:2*pairs(i,1)) Nim(:,3*pairs(i,1)+1:end)];
            tau_all = [tau_xy(1:end/2); stress_x; pot_x; tau_xy(end/2+1:end); stress_y; pot_y];
            uii = Nimod*tau_all; 
        end
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


end