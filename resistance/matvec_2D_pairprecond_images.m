function res = matvec_2D_pairprecond_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rvec_out,q,U,Y,pairs,Upf,Ypf,s)


%Transform data at collocation points to data at source points. 
[rvec_in,rimage_in,nimage_in,coarse_ind,tau_stokes_x,tau_stokes_y,tau_stress_all_x,tau_stress_all_y,tau_pot_all_x,tau_pot_all_y,u_corr] = getPairTransformation(tau,rbase_in_c,rbase_in_f,refine,...
    rimage_vec,nimage,opt,rvec_out,q,U,Y,pairs,Upf,Ypf); 
 
P = length(q); 
N_large = size(rvec_out,1)/P;
mu = 1; 
PM = length(rvec_out);


%% Get velocity field from all source types 

rot = []; %for now, no rotlets are in play at image points

res = getVelocityField(rvec_in,rvec_out,tau_stokes_x,tau_stokes_y,rimage_in,nimage_in,rot,...
    tau_stress_all_x,tau_stress_all_y,tau_pot_all_x,tau_pot_all_y);

two_corr = 1; 


%Need to subract off self-interactions and replace with identity diagonals
rout = rvec_out(1:N_large)-q(1,:);
rin = rbase_in_c;
Nii = singleLayer(rin,rout,mu);

for i = 1:P
    % Get sources on this particle from single layer
    tau_xy = [tau_stokes_x(coarse_ind{i}); tau_stokes_y(coarse_ind{i})];
    
    %rout = rvec_out((i-1)*N_large+1:i*N_large,:);
    %rin = rbase_in_c+q(i);    
    %Nii = singleLayer(rin,rout,mu); %why recomputed?
    uii = Nii*tau_xy;

    %subract contribution in x
    res((i-1)*N_large+1:i*N_large) = res((i-1)*N_large+1:i*N_large)-uii(1:end/2);%+tau((i-1)*N_large+1:i*N_large);

    %subract contribution in y
    res((i-1)*N_large+1+PM:i*N_large+PM) = res((i-1)*N_large+1+PM:i*N_large+PM)-...
        uii(end/2+1:end);%+tau((i-1)*N_large+PM+1:i*N_large+PM);
   
end

if two_corr
    res = res-u_corr;
end 


res = res+tau;





end

