function opt = get2Dparams()
%get2Dparams returns struct with default parameters for a Stokes MFS 2D problem

N_c = 60; % sources on coarse proxy grid
opt.N_c = N_c; 

%set separation between true boundary and proxy surface
tol = 1e-12;
sep = (1/N_c)*log(1/tol); %separation between surfaces based on rule of thumb in Stein & Barnett QFS paper (2021) 
opt.Rp_c = max([1-sep,0.01]); %radius of proxy surface


opt.a_c = 1.2; % upsampling factor for collocation points so that M_c = a_c*N_c; 
opt.image = 1;
opt.M_image = 35;  % sets extra collocation points for close to touching region

%Set fine grid
N_f = 150;
opt.N_f = N_f;
tol = 1e-12;
sep = (1/N_f)*log(1/tol); %separation between surfaces based on rule of thumb in Stein & Barnett QFS paper (2021) 

opt.Rp_f = max([1-sep,0.01]);
opt.a_f = 1.2; % upsampling factor for collocation points so that M_f = a_f*N_f; 
opt.delta_pair = 0.2;


end