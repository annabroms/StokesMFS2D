function opt = get2Dparams()
%get2Dparams returns struct with default parameters for a Stokes MFS 2D problem

N_c = 60; % sources on coarse proxy grid
opt.N_c = N_c; 

%set separation between true boundary and proxy surface
tol = 1e-12;
sep = (1/N_c)*log(1/tol); %separation between surfaces based on rule of thumb in Stein & Barnett QFS paper (2021) 
opt.Rp_c = max([1-sep,0.01]); %radius of proxy surface

opt.a_c = 1.2; % upsampling factor for collocation points so that M_c = a_c*N_c; 
opt.image = 1;%use images
opt.M_image = 35;  % Not always in use - sets extra collocation points for close to touching region

opt.s = [0 0 1 1]; %source types at clustered nodes
opt.alpha = 0; %use two lines?

opt.proj_all = 0; %only for the mobility problem

% Pair corrections
opt.pc = 0; %use pair corrections?
opt.delta_pair = 0.2; %activate pair_corrections
N_f = 150; %Set fine grid
opt.N_f = N_f;
tol = 1e-12;
sep = (1/N_f)*log(1/tol); %separation between surfaces based on rule of thumb in Stein & Barnett QFS paper (2021) 
opt.Rp_f = max([1-sep,0.01]);
opt.a_f = 1.2; % upsampling factor for collocation points so that M_f = a_f*N_f; 
 
opt.beta = 0.3; %determines distance between focus and tip of ellipse for enhancing nodes
opt.cmap = 0; %use coarse-to-coarse mapping. Only applicable with peanut compression. 
opt.gmres_verbose = 0; % 0=silent, 1=final summary, 2=per-iteration

end
